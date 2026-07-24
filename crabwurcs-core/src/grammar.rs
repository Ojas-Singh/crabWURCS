use crate::error::{CoreError, CoreResult};
use crate::model::{
    AnomericSymbol, CarbonPosition, Linkage, Modification, Monosaccharide, Probability,
    RepeatCount, ResidueGraph, RingClosure, UndefinedLinkage, UndefinedModification,
    UndefinedParent,
};
use nom::bytes::complete::{tag, take_until};
use nom::character::complete::char;
use nom::sequence::delimited;
use nom::IResult;
use petgraph::visit::EdgeRef;

type ParsedWurcsParts = (WURCSHeader, Vec<String>, Vec<u32>, String);
type ParsedEndpoint = (u32, Vec<u8>, Option<String>, Option<u8>);

pub fn parse_wurcs(input: &str) -> CoreResult<ResidueGraph> {
    let input = input.trim();
    let (_, (header, unique_residues, sequence, linkages)) =
        parse_wurcs_parts(input).map_err(|e| CoreError::ParseError {
            offset: 0,
            message: format!("{:?}", e),
        })?;

    let mut graph = assemble_graph(&header, &unique_residues, &sequence, &linkages)?;
    graph.set_source_wurcs(input.to_string());
    Ok(graph)
}

pub fn write_wurcs(graph: &ResidueGraph) -> CoreResult<String> {
    if let Some(name) = graph
        .inner()
        .node_weights()
        .find_map(|residue| residue.display_name.as_deref())
    {
        return Err(CoreError::UnrepresentableResidue(name.to_string()));
    }
    if let Some(kind) = graph.inner().node_weights().find_map(|residue| {
        residue
            .residue_kind
            .filter(|kind| kind.unique_residue().is_none())
    }) {
        return Err(CoreError::UnrepresentableResidue(
            kind.canonical_name().to_string(),
        ));
    }
    if let Some(source) = graph.source_wurcs() {
        return Ok(source.to_string());
    }
    let inner = graph.inner();
    let node_count = inner.node_count();
    let edge_count = inner.edge_count()
        + graph.undefined_linkages().len()
        + graph.undefined_modifications().len();

    let mut residue_map: std::collections::HashMap<String, usize> =
        std::collections::HashMap::new();
    let mut unique_residues: Vec<String> = Vec::new();
    let mut sequence: Vec<usize> = Vec::new();

    let mut node_to_seq_pos: std::collections::HashMap<usize, usize> =
        std::collections::HashMap::new();
    let mut ordered_nodes = Vec::with_capacity(node_count);
    let mut visited = std::collections::HashSet::new();
    let mut queue = std::collections::VecDeque::new();
    if let Some(root) = graph.root().or_else(|| inner.node_indices().next()) {
        queue.push_back(root);
    }
    while let Some(node) = queue.pop_front() {
        if !visited.insert(node) {
            continue;
        }
        ordered_nodes.push(node);
        let mut children: Vec<_> = inner
            .edges_directed(node, petgraph::Direction::Outgoing)
            .map(|edge| (edge.target(), edge.weight().parent_position.0))
            .collect();
        children.sort_by_key(|(_, position)| *position);
        for (child, _) in children {
            if !visited.contains(&child) {
                queue.push_back(child);
            }
        }
    }
    ordered_nodes.extend(inner.node_indices().filter(|node| !visited.contains(node)));

    for (seq_pos, ni) in ordered_nodes.into_iter().enumerate() {
        if let Some(residue) = inner.node_weight(ni) {
            let ur = write_unique_residue(residue);
            let idx = *residue_map.entry(ur.clone()).or_insert_with(|| {
                unique_residues.push(ur);
                unique_residues.len()
            });
            sequence.push(idx);
            node_to_seq_pos.insert(ni.index(), seq_pos);
        }
    }

    let unique_count = unique_residues.len();
    let seq_str = sequence
        .iter()
        .map(|i| i.to_string())
        .collect::<Vec<_>>()
        .join("-");

    let mut edge_parts = Vec::new();
    for edge_idx in inner.edge_indices() {
        if let Some((parent, child)) = inner.edge_endpoints(edge_idx) {
            if let Some(linkage) = inner.edge_weight(edge_idx) {
                // Use sequence position (0-based) instead of node index
                let parent_seq = *node_to_seq_pos.get(&parent.index()).unwrap_or(&0);
                let child_seq = *node_to_seq_pos.get(&child.index()).unwrap_or(&0);
                let parent_letter = (b'a' + parent_seq as u8) as char;
                let child_letter = (b'a' + child_seq as u8) as char;

                let mut parent_endpoint = linkage
                    .parent_positions()
                    .map(|position| {
                        format_endpoint_with_map(
                            parent_letter,
                            position,
                            linkage.parent_direction.as_deref(),
                            linkage.parent_modification_position,
                        )
                    })
                    .collect::<Vec<_>>()
                    .join("|");
                if let Some(probability) = linkage.parent_probability {
                    parent_endpoint.push_str(&format!("%{}%", probability.to_wurcs()));
                }
                let mut child_endpoint = linkage
                    .child_positions()
                    .map(|position| {
                        format_endpoint_with_map(
                            child_letter,
                            position,
                            linkage.child_direction.as_deref(),
                            linkage.child_modification_position,
                        )
                    })
                    .collect::<Vec<_>>()
                    .join("|");
                if let Some(probability) = linkage.child_probability {
                    child_endpoint.push_str(&format!("%{}%", probability.to_wurcs()));
                }
                let repeat_suffix = linkage
                    .repeat
                    .as_ref()
                    .map(|repeat| format!("~{}", repeat.to_wurcs()))
                    .unwrap_or_default();
                edge_parts.push(format!(
                    "{}-{}{}{}",
                    parent_endpoint,
                    child_endpoint,
                    linkage.map_code.as_deref().unwrap_or_default(),
                    repeat_suffix
                ));
            }
        }
    }

    for undefined in graph.undefined_linkages() {
        let Some(child_seq) = node_to_seq_pos.get(&undefined.child.index()) else {
            continue;
        };
        let child_letter = (b'a' + *child_seq as u8) as char;
        let child_endpoint = undefined
            .child_positions
            .iter()
            .copied()
            .map(|position| format_endpoint(child_letter, position))
            .collect::<Vec<_>>()
            .join("|");
        let mut parent_endpoints = Vec::new();
        for parent in &undefined.parents {
            let Some(parent_seq) = node_to_seq_pos.get(&parent.residue.index()) else {
                continue;
            };
            let parent_letter = (b'a' + *parent_seq as u8) as char;
            parent_endpoints.extend(
                parent
                    .positions
                    .iter()
                    .copied()
                    .map(|position| format_endpoint(parent_letter, position)),
            );
        }
        edge_parts.push(format!(
            "{}-{}{}",
            child_endpoint,
            parent_endpoints.join("|"),
            '}'
        ));
    }

    for undefined in graph.undefined_modifications() {
        let mut parent_endpoints = Vec::new();
        for parent in &undefined.parents {
            let Some(parent_seq) = node_to_seq_pos.get(&parent.residue.index()) else {
                continue;
            };
            let parent_letter = (b'a' + *parent_seq as u8) as char;
            parent_endpoints.extend(
                parent
                    .positions
                    .iter()
                    .copied()
                    .map(|position| format_endpoint(parent_letter, position)),
            );
        }
        edge_parts.push(format!(
            "{}}}{}",
            parent_endpoints.join("|"),
            undefined.map_code
        ));
    }

    edge_parts.sort_by_key(|part| {
        part.chars()
            .skip(1)
            .take_while(char::is_ascii_digit)
            .collect::<String>()
            .parse::<u8>()
            .unwrap_or(u8::MAX)
    });

    let linkage_str = edge_parts.join("_");

    let edge_count_field = if graph.is_composition() {
        format!("{}+", edge_count)
    } else {
        edge_count.to_string()
    };
    let header = format!(
        "WURCS=2.0/{},{},{}",
        unique_count, node_count, edge_count_field
    );
    let ur_str = unique_residues
        .iter()
        .map(|ur| format!("[{}]", ur))
        .collect::<Vec<_>>()
        .join("");

    Ok(format!("{}/{}/{}/{}", header, ur_str, seq_str, linkage_str))
}

fn format_endpoint(residue: char, position: CarbonPosition) -> String {
    if position.0 == 0 {
        format!("{residue}?")
    } else {
        format!("{residue}{}", position.0)
    }
}

fn format_endpoint_with_map(
    residue: char,
    position: CarbonPosition,
    direction: Option<&str>,
    modification_position: Option<u8>,
) -> String {
    let mut endpoint = format_endpoint(residue, position);
    if let Some(direction) = direction {
        endpoint.push_str(direction);
    }
    if let Some(position) = modification_position {
        endpoint.push_str(&position.to_string());
    }
    endpoint
}

fn write_unique_residue(residue: &Monosaccharide) -> String {
    let anomeric_char = residue.anomeric_prefix.as_str();

    // The terminal h/m/x is part of the carbon backbone, not the ring size.
    // Remove the single terminal-carbon descriptor. `trim_end_matches` would
    // also erase an unknown `x` stereochemical backbone such as `xxxxm`.
    let skeleton_base = residue
        .skeleton_code
        .strip_suffix(['h', 'm', 'a', 'x'])
        .unwrap_or(&residue.skeleton_code);

    let terminal = residue
        .skeleton_code
        .chars()
        .last()
        .filter(|character| matches!(character, 'h' | 'm' | 'a' | 'x'))
        .unwrap_or(match residue.ring {
            RingClosure::Pyranose => 'h',
            RingClosure::Furanose => 'm',
            RingClosure::Open | RingClosure::Unknown => 'x',
        });

    // For unknown anomeric (reducing end), don't write the anomeric position
    // Uronic-acid skeletons encode the terminal carbon as `A` and do not use
    // the aldose `h`/`m` suffix in WURCS unique-residue notation.
    let ring_suffix = if skeleton_base.ends_with('A') {
        String::new()
    } else {
        terminal.to_string()
    };
    let header_part = if residue.anomeric_position > 0 {
        format!(
            "{}{}{}-{}{}",
            anomeric_char,
            skeleton_base,
            ring_suffix,
            residue.anomeric_position,
            residue.anomeric_symbol.to_char()
        )
    } else {
        format!("{}{}{}", anomeric_char, skeleton_base, ring_suffix)
    };

    let mut segments = vec![header_part];

    // Only write ring position for non-reducing end residues
    if residue.anomeric_position > 0 && (residue.ring_start.is_some() || residue.ring_end.is_some())
    {
        let ring_start = residue
            .ring_start
            .map(|p| p.to_string())
            .unwrap_or_else(|| "?".to_string());
        let ring_end = residue
            .ring_end
            .map(|p| p.to_string())
            .unwrap_or_else(|| "?".to_string());
        segments.push(format!("{}-{}", ring_start, ring_end));
    }

    for m in &residue.modifications {
        let position = if m.position.0 == 0 {
            "?".to_string()
        } else {
            m.position.0.to_string()
        };
        let probability = m
            .probability
            .map(|value| format!("%{}%", value.to_wurcs()))
            .unwrap_or_default();
        segments.push(format!("{}{}*{}", position, probability, m.descriptor));
    }

    segments.join("_")
}

pub fn standardize_wurcs(input: &str) -> CoreResult<String> {
    let graph = parse_wurcs(input)?;
    write_wurcs(&graph)
}

#[derive(Debug)]
struct WURCSHeader {
    #[allow(dead_code)]
    unique_residue_count: u32,
    #[allow(dead_code)]
    total_residue_count: u32,
    #[allow(dead_code)]
    linkage_count: u32,
    composition: bool,
}

fn parse_wurcs_parts(input: &str) -> IResult<&str, ParsedWurcsParts> {
    let (input, _) = tag("WURCS=")(input)?;
    let (input, _) = tag("2.0/")(input)?;
    let (input, counts_str) = take_until("/")(input)?;
    let (input, _) = tag("/")(input)?;

    let parts: Vec<&str> = counts_str.split(',').collect();
    if parts.len() < 3 {
        return Err(nom::Err::Failure(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Fail,
        )));
    }

    let unique_count: u32 = parts[0].parse().map_err(|_| {
        nom::Err::Failure(nom::error::Error::new(input, nom::error::ErrorKind::Fail))
    })?;
    let total_count: u32 = parts[1].parse().map_err(|_| {
        nom::Err::Failure(nom::error::Error::new(input, nom::error::ErrorKind::Fail))
    })?;

    let third_part = parts[2];
    let linkage_count: u32 = if third_part.contains('+') {
        third_part
            .split('+')
            .next()
            .unwrap_or("0")
            .parse()
            .map_err(|_| {
                nom::Err::Failure(nom::error::Error::new(input, nom::error::ErrorKind::Fail))
            })?
    } else {
        third_part.parse().map_err(|_| {
            nom::Err::Failure(nom::error::Error::new(input, nom::error::ErrorKind::Fail))
        })?
    };

    let header = WURCSHeader {
        unique_residue_count: unique_count,
        total_residue_count: total_count,
        linkage_count,
        composition: third_part.ends_with('+'),
    };

    let (input, unique_residues) = parse_bracketed_residues(input, unique_count)?;
    let (input, _) = tag("/")(input)?;

    let (input, sequence_str) = take_until("/")(input)?;
    let (input, _) = tag("/")(input)?;

    let sequence: Vec<u32> = sequence_str
        .split('-')
        .filter(|s| !s.is_empty())
        .map(|s| s.parse::<u32>().unwrap_or(0))
        .collect();

    let linkage_str = input;

    Ok((
        "",
        (header, unique_residues, sequence, linkage_str.to_string()),
    ))
}

fn parse_bracketed_residues(input: &str, count: u32) -> IResult<&str, Vec<String>> {
    let mut residues = Vec::new();
    let mut remaining = input;
    let mut parsed = 0u32;

    while parsed < count && !remaining.is_empty() {
        if !remaining.starts_with('[') {
            break;
        }
        let (rest, content) = delimited(char('['), take_until("]"), char(']'))(remaining)?;
        residues.push(content.to_string());
        remaining = rest;
        parsed += 1;
    }

    Ok((remaining, residues))
}

fn parse_single_residue(input: &str) -> Monosaccharide {
    let parts: Vec<&str> = input.split('_').collect();
    if parts.is_empty() {
        return Monosaccharide::new(
            0,
            String::new(),
            vec![],
            RingClosure::Open,
            None,
            None,
            0,
            AnomericSymbol::Unknown,
            String::from("x"),
            vec![],
        );
    }

    let header_part = parts[0];
    let (anomeric_prefix, rest) = parse_anomeric_prefix(header_part);
    let (backbone, ring_char, rest_after_ring) = parse_backbone_and_ring(rest);

    let (anomeric_pos, anomeric_sym) =
        if let Some(rest_after_dash) = rest_after_ring.strip_prefix('-') {
            let (pos, sym, _rest) = parse_anomeric_info(rest_after_dash);
            (pos, sym)
        } else {
            (0, AnomericSymbol::Unknown)
        };

    let mut ring_start: Option<u8> = None;
    let mut ring_end: Option<u8> = None;
    let mut mod_start_idx = 1usize;

    if parts.len() > 1 {
        let second = parts[1];
        if !second.contains('*')
            && (second.contains('-')
                || second == "?"
                || second
                    .chars()
                    .all(|c| c.is_ascii_digit() || c == '-' || c == '?'))
        {
            if let Some((start_str, end_str)) = second.split_once('-') {
                ring_start = start_str.parse::<u8>().ok();
                ring_end = end_str.parse::<u8>().ok();
            }
            mod_start_idx = 2;
        }
    }

    let ring = if let (Some(start), Some(end)) = (ring_start, ring_end) {
        match end.saturating_sub(start) {
            4 => RingClosure::Pyranose,
            3 => RingClosure::Furanose,
            _ => RingClosure::Open,
        }
    } else if anomeric_sym == AnomericSymbol::OpenChain || ring_char == 'x' {
        RingClosure::Open
    } else {
        RingClosure::Unknown
    };

    let mut modifications = Vec::new();
    for part in parts.iter().skip(mod_start_idx) {
        if let Some((location, desc)) = part.split_once('*') {
            if let Some((pos, probability)) = parse_modification_location(location) {
                modifications.push(Modification {
                    position: CarbonPosition(pos),
                    descriptor: desc.to_string(),
                    probability,
                });
            }
        }
    }

    // Keep the terminal carbon descriptor; the anomeric prefix remains a
    // separate field and ring closure is represented by ring_start/end.
    let skeleton_code = if ring_char == '\0' {
        backbone
    } else {
        format!("{backbone}{ring_char}")
    };
    let backbone_length = skeleton_code.chars().filter(|c| c.is_ascii_digit()).count() as u8;

    Monosaccharide {
        backbone_length,
        skeleton_code,
        stereo: vec![],
        ring,
        ring_start,
        ring_end,
        anomeric_position: anomeric_pos,
        anomeric_symbol: anomeric_sym,
        anomeric_prefix,
        modifications,
        display_name: None,
        residue_kind: None,
    }
}

fn parse_modification_location(value: &str) -> Option<(u8, Option<Probability>)> {
    let (position, probability) = if let Some((position, probability)) = value.split_once('%') {
        (
            position,
            Some(Probability::parse_wurcs(probability.strip_suffix('%')?)?),
        )
    } else {
        (value, None)
    };
    let position = if position == "?" {
        0
    } else {
        position.parse().ok()?
    };
    Some((position, probability))
}

fn parse_anomeric_prefix(input: &str) -> (String, &str) {
    if input.is_empty() {
        return (String::from("x"), input);
    }
    // Prefixes have a small grammar; skeletons may themselves begin with
    // `x` or `d`, so consuming an arbitrary letter run erases generic and
    // dideoxy backbones such as `axxxxh` and `adxxxm`.
    for prefix in ["Aad", "AUd", "AOd", "Ad", "ha", "hU", "hO"] {
        if let Some(rest) = input.strip_prefix(prefix) {
            return (prefix.to_string(), rest);
        }
    }
    if let Some(prefix) = input.chars().next().filter(|character| {
        matches!(character, 'u' | 'a' | 'o' | 'U' | 'A')
    }) {
        let end = prefix.len_utf8();
        (input[..end].to_string(), &input[end..])
    } else {
        (String::from("x"), input)
    }
}

fn parse_backbone_and_ring(input: &str) -> (String, char, &str) {
    if input.is_empty() {
        return (String::new(), 'x', input);
    }

    let suffix_start = input.find('-').unwrap_or(input.len());
    let head = &input[..suffix_start];
    let rest = &input[suffix_start..];
    if let Some(ring_char @ ('h' | 'm' | 'x')) = head.chars().last() {
        let backbone_end = head.len() - ring_char.len_utf8();
        (head[..backbone_end].to_string(), ring_char, rest)
    } else {
        (head.to_string(), '\0', rest)
    }
}

fn parse_anomeric_info(input: &str) -> (u8, AnomericSymbol, &str) {
    if input.len() < 2 {
        return (0, AnomericSymbol::Unknown, input);
    }

    let pos = input
        .chars()
        .next()
        .and_then(|c| c.to_digit(10))
        .unwrap_or(0) as u8;
    let sym_char = input.chars().nth(1).unwrap_or('x');
    let rest = &input[2..];

    let sym = match sym_char {
        'a' => AnomericSymbol::Alpha,
        'b' => AnomericSymbol::Beta,
        'o' => AnomericSymbol::OpenChain,
        _ => AnomericSymbol::Unknown,
    };

    (pos, sym, rest)
}

fn assemble_graph(
    header: &WURCSHeader,
    unique_residues: &[String],
    sequence: &[u32],
    linkage_str: &str,
) -> CoreResult<ResidueGraph> {
    let mut graph = ResidueGraph::new();
    if header.composition {
        graph.set_composition(true);
    }

    let mut parsed_residues: Vec<Monosaccharide> = Vec::new();
    for ur_str in unique_residues {
        let mono = parse_single_residue(ur_str);
        parsed_residues.push(mono);
    }

    // A reducing-end unique residue may omit `_1-4`/`_1-5` (for example
    // `[u211h]`) even when another anomer of the same backbone declares
    // the ring closure. Reuse that unambiguous evidence instead of
    // guessing from the terminal carbon descriptor.
    for index in 0..parsed_residues.len() {
        if parsed_residues[index].ring_start.is_some() {
            continue;
        }
        let skeleton = parsed_residues[index].skeleton_code.clone();
        let declared = parsed_residues
            .iter()
            .filter(|candidate| {
                candidate.skeleton_code == skeleton && candidate.ring_start.is_some()
            })
            .map(|candidate| candidate.ring)
            .collect::<Vec<_>>();
        if let Some(first) = declared.first().copied() {
            if declared.iter().all(|ring| *ring == first) {
                parsed_residues[index].ring = first;
            }
        }
    }

    let mut node_map: Vec<petgraph::graph::NodeIndex> = Vec::new();
    for seq_idx in sequence {
        let residue_idx = (*seq_idx as usize).wrapping_sub(1);
        if residue_idx < parsed_residues.len() {
            let mono = parsed_residues[residue_idx].clone();
            let node_idx = graph.add_residue(mono);
            node_map.push(node_idx);
        }
    }

    let linkages = parse_linkages(linkage_str);

    for parsed_linkage in linkages {
        match parsed_linkage {
            ParsedLinkage::Defined {
                from_seq,
                from_positions,
                to_seq,
                to_positions,
                repeat,
                from_probability,
                to_probability,
                map_code,
                from_direction,
                from_modification_position,
                to_direction,
                to_modification_position,
            } => {
                let from_idx = (from_seq as usize).wrapping_sub(1);
                let to_idx = (to_seq as usize).wrapping_sub(1);
                if from_idx >= node_map.len() || to_idx >= node_map.len() {
                    continue;
                }

                let from_anomer = graph
                    .residue(node_map[from_idx])
                    .map(|residue| residue.anomeric_position)
                    .unwrap_or(0);
                let to_anomer = graph
                    .residue(node_map[to_idx])
                    .map(|residue| residue.anomeric_position)
                    .unwrap_or(0);
                let from_is_donor = from_anomer > 0 && from_positions.contains(&from_anomer);
                let to_is_donor = to_anomer > 0 && to_positions.contains(&to_anomer);
                let (
                    parent,
                    child,
                    parent_positions,
                    child_positions,
                    parent_probability,
                    child_probability,
                    parent_direction,
                    parent_modification_position,
                    child_direction,
                    child_modification_position,
                ) = if from_is_donor && !to_is_donor {
                    (
                        node_map[to_idx],
                        node_map[from_idx],
                        to_positions,
                        from_positions,
                        to_probability,
                        from_probability,
                        to_direction,
                        to_modification_position,
                        from_direction,
                        from_modification_position,
                    )
                } else {
                    (
                        node_map[from_idx],
                        node_map[to_idx],
                        from_positions,
                        to_positions,
                        from_probability,
                        to_probability,
                        from_direction,
                        from_modification_position,
                        to_direction,
                        to_modification_position,
                    )
                };
                let mut linkage = Linkage::with_alternatives(
                    parent_positions.into_iter().map(CarbonPosition).collect(),
                    child_positions.into_iter().map(CarbonPosition).collect(),
                );
                linkage.repeat = repeat;
                linkage.parent_probability = parent_probability;
                linkage.child_probability = child_probability;
                linkage.map_code = map_code;
                linkage.parent_direction = parent_direction;
                linkage.parent_modification_position = parent_modification_position;
                linkage.child_direction = child_direction;
                linkage.child_modification_position = child_modification_position;
                linkage.cyclic =
                    linkage.repeat.is_none() && node_map.len() > 1 && child == node_map[0];
                graph.add_linkage(parent, child, linkage);
            }
            ParsedLinkage::Undefined {
                child_seq,
                child_positions,
                parent_candidates,
            } => {
                let child_idx = (child_seq as usize).wrapping_sub(1);
                if child_idx >= node_map.len() {
                    continue;
                }
                let parents = parent_candidates
                    .into_iter()
                    .filter_map(|(sequence, positions)| {
                        let index = (sequence as usize).wrapping_sub(1);
                        node_map.get(index).copied().map(|residue| UndefinedParent {
                            residue,
                            positions: positions.into_iter().map(CarbonPosition).collect(),
                        })
                    })
                    .collect();
                graph.add_undefined_linkage(UndefinedLinkage {
                    child: node_map[child_idx],
                    child_positions: child_positions.into_iter().map(CarbonPosition).collect(),
                    parents,
                });
            }
            ParsedLinkage::UndefinedModification {
                parent_candidates,
                map_code,
            } => {
                let parents = parent_candidates
                    .into_iter()
                    .filter_map(|(sequence, positions)| {
                        let index = (sequence as usize).wrapping_sub(1);
                        node_map.get(index).copied().map(|residue| UndefinedParent {
                            residue,
                            positions: positions.into_iter().map(CarbonPosition).collect(),
                        })
                    })
                    .collect();
                graph.add_undefined_modification(UndefinedModification { parents, map_code });
            }
        }
    }

    if !node_map.is_empty() {
        graph.set_root(node_map[0]);
    }

    Ok(graph)
}

#[derive(Debug)]
enum ParsedLinkage {
    Defined {
        from_seq: u32,
        from_positions: Vec<u8>,
        to_seq: u32,
        to_positions: Vec<u8>,
        repeat: Option<RepeatCount>,
        from_probability: Option<Probability>,
        to_probability: Option<Probability>,
        map_code: Option<String>,
        from_direction: Option<String>,
        from_modification_position: Option<u8>,
        to_direction: Option<String>,
        to_modification_position: Option<u8>,
    },
    Undefined {
        child_seq: u32,
        child_positions: Vec<u8>,
        parent_candidates: Vec<(u32, Vec<u8>)>,
    },
    UndefinedModification {
        parent_candidates: Vec<(u32, Vec<u8>)>,
        map_code: String,
    },
}

fn parse_linkages(linkage_str: &str) -> Vec<ParsedLinkage> {
    let mut result = Vec::new();
    if linkage_str.is_empty() {
        return result;
    }

    for segment in linkage_str.split('_') {
        if let Some(undefined) = parse_undefined_linkage(segment) {
            result.push(undefined);
            continue;
        }
        let (segment, repeat) = match segment.rsplit_once('~') {
            Some((linkage, count)) => (linkage, RepeatCount::parse(count)),
            None => (segment, None),
        };
        let (segment, map_code) = match segment.find('*') {
            Some(index) => (&segment[..index], Some(segment[index..].to_string())),
            None => (segment, None),
        };
        if segment.contains('-') {
            if let Some((from_part, to_part)) = segment.split_once('-') {
                let (from_part, from_probability) = strip_endpoint_probability(from_part);
                let (to_part, to_probability) = strip_endpoint_probability(to_part);
                let from = parse_endpoint(from_part);
                let to = parse_endpoint(to_part);
                if let (
                    Some((fr, fp, from_direction, from_modification_position)),
                    Some((tr, tp, to_direction, to_modification_position)),
                ) = (from, to)
                {
                    result.push(ParsedLinkage::Defined {
                        from_seq: fr,
                        from_positions: fp,
                        to_seq: tr,
                        to_positions: tp,
                        repeat,
                        from_probability,
                        to_probability,
                        map_code,
                        from_direction,
                        from_modification_position,
                        to_direction,
                        to_modification_position,
                    });
                }
            }
        }
    }

    result
}

fn strip_endpoint_probability(endpoint: &str) -> (&str, Option<Probability>) {
    let Some(start) = endpoint.find('%') else {
        return (endpoint, None);
    };
    let base = &endpoint[..start];
    let probability = endpoint[start + 1..]
        .strip_suffix('%')
        .and_then(Probability::parse_wurcs);
    (base, probability)
}

fn parse_undefined_linkage(segment: &str) -> Option<ParsedLinkage> {
    if let Some((parents, map)) = segment.split_once("}*") {
        let parent_candidates = parse_endpoint_groups(parents);
        if parent_candidates.is_empty() {
            return None;
        }
        return Some(ParsedLinkage::UndefinedModification {
            parent_candidates,
            map_code: format!("*{map}"),
        });
    }
    let segment = segment.strip_suffix('}')?;
    let (child, parents) = segment.split_once('-')?;
    let mut child_groups = parse_endpoint_groups(child);
    let parent_groups = parse_endpoint_groups(parents);
    if child_groups.len() != 1 || parent_groups.is_empty() {
        return None;
    }
    let (child_seq, child_positions) = child_groups.pop()?;
    Some(ParsedLinkage::Undefined {
        child_seq,
        child_positions,
        parent_candidates: parent_groups,
    })
}

fn parse_endpoint_groups(endpoint: &str) -> Vec<(u32, Vec<u8>)> {
    let mut groups: Vec<(u32, Vec<u8>)> = Vec::new();
    for alternative in endpoint.split('|') {
        let Some(sequence) = extract_seq_index(alternative) else {
            continue;
        };
        let position = extract_last_number(alternative).unwrap_or(0);
        if let Some((_, positions)) = groups.iter_mut().find(|(value, _)| *value == sequence) {
            positions.push(position);
        } else {
            groups.push((sequence, vec![position]));
        }
    }
    groups
}

fn parse_endpoint(endpoint: &str) -> Option<ParsedEndpoint> {
    let mut residue = None;
    let mut positions = Vec::new();
    let mut direction = None;
    let mut modification_position = None;
    for alternative in endpoint.split('|') {
        residue = residue.or_else(|| extract_seq_index(alternative));
        let remainder = &alternative[alternative.chars().next()?.len_utf8()..];
        let position_end = remainder
            .find(|character: char| !character.is_ascii_digit() && character != '?')
            .unwrap_or(remainder.len());
        let position = match &remainder[..position_end] {
            "" | "?" => 0,
            value => value.parse().ok()?,
        };
        positions.push(position);
        let annotation = &remainder[position_end..];
        if !annotation.is_empty() {
            let direction_end = annotation
                .find(char::is_numeric)
                .unwrap_or(annotation.len());
            if direction_end > 0 {
                direction = Some(annotation[..direction_end].to_string());
            }
            if direction_end < annotation.len() {
                modification_position = annotation[direction_end..].parse().ok();
            }
        }
    }
    Some((residue?, positions, direction, modification_position))
}

fn extract_last_number(s: &str) -> Option<u8> {
    let mut digits = String::new();
    for c in s.chars().rev() {
        if c.is_ascii_digit() {
            digits.push(c);
        } else {
            break;
        }
    }
    digits.chars().rev().collect::<String>().parse().ok()
}

fn extract_seq_index(s: &str) -> Option<u32> {
    let first_char = s.chars().next()?;
    if first_char.is_ascii_lowercase() {
        Some((first_char as u32).wrapping_sub('a' as u32) + 1)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_simple_monosaccharide() {
        let result = parse_wurcs("WURCS=2.0/1,1,0/[a2122h-1x_1-5_6*OPO/3O/3=O]/1/");
        assert!(result.is_ok(), "Failed: {:?}", result.err());
        let graph = result.unwrap();
        assert_eq!(graph.node_count(), 1);
        assert_eq!(graph.edge_count(), 0);
        let root = graph.root().unwrap();
        let residue = graph.residue(root).unwrap();
        assert_eq!(residue.backbone_length, 4);
        assert_eq!(residue.anomeric_symbol, AnomericSymbol::Unknown);
        assert_eq!(residue.anomeric_position, 1);
        assert_eq!(residue.modifications.len(), 1);
        assert_eq!(residue.modifications[0].position.0, 6);
    }

    #[test]
    fn test_parse_simple_disaccharide() {
        let result = parse_wurcs(
            "WURCS=2.0/4,6,5/[a2112h-1b_1-?_2*NCC/3=O][a2112h-1b_1-5][a1221m-1a_1-5][a2122h-1b_1-5_2*NCC/3=O]/1-2-3-4-2-3/a3-b1_a6-d1_b2-c1_d4-e1_e2-f1",
        );
        assert!(result.is_ok(), "Failed: {:?}", result.err());
        let graph = result.unwrap();
        assert_eq!(graph.node_count(), 6);
        assert_eq!(graph.edge_count(), 5);
    }

    #[test]
    fn test_parse_monosaccharide_without_mods() {
        let result = parse_wurcs("WURCS=2.0/1,1,0/[a2122h-1b_1-5]/1/");
        assert!(result.is_ok(), "Failed: {:?}", result.err());
        let graph = result.unwrap();
        assert_eq!(graph.node_count(), 1);
    }

    #[test]
    fn test_parse_linear_chain() {
        let result = parse_wurcs(
            "WURCS=2.0/4,4,3/[a2122h-1b_1-5][a2112h-1b_1-5][a2122h-1a_1-5_2*NCC/3=O][Aad21122h-2a_2-6_5*NCC/3=O]/1-2-3-4/a4-b1_b4-c1_c4-d1",
        );
        assert!(result.is_ok(), "Failed: {:?}", result.err());
        let graph = result.unwrap();
        assert_eq!(graph.node_count(), 4);
        assert_eq!(graph.edge_count(), 3);
    }

    #[test]
    fn test_roundtrip_monosaccharide() {
        let input = "WURCS=2.0/1,1,0/[a2122h-1b_1-5]/1/";
        let graph = parse_wurcs(input).unwrap();
        let output = write_wurcs(&graph).unwrap();
        assert!(output.contains("WURCS=2.0/"));
    }

    #[test]
    fn test_parse_cyclic() {
        let result = parse_wurcs(
            "WURCS=2.0/1,7,7/[a2122h-1a_1-5_2*OC_3*OC_6*N]/1-1-1-1-1-1-1/a1-g4_a4-b1_b4-c1_c4-d1_d4-e1_e4-f1_f4-g1",
        );
        assert!(result.is_ok(), "Failed: {:?}", result.err());
        let graph = result.unwrap();
        assert_eq!(graph.node_count(), 7);
        assert_eq!(graph.edge_count(), 7);
    }

    #[test]
    fn test_parse_fuzzy_linkage() {
        let result = parse_wurcs(
            "WURCS=2.0/4,5,4/[h2112h_2*NCC/3=O][a2122h-1b_1-5_2*NCC/3=O][a2112h-1b_1-5][a1221m-1a_1-5]/1-2-3-3-4/b4-c1_d2-e1_b1-a3|a6_d1-a3|a6",
        );
        assert!(result.is_ok(), "Failed: {:?}", result.err());
    }

    #[test]
    fn test_parse_substituent() {
        let result = parse_wurcs(
            "WURCS=2.0/3,4,3/[u2122h_2*NCC/3=O_6*OSO/3=O/3=O][a2112h-1b_1-5][a2122h-1b_1-5_2*NCC/3=O_6*OSO/3=O/3=O]/1-2-3-2/a4-b1_b3-c1_c4-d1",
        );
        assert!(result.is_ok(), "Failed: {:?}", result.err());
        let graph = result.unwrap();
        let root = graph.root().unwrap();
        let residue = graph.residue(root).unwrap();
        assert_eq!(residue.modifications.len(), 2);
    }

    #[test]
    fn test_probability_annotations_survive_graph_edits() {
        let mut modification_graph =
            parse_wurcs("WURCS=2.0/1,1,0/[a2112A-1a_1-5_6%?%*OC]/1/").unwrap();
        let modification = &modification_graph
            .residue(modification_graph.root().unwrap())
            .unwrap()
            .modifications[0];
        assert_eq!(
            modification.probability,
            Some(Probability {
                lower: crate::ProbabilityValue::Unknown,
                upper: crate::ProbabilityValue::Unknown,
            })
        );
        let _ = modification_graph.inner_mut();
        let generated = write_wurcs(&modification_graph).unwrap();
        assert!(generated.contains("6%?%*OC"), "{generated}");

        let mut linkage_graph =
            parse_wurcs("WURCS=2.0/1,2,1/[a2122h-1b_1-5]/1-1/b1-a4%.55%").unwrap();
        let linkage_probability = linkage_graph
            .inner()
            .edge_weights()
            .next()
            .unwrap()
            .parent_probability;
        assert_eq!(
            linkage_probability,
            Some(Probability {
                lower: crate::ProbabilityValue::Known(5500),
                upper: crate::ProbabilityValue::Known(5500),
            })
        );
        let _ = linkage_graph.inner_mut();
        let generated = write_wurcs(&linkage_graph).unwrap();
        assert!(generated.contains("%.55%"), "{generated}");
        let reparsed = parse_wurcs(&generated).unwrap();
        assert_eq!(
            reparsed
                .inner()
                .edge_weights()
                .next()
                .unwrap()
                .parent_probability,
            linkage_probability
        );
    }

    #[test]
    fn test_map_bridge_survives_graph_edits() {
        let input = "WURCS=2.0/3,3,2/[a261m-1a_1-4_3*C=O][a1211h-1a_1-5_2*NC][a222h-1b_1-4_1*N]/1-2-3/a2-b1_b3-c5*OPO*/3O/3=O";
        let mut graph = parse_wurcs(input).unwrap();
        let bridge = graph
            .inner()
            .edge_weights()
            .find(|linkage| linkage.map_code.is_some())
            .unwrap();
        assert_eq!(bridge.map_code.as_deref(), Some("*OPO*/3O/3=O"));
        assert_eq!(bridge.parent_position, CarbonPosition(3));
        assert_eq!(bridge.child_position, CarbonPosition(5));

        let _ = graph.inner_mut();
        let generated = write_wurcs(&graph).unwrap();
        assert!(generated.contains("*OPO*/3O/3=O"), "{generated}");
        let reparsed = parse_wurcs(&generated).unwrap();
        assert!(reparsed
            .inner()
            .edge_weights()
            .any(|linkage| linkage.map_code.as_deref() == Some("*OPO*/3O/3=O")));
    }

    #[test]
    fn test_map_endpoint_directions_survive_graph_edits() {
        let input = "WURCS=2.0/2,2,1/[hxh][a2122h-1b_1-5]/1-2/a3n2-b1n1*1NCCOP^XO*2/6O/6=O";
        let mut graph = parse_wurcs(input).unwrap();
        let bridge = graph.inner().edge_weights().next().unwrap();
        assert_eq!(bridge.parent_direction.as_deref(), Some("n"));
        assert_eq!(bridge.parent_modification_position, Some(2));
        assert_eq!(bridge.child_direction.as_deref(), Some("n"));
        assert_eq!(bridge.child_modification_position, Some(1));
        assert_eq!(bridge.map_code.as_deref(), Some("*1NCCOP^XO*2/6O/6=O"));

        let _ = graph.inner_mut();
        let generated = write_wurcs(&graph).unwrap();
        assert!(
            generated.contains("a3n2-b1n1*1NCCOP^XO*2/6O/6=O"),
            "{generated}"
        );
    }

    #[test]
    fn test_undefined_modification_survives_graph_edits() {
        let input = "WURCS=2.0/2,2,1/[u2122h][u2112h]/1-2/a?|b?}*OCC/3=O";
        let mut graph = parse_wurcs(input).unwrap();
        assert_eq!(graph.undefined_modifications().len(), 1);
        let modification = &graph.undefined_modifications()[0];
        assert_eq!(modification.parents.len(), 2);
        assert_eq!(modification.map_code, "*OCC/3=O");

        let _ = graph.inner_mut();
        let generated = write_wurcs(&graph).unwrap();
        assert!(generated.contains("a?|b?}*OCC/3=O"), "{generated}");
        let reparsed = parse_wurcs(&generated).unwrap();
        assert_eq!(
            reparsed.undefined_modifications(),
            graph.undefined_modifications()
        );
    }

    #[test]
    fn test_parse_repeat() {
        let mut graph = parse_wurcs(
            "WURCS=2.0/4,5,5/[a2112h-1b_1-5_2*NCC/3=O][a2122A-1b_1-5][a2112h-1a_1-5][a2122h-1b_1-5]/1-2-1-3-4/a3-b1_b4-c1_c4-d1_d3-e1_a1-e4~n",
        )
        .unwrap();
        let repeat = graph
            .inner()
            .edge_weights()
            .find_map(|linkage| linkage.repeat.as_ref());
        assert_eq!(repeat, Some(&RepeatCount::Unknown));

        // Invalidate verbatim provenance and prove repeat semantics survive
        // editable graph serialization.
        let _ = graph.inner_mut();
        let generated = write_wurcs(&graph).unwrap();
        assert!(generated.contains("~n"), "{generated}");
        let reparsed = parse_wurcs(&generated).unwrap();
        assert!(reparsed
            .inner()
            .edge_weights()
            .any(|linkage| linkage.repeat == Some(RepeatCount::Unknown)));
    }

    #[test]
    fn test_repeat_count_ranges() {
        assert_eq!(RepeatCount::parse("3"), Some(RepeatCount::Exact(3)));
        assert_eq!(
            RepeatCount::parse("2-5"),
            Some(RepeatCount::Range {
                min: Some(2),
                max: Some(5)
            })
        );
        assert_eq!(RepeatCount::parse("3-n").unwrap().to_wurcs(), "3-n");
    }

    #[test]
    fn test_parse_fragment() {
        let mut graph = parse_wurcs(
            "WURCS=2.0/6,11,10/[a2122h-1x_1-5_2*NCC/3=O][a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5][a2112h-1b_1-5][a1221m-1a_1-5]/1-2-3-4-2-5-4-2-6-2-5/a4-b1_a6-i1_b4-c1_c3-d1_c6-g1_d2-e1_e4-f1_g2-h1_j4-k1_j1-d4|d6|g4|g6}",
        )
        .unwrap();
        assert_eq!(graph.edge_count(), 9);
        assert_eq!(graph.undefined_linkages().len(), 1);
        let undefined = &graph.undefined_linkages()[0];
        assert_eq!(undefined.parents.len(), 2);
        assert_eq!(undefined.parents[0].positions.len(), 2);
        assert_eq!(undefined.parents[1].positions.len(), 2);

        let _ = graph.inner_mut();
        let generated = write_wurcs(&graph).unwrap();
        assert!(generated.contains('}'), "{generated}");
        assert!(generated.starts_with("WURCS=2.0/6,11,10/"), "{generated}");
        let reparsed = parse_wurcs(&generated).unwrap();
        assert_eq!(reparsed.edge_count(), 9);
        assert_eq!(reparsed.undefined_linkages().len(), 1);
    }

    #[test]
    fn test_parse_special_backbone() {
        let cases = vec![
            "WURCS=2.0/1,1,0/[a26h-1b_1-4_3*CO]/1/",
            "WURCS=2.0/1,1,0/[o261m_3*CO]/1/",
            "WURCS=2.0/3,3,2/[a261m-1a_1-4_3*C=O][a1211h-1a_1-5_2*NC][a222h-1b_1-4_1*N]/1-2-3/a2-b1_b3-c5*OPO*/3O/3=O",
            "WURCS=2.0/1,1,0/[Ad2dd22h_3-7_1*OC_6*N]/1/",
        ];
        for case in &cases {
            let result = parse_wurcs(case);
            assert!(result.is_ok(), "Failed on '{}': {:?}", case, result.err());
        }
    }

    #[test]
    fn test_parse_composition() {
        let mut graph = parse_wurcs(
            "WURCS=2.0/4,15,0+/[AUd21122h_5*NCC/3=O][uxxxxh_2*NCC/3=O][uxxxxh][u1221m]/1-2-2-2-2-2-2-3-3-3-4-4-4-4-4/",
        )
        .unwrap();
        assert_eq!(graph.node_count(), 15);
        assert_eq!(graph.edge_count(), 0);
        assert!(graph.is_composition());
        let _ = graph.inner_mut();
        let generated = write_wurcs(&graph).unwrap();
        assert!(generated.starts_with("WURCS=2.0/4,15,0+/"), "{generated}");
    }

    #[test]
    fn test_cyclic_closure_is_structural() {
        let mut graph = parse_wurcs(
            "WURCS=2.0/1,7,7/[a2122h-1a_1-5_2*OC_3*OC_6*N]/1-1-1-1-1-1-1/a1-g4_a4-b1_b4-c1_c4-d1_d4-e1_e4-f1_f4-g1",
        )
        .unwrap();
        assert_eq!(
            graph
                .inner()
                .edge_weights()
                .filter(|linkage| linkage.cyclic)
                .count(),
            1
        );
        let _ = graph.inner_mut();
        let generated = write_wurcs(&graph).unwrap();
        let reparsed = parse_wurcs(&generated).unwrap();
        assert_eq!(reparsed.node_count(), 7);
        assert_eq!(reparsed.edge_count(), 7);
        assert_eq!(
            reparsed
                .inner()
                .edge_weights()
                .filter(|linkage| linkage.cyclic)
                .count(),
            1
        );
    }

    #[test]
    fn test_write_wurcs_roundtrip_parsable() {
        let inputs = vec![
            "WURCS=2.0/1,1,0/[a2122h-1x_1-5_6*OPO/3O/3=O]/1/",
            "WURCS=2.0/1,1,0/[a2122h-1b_1-5]/1/",
            "WURCS=2.0/4,4,3/[a2122h-1b_1-5][a2112h-1b_1-5][a2122h-1a_1-5_2*NCC/3=O][Aad21122h-2a_2-6_5*NCC/3=O]/1-2-3-4/a4-b1_b4-c1_c4-d1",
        ];
        for input in inputs {
            let graph = match parse_wurcs(input) {
                Ok(g) => g,
                Err(e) => {
                    panic!("Failed to parse WURCS '{}': {:?}", input, e);
                }
            };
            let output = match write_wurcs(&graph) {
                Ok(o) => o,
                Err(e) => {
                    panic!("Failed to write WURCS for '{}': {:?}", input, e);
                }
            };
            let reparse = parse_wurcs(&output);
            assert!(
                reparse.is_ok(),
                "Roundtrip failed for '{}': output='{}', error={:?}",
                input,
                output,
                reparse.err()
            );
        }
    }

    #[test]
    fn test_parse_hlose_substituent() {
        let result = parse_wurcs("WURCS=2.0/1,1,0/[a26h-1b_1-4_3*CO]/1/");
        assert!(result.is_ok());
    }
}
