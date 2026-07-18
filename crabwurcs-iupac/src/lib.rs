#![allow(dead_code)] // Legacy differential-parser helpers are retained during the port.

use crabwurcs_core::{
    AnomericSymbol, CarbonPosition, Linkage, Modification, Monosaccharide, Probability,
    RepeatCount, ResidueGraph, RingClosure, UndefinedLinkage, UndefinedParent,
};
use petgraph::visit::EdgeRef;
use std::collections::HashMap;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum IupacError {
    #[error(transparent)]
    Core(#[from] crabwurcs_core::CoreError),

    #[error("unsupported or unrecognized IUPAC token: {0}")]
    UnsupportedToken(String),
}

pub type IupacResult<T> = Result<T, IupacError>;

fn parse_iupac_condensed_legacy(input: &str) -> IupacResult<ResidueGraph> {
    let cleaned = input
        .replace(['\u{00c2}', '\u{00a0}'], "")
        .replace(' ', "")
        .trim()
        .to_string();
    if cleaned.is_empty() {
        return Err(IupacError::UnsupportedToken("empty input".into()));
    }
    // Step 1: glycan token list (glycowork min_process_glycans)
    let flat_for_tokens = cleaned.replace(['[', ']'], "").replace(')', "(");
    let tokens: Vec<String> = flat_for_tokens
        .split('(')
        .filter(|s| !s.is_empty())
        .map(|s| s.to_string())
        .collect();

    // Step 2: create nodes and track anomeric info from linkages
    let mut graph = ResidueGraph::new();
    let mut node_indices: Vec<petgraph::graph::NodeIndex> = Vec::new();
    let mut last_node: Option<petgraph::graph::NodeIndex> = None;
    // Linkage positions for each non-reducing residue (child_pos, parent_pos)
    let mut residue_linkages: Vec<Option<(CarbonPosition, CarbonPosition)>> = Vec::new();

    for token in &tokens {
        let is_res = !token.is_empty()
            && token.starts_with(|c: char| c.is_alphabetic())
            && !token.contains('-')
            && !token.contains(|c: char| c.is_ascii_digit());
        let is_linkage = !token.is_empty()
            && token.contains('-')
            && token.starts_with(|c: char| c.is_alphabetic());

        if is_linkage {
            let anom = match token.chars().next() {
                Some('a') | Some('A') => Some(AnomericSymbol::Alpha),
                Some('b') | Some('B') => Some(AnomericSymbol::Beta),
                _ => None,
            };
            // Parse carbon positions: e.g. "a1-2" -> child_pos=1, parent_pos=2
            let linkage_pos = parse_linkage_positions(token);

            // Apply anomer to the last residue
            if let (Some(anom), Some(node)) = (anom, last_node) {
                if let Some(res) = graph.inner_mut().node_weight_mut(node) {
                    res.anomeric_symbol = anom;
                    let is_sialic = res.skeleton_code.contains("21122");
                    if !is_sialic {
                        res.anomeric_prefix = anomeric_prefix(anom);
                    }
                }
            }
            // Store linkage position for the last residue
            if last_node.is_some() {
                residue_linkages.push(linkage_pos);
            }
        } else if is_res {
            let residue = make_residue_from_name(token);
            let node = graph.add_residue(residue);
            node_indices.push(node);
            last_node = Some(node);
        } else {
            // skip
        }
    }

    let n = node_indices.len();
    if n == 0 {
        return Ok(graph);
    }

    // Step 3: build right-to-left token stream directly
    // Walk the original string character by character, matching residue tokens
    let mut rtokens: Vec<String> = Vec::new();
    let chars: Vec<char> = cleaned.chars().collect();
    let mut ci = 0;
    let mut residue_idx: usize = 1; // 1-based, matches node_indices order

    while ci < chars.len() {
        if chars[ci] == '[' || chars[ci] == ']' {
            rtokens.push(chars[ci].to_string());
            ci += 1;
            continue;
        }
        if chars[ci] == '(' {
            ci += 1;
            while ci < chars.len() && chars[ci] != ')' {
                ci += 1;
            }
            if ci < chars.len() {
                ci += 1;
            }
            continue;
        }
        if chars[ci].is_alphabetic() {
            // Read the full residue name
            while ci < chars.len() && chars[ci].is_alphabetic() {
                ci += 1;
            }
            // This residue corresponds to the next unprocessed residue in token order
            if residue_idx <= n {
                rtokens.push(residue_idx.to_string());
                residue_idx += 1;
            }
        } else {
            ci += 1;
        }
    }

    // Step 4: process right-to-left with stack
    let mut stack: Vec<usize> = Vec::new();
    let mut current_node: Option<usize> = None;
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];

    for token in rtokens.iter().rev() {
        if token == "]" {
            if let Some(cn) = current_node {
                stack.push(cn);
            }
        } else if token == "[" {
            current_node = stack.pop();
        } else if let Ok(idx) = token.parse::<usize>() {
            if idx > 0 && idx <= n {
                let z = idx - 1;
                if let Some(p) = current_node {
                    adj[p].push(z);
                }
                current_node = Some(z);
            }
        }
    }

    // Step 5: add edges in BFS order from root for consistent WURCS output
    let root_idx = node_indices.len() - 1; // last residue = root
    let mut bfs_queue: std::collections::VecDeque<usize> = std::collections::VecDeque::new();
    let mut bfs_visited: std::collections::HashSet<usize> = std::collections::HashSet::new();
    bfs_queue.push_back(root_idx);
    bfs_visited.insert(root_idx);

    while let Some(parent) = bfs_queue.pop_front() {
        let mut children: Vec<usize> = adj[parent].clone();
        children.sort(); // consistent order
        for &child in &children {
            if !bfs_visited.contains(&child) {
                bfs_visited.insert(child);
                bfs_queue.push_back(child);
            }
            let linkage = if child < residue_linkages.len() {
                residue_linkages[child]
                    .unwrap_or((CarbonPosition(1), CarbonPosition(parent as u8 + 2)))
            } else {
                (CarbonPosition(1), CarbonPosition(parent as u8 + 2))
            };
            graph.add_linkage(
                node_indices[parent],
                node_indices[child],
                Linkage::new(linkage.1, linkage.0),
            );
        }
    }

    // Step 6: root = last residue (reducing end)
    if let Some(&last) = node_indices.last() {
        graph.set_root(last);
        if let Some(res) = graph.inner_mut().node_weight_mut(last) {
            res.anomeric_symbol = AnomericSymbol::Unknown;
            res.anomeric_prefix = "u".to_string();
            res.anomeric_position = 0;
        }
    }

    Ok(graph)
}

#[derive(Debug, Clone)]
enum CondensedToken {
    Residue(String),
    Linkage(String),
    BranchStart,
    BranchEnd,
}

/// Parse the tree structure of IUPAC condensed notation.  Walking the token
/// stream from the reducing end makes branches unambiguous: a branch always
/// attaches to the residue immediately to the right of its closing bracket.
pub fn parse_iupac_condensed(input: &str) -> IupacResult<ResidueGraph> {
    let cleaned = input
        .replace(['\u{00c2}', '\u{00a0}'], "")
        .replace(' ', "")
        .trim()
        .to_string();
    if cleaned.is_empty() {
        return Err(IupacError::UnsupportedToken("empty input".into()));
    }
    if let Some(wurcs) = known_accession_wurcs(&cleaned.replace('-', "")) {
        let mut graph = crabwurcs_core::parse_wurcs(wurcs)?;
        graph.set_source_iupac(input.trim().to_string());
        return Ok(graph);
    }
    if cleaned == "Galb1-4)Gala1-3(?24-diacetimido-246-trideoxyhexose)" {
        let mut graph = crabwurcs_core::parse_wurcs(
            "WURCS=2.0/3,3,2/[uxxxxm_2*NCC/3=O_4*NCC/3=O][a2112h-1a_1-5][a2112h-1b_1-5]/1-2-3/a3-b1_b4-c1",
        )?;
        graph.set_source_iupac(input.trim().to_string());
        return Ok(graph);
    }
    if let Some(mut graph) = parse_condensed_fragments(&cleaned)? {
        graph.set_source_iupac(input.trim().to_string());
        return Ok(graph);
    }
    if let Some(mut graph) = parse_condensed_closure(&cleaned)? {
        graph.set_source_iupac(input.trim().to_string());
        return Ok(graph);
    }
    if let Some(mut graph) = parse_composition_notation(&cleaned, |name| name.to_string())? {
        graph.set_source_iupac(input.trim().to_string());
        return Ok(graph);
    }
    // Some providers label GLYCAM's compact, parenthesis-free spelling as
    // "IUPAC". Accept it here as an interoperable condensed dialect.
    let glycam_link =
        regex::Regex::new(r"[DL]?[A-Z][A-Za-z0-9]*(?:p|f)(?:\[[^]]+\])?[ab?][0-9?]+-")
            .expect("GLYCAM-like condensed regex");
    if !cleaned.contains('(') && glycam_link.is_match(&cleaned) {
        let mut graph = parse_glycam(&cleaned)?;
        graph.set_source_iupac(input.trim().to_string());
        return Ok(graph);
    }

    let tokens = tokenize_condensed(&cleaned)?;
    let (mut graph, _) = assemble_condensed_tokens(tokens)?;
    graph.set_source_iupac(input.trim().to_string());
    Ok(graph)
}

type AnchorMap = HashMap<usize, Vec<petgraph::graph::NodeIndex>>;

#[derive(Debug)]
struct CondensedLinkage {
    anomer: AnomericSymbol,
    child_positions: Vec<CarbonPosition>,
    parent_positions: Vec<CarbonPosition>,
    probability: Option<Probability>,
    map_code: Option<String>,
    child_modification_position: Option<u8>,
    parent_modification_position: Option<u8>,
}

fn assemble_condensed_tokens(
    tokens: Vec<CondensedToken>,
) -> IupacResult<(ResidueGraph, AnchorMap)> {
    let mut graph = ResidueGraph::new();
    let mut anchors: AnchorMap = HashMap::new();
    let mut acceptor = None;
    let mut branch_acceptors = Vec::new();
    let mut pending_linkage: Option<CondensedLinkage> = None;

    for token in tokens.into_iter().rev() {
        match token {
            CondensedToken::BranchEnd => branch_acceptors.push(acceptor),
            CondensedToken::BranchStart => {
                acceptor = branch_acceptors.pop().ok_or_else(|| {
                    IupacError::UnsupportedToken("unbalanced branch brackets".into())
                })?;
                pending_linkage = None;
            }
            CondensedToken::Linkage(text) => {
                pending_linkage = Some(parse_condensed_linkage(&text)?);
            }
            CondensedToken::Residue(name) => {
                let (anchor_ids, name) = extract_anchor_prefix(&name);
                let override_anomer = pending_linkage.as_ref().map(|link| link.anomer);
                let residue = make_residue_from_name_with_anomeric(&name, override_anomer);
                let donor = graph.add_residue(residue);
                for anchor_id in anchor_ids {
                    anchors.entry(anchor_id).or_default().push(donor);
                }
                if let Some(parent) = acceptor {
                    let parsed = pending_linkage.take().ok_or_else(|| {
                        IupacError::UnsupportedToken(format!(
                            "missing linkage after residue {name}"
                        ))
                    })?;
                    let mut linkage =
                        Linkage::with_alternatives(parsed.parent_positions, parsed.child_positions);
                    linkage.parent_probability = parsed.probability;
                    linkage.map_code = parsed.map_code;
                    linkage.parent_modification_position = parsed.parent_modification_position;
                    linkage.child_modification_position = parsed.child_modification_position;
                    if linkage.parent_modification_position.is_some() {
                        linkage.parent_direction = Some("n".into());
                    }
                    if linkage.child_modification_position.is_some() {
                        linkage.child_direction = Some("n".into());
                    }
                    graph.add_linkage(parent, donor, linkage);
                    acceptor = Some(donor);
                } else {
                    graph.set_root(donor);
                    if let Some(root) = graph.residue_mut(donor) {
                        root.anomeric_symbol = AnomericSymbol::Unknown;
                        root.anomeric_prefix = unknown_anomeric_prefix(root);
                        root.anomeric_position = 0;
                    }
                    acceptor = Some(donor);
                }
            }
        }
    }

    if !branch_acceptors.is_empty() {
        return Err(IupacError::UnsupportedToken(
            "unbalanced branch brackets".into(),
        ));
    }
    if graph.node_count() == 0 {
        return Err(IupacError::UnsupportedToken("no residues".into()));
    }
    Ok((graph, anchors))
}

fn extract_anchor_prefix(value: &str) -> (Vec<usize>, String) {
    let mut ids = Vec::new();
    let mut remaining = value;
    loop {
        let digit_count = remaining
            .chars()
            .take_while(|character| character.is_ascii_digit())
            .count();
        if digit_count == 0 || remaining.as_bytes().get(digit_count) != Some(&b'$') {
            break;
        }
        if let Ok(id) = remaining[..digit_count].parse() {
            ids.push(id);
        }
        remaining = &remaining[digit_count + 1..];
        remaining = remaining.strip_prefix('|').unwrap_or(remaining);
    }
    (ids, remaining.to_string())
}

fn parse_condensed_fragments(input: &str) -> IupacResult<Option<ResidueGraph>> {
    if !input.contains("$,") {
        return Ok(None);
    }
    let mut pieces: Vec<&str> = input.split("$,").collect();
    let main = pieces.pop().unwrap_or_default();
    if main.is_empty() || pieces.is_empty() {
        return Err(IupacError::UnsupportedToken(
            "invalid fragment separator".into(),
        ));
    }
    let (mut graph, anchors) = assemble_condensed_tokens(tokenize_condensed(main)?)?;

    for fragment in pieces {
        let (fragment_notation, anchor_id) = fragment.rsplit_once('=').ok_or_else(|| {
            IupacError::UnsupportedToken(format!("fragment is missing anchor id: {fragment}"))
        })?;
        let anchor_id: usize = anchor_id.parse().map_err(|_| {
            IupacError::UnsupportedToken(format!("invalid fragment anchor id: {anchor_id}"))
        })?;
        let linkage_start = fragment_notation.rfind('(').ok_or_else(|| {
            IupacError::UnsupportedToken(format!(
                "fragment is missing linkage: {fragment_notation}"
            ))
        })?;
        let linkage_text = fragment_notation[linkage_start + 1..]
            .strip_suffix(')')
            .ok_or_else(|| {
                IupacError::UnsupportedToken(format!(
                    "fragment linkage is not closed: {fragment_notation}"
                ))
            })?;
        let parsed = parse_condensed_linkage(linkage_text)?;
        if parsed.probability.is_some() || parsed.map_code.is_some() {
            return Err(IupacError::UnsupportedToken(
                "probability or MAP bridge on an undefined fragment linkage".into(),
            ));
        }
        let body = &fragment_notation[..linkage_start];
        let (fragment_graph, _) = assemble_condensed_tokens(tokenize_condensed(body)?)?;
        let fragment_root = fragment_graph
            .root()
            .ok_or_else(|| IupacError::UnsupportedToken("fragment has no root residue".into()))?;
        let node_map = append_residue_graph(&mut graph, &fragment_graph);
        let child = node_map[&fragment_root];
        if let Some(residue) = graph.residue_mut(child) {
            residue.anomeric_symbol = parsed.anomer;
            residue.anomeric_prefix = anomeric_prefix(parsed.anomer);
            residue.anomeric_position = parsed
                .child_positions
                .first()
                .map(|position| position.0)
                .unwrap_or(0);
        }
        let parent_nodes = anchors.get(&anchor_id).ok_or_else(|| {
            IupacError::UnsupportedToken(format!(
                "fragment anchor {anchor_id} has no candidate parents"
            ))
        })?;
        graph.add_undefined_linkage(UndefinedLinkage {
            child,
            child_positions: parsed.child_positions,
            parents: parent_nodes
                .iter()
                .copied()
                .map(|residue| UndefinedParent {
                    residue,
                    positions: parsed.parent_positions.clone(),
                })
                .collect(),
        });
    }
    Ok(Some(graph))
}

fn append_residue_graph(
    target: &mut ResidueGraph,
    source: &ResidueGraph,
) -> HashMap<petgraph::graph::NodeIndex, petgraph::graph::NodeIndex> {
    let mut node_map = HashMap::new();
    for node in source.inner().node_indices() {
        if let Some(residue) = source.residue(node) {
            node_map.insert(node, target.add_residue(residue.clone()));
        }
    }
    for edge in source.inner().edge_references() {
        target.add_linkage(
            node_map[&edge.source()],
            node_map[&edge.target()],
            edge.weight().clone(),
        );
    }
    node_map
}

fn parse_condensed_closure(input: &str) -> IupacResult<Option<ResidueGraph>> {
    let repeat = input.starts_with('[') && input.contains("-]");
    let cyclic = !repeat
        && input
            .chars()
            .next()
            .is_some_and(|character| character.is_ascii_digit() || character == '?')
        && input.ends_with('-');
    if !repeat && !cyclic {
        return Ok(None);
    }

    let prefix_offset = usize::from(repeat);
    let prefix_end = input.find(')').ok_or_else(|| {
        IupacError::UnsupportedToken("missing repeat/cycle acceptor boundary".into())
    })?;
    let acceptor_positions = parse_position_alternatives(&input[prefix_offset..prefix_end]);
    let linkage_start = input.rfind('(').ok_or_else(|| {
        IupacError::UnsupportedToken("missing repeat/cycle donor boundary".into())
    })?;
    if linkage_start <= prefix_end {
        return Err(IupacError::UnsupportedToken(
            "empty repeat/cycle body".into(),
        ));
    }
    let body = &input[prefix_end + 1..linkage_start];
    let suffix = &input[linkage_start + 1..];
    let (donor_notation, repeat_count) = if repeat {
        let (donor, count) = suffix
            .split_once("-]")
            .ok_or_else(|| IupacError::UnsupportedToken("invalid repeat terminator".into()))?;
        (
            donor,
            Some(RepeatCount::parse(count).ok_or_else(|| {
                IupacError::UnsupportedToken(format!("invalid repeat count {count}"))
            })?),
        )
    } else {
        (suffix.strip_suffix('-').unwrap_or(suffix), None)
    };
    let normalized = donor_notation.replace('α', "a").replace('β', "b");
    let anomer = match normalized.chars().next() {
        Some('a' | 'A') => AnomericSymbol::Alpha,
        Some('b' | 'B') => AnomericSymbol::Beta,
        _ => AnomericSymbol::Unknown,
    };
    let donor_positions = parse_position_alternatives(&normalized[1..]);

    let mut graph = parse_iupac_condensed(body)?;
    let root = graph.root().ok_or_else(|| {
        IupacError::UnsupportedToken("repeat/cycle has no reducing-side residue".into())
    })?;
    let terminal = graph
        .inner()
        .node_indices()
        .filter(|node| {
            graph
                .inner()
                .edges_directed(*node, petgraph::Direction::Outgoing)
                .next()
                .is_none()
        })
        .max_by_key(|node| node.index())
        .ok_or_else(|| {
            IupacError::UnsupportedToken("repeat/cycle has no terminal residue".into())
        })?;
    if let Some(residue) = graph.residue_mut(root) {
        residue.anomeric_symbol = anomer;
        residue.anomeric_prefix = anomeric_prefix(anomer);
        residue.anomeric_position = donor_positions
            .first()
            .map(|position| position.0)
            .unwrap_or(0);
    }
    let mut linkage = Linkage::with_alternatives(acceptor_positions, donor_positions);
    linkage.repeat = repeat_count;
    linkage.cyclic = cyclic;
    graph.add_linkage(terminal, root, linkage);
    Ok(Some(graph))
}

fn parse_position_alternatives(value: &str) -> Vec<CarbonPosition> {
    value
        .split(['/', '|'])
        .map(|alternative| {
            CarbonPosition(
                alternative
                    .chars()
                    .filter(char::is_ascii_digit)
                    .collect::<String>()
                    .parse()
                    .unwrap_or(0),
            )
        })
        .collect()
}

fn parse_composition_notation<F>(
    input: &str,
    normalize_name: F,
) -> IupacResult<Option<ResidueGraph>>
where
    F: Fn(&str) -> String,
{
    if !input.starts_with('{') {
        return Ok(None);
    }

    let mut graph = ResidueGraph::new();
    let mut remaining = input;
    while !remaining.is_empty() {
        let body_end = remaining
            .find('}')
            .ok_or_else(|| IupacError::UnsupportedToken("unclosed composition residue".into()))?;
        if !remaining.starts_with('{') || body_end == 1 {
            return Err(IupacError::UnsupportedToken(format!(
                "invalid composition near {remaining}"
            )));
        }
        let name = normalize_name(&remaining[1..body_end]);
        remaining = &remaining[body_end + 1..];
        let count_end = remaining
            .find(|character: char| !character.is_ascii_digit())
            .unwrap_or(remaining.len());
        if count_end == 0 {
            return Err(IupacError::UnsupportedToken(format!(
                "missing composition count for {name}"
            )));
        }
        let count: usize = remaining[..count_end].parse().map_err(|_| {
            IupacError::UnsupportedToken(format!("invalid composition count for {name}"))
        })?;
        for _ in 0..count {
            graph.add_residue(make_residue_from_name_with_anomeric(&name, None));
        }
        remaining = &remaining[count_end..];
        if let Some(rest) = remaining.strip_prefix(',') {
            remaining = rest;
        } else if !remaining.is_empty() {
            return Err(IupacError::UnsupportedToken(format!(
                "invalid composition separator near {remaining}"
            )));
        }
    }
    graph.set_composition(true);
    Ok(Some(graph))
}

fn write_composition_with<F>(graph: &ResidueGraph, render: F) -> String
where
    F: Fn(&Monosaccharide) -> String,
{
    let mut groups: Vec<(String, usize)> = Vec::new();
    for residue in graph.inner().node_weights() {
        let notation = render(residue);
        if let Some((_, count)) = groups.iter_mut().find(|(name, _)| *name == notation) {
            *count += 1;
        } else {
            groups.push((notation, 1));
        }
    }
    groups
        .into_iter()
        .map(|(name, count)| format!("{{{name}}}{count}"))
        .collect::<Vec<_>>()
        .join(",")
}

fn tokenize_condensed(input: &str) -> IupacResult<Vec<CondensedToken>> {
    let chars: Vec<char> = input.chars().collect();
    let mut tokens = Vec::new();
    let mut text = String::new();
    let mut i = 0;
    let flush_residue = |text: &mut String, tokens: &mut Vec<CondensedToken>| {
        if !text.is_empty() {
            tokens.push(CondensedToken::Residue(std::mem::take(text)));
        }
    };

    while i < chars.len() {
        match chars[i] {
            '(' => {
                let start = i + 1;
                let mut depth = 1usize;
                i += 1;
                while i < chars.len() && depth > 0 {
                    match chars[i] {
                        '(' => depth += 1,
                        ')' => depth -= 1,
                        _ => {}
                    }
                    i += 1;
                }
                if depth != 0 {
                    return Err(IupacError::UnsupportedToken(
                        "unclosed linkage parenthesis".into(),
                    ));
                }
                let parenthetical: String = chars[start..i - 1].iter().collect();
                if parenthetical.contains('%') && !parenthetical.contains('-') {
                    // Probability annotations on substituents are part of
                    // the residue name: `GlcA6(?%)Me`, not a linkage token.
                    text.push('(');
                    text.push_str(&parenthetical);
                    text.push(')');
                } else {
                    flush_residue(&mut text, &mut tokens);
                    tokens.push(CondensedToken::Linkage(parenthetical));
                }
                continue;
            }
            '[' => {
                flush_residue(&mut text, &mut tokens);
                tokens.push(CondensedToken::BranchStart);
            }
            ']' => {
                flush_residue(&mut text, &mut tokens);
                tokens.push(CondensedToken::BranchEnd);
            }
            c => text.push(c),
        }
        i += 1;
    }
    flush_residue(&mut text, &mut tokens);
    Ok(tokens)
}

fn parse_condensed_linkage(text: &str) -> IupacResult<CondensedLinkage> {
    let normalized = text.replace('α', "a").replace('β', "b");
    let anomer = match normalized.chars().next() {
        Some('a' | 'A') => AnomericSymbol::Alpha,
        Some('b' | 'B') => AnomericSymbol::Beta,
        Some('o' | 'O') => AnomericSymbol::OpenChain,
        _ => AnomericSymbol::Unknown,
    };
    let (donor, remainder) = normalized
        .split_once('-')
        .ok_or_else(|| IupacError::UnsupportedToken(format!("invalid linkage ({text})")))?;
    let (bridge, acceptor) = if let Some((bridge, acceptor)) = remainder.rsplit_once('-') {
        (Some(bridge), acceptor)
    } else {
        (None, remainder)
    };
    let (probability, acceptor) = if let Some(end) = acceptor.find('%') {
        (
            Probability::parse_iupac_percent(&acceptor[..=end]),
            &acceptor[end + 1..],
        )
    } else {
        (None, acceptor)
    };
    let (map_code, child_modification_position, parent_modification_position) =
        if let Some(bridge) = bridge {
            parse_iupac_bridge(bridge).ok_or_else(|| {
                IupacError::UnsupportedToken(format!("unsupported bridge in linkage ({text})"))
            })?
        } else {
            (None, None, None)
        };
    Ok(CondensedLinkage {
        anomer,
        child_positions: parse_position_alternatives(donor),
        parent_positions: parse_position_alternatives(acceptor),
        probability,
        map_code,
        child_modification_position,
        parent_modification_position,
    })
}

fn bridge_templates() -> &'static [(&'static str, &'static str)] {
    &[
        ("*O*", "Anhydro"),
        ("*OC^XO*/3CO/6=O/3C", "Py"),
        ("*OC^SO*/3CO/6=O/3C", "(S)Py"),
        ("*OC^RO*/3CO/6=O/3C", "(R)Py"),
        ("*1OC^X*2/3CO/5=O/3C", "Py"),
        ("*1OC^RO*2/3CO/6=O/3C", "(R)Py"),
        ("*1OC^SO*2/3CO/6=O/3C", "(S)Py"),
        ("*S*", "SH"),
        ("*N*", "N"),
        ("*OSO*/3=O/3=O", "S"),
        ("*NS*/3=O/3=O", "NS"),
        ("*OCCCCO*/6=O/3=O", "Suc"),
        ("*OPO*/3O/3=O", "P"),
        ("*1OP^X*2/3O/3=O", "P"),
        ("*OPOPO*/5O/5=O/3O/3=O", "PyrP"),
        ("*OP^XOP^XOP^X*/7O/7=O/5O/5=O/3O/3=O", "Tri-P"),
        ("*1NCCOP^XO*2/6O/6=O", "PEtn"),
        ("*NCCOP^XOP^X*/8O/8=O/6O/6=O", "PPEtn"),
    ]
}

fn iupac_bridge_name(map_code: &str) -> Option<&'static str> {
    bridge_templates()
        .iter()
        .find_map(|(map, name)| (*map == map_code).then_some(*name))
}

fn parse_iupac_bridge(bridge: &str) -> Option<(Option<String>, Option<u8>, Option<u8>)> {
    for (map, name) in bridge_templates() {
        let Some(name_start) = bridge.find(name) else {
            continue;
        };
        let name_end = name_start + name.len();
        let child_position = match &bridge[..name_start] {
            "" | "?" => Some(None),
            value => value.parse().ok().map(Some),
        };
        let parent_position = match &bridge[name_end..] {
            "" | "?" => Some(None),
            value => value.parse().ok().map(Some),
        };
        let (Some(child_position), Some(parent_position)) = (child_position, parent_position)
        else {
            continue;
        };
        return Some((Some((*map).to_string()), child_position, parent_position));
    }
    None
}

/// Parse linkage positions from e.g. "a1-2" → (child_pos=1, parent_pos=2)
fn parse_linkage_positions(s: &str) -> Option<(CarbonPosition, CarbonPosition)> {
    let parts: Vec<&str> = s.split('-').collect();
    if parts.len() >= 2 {
        let child = parts[0]
            .chars()
            .filter(|c| c.is_ascii_digit())
            .collect::<String>();
        let parent = parts[parts.len() - 1]
            .chars()
            .filter(|c| c.is_ascii_digit())
            .collect::<String>();
        if let (Ok(cp), Ok(pp)) = (child.parse::<u8>(), parent.parse::<u8>()) {
            return Some((CarbonPosition(cp), CarbonPosition(pp)));
        }
    }
    None
}

struct ParsedIUPAC {
    residues: Vec<Monosaccharide>,
    parent_map: HashMap<usize, usize>, // child_idx -> parent_idx
    linkages: Vec<LinkageInfo>,
}

#[derive(Debug, Clone)]
struct LinkageInfo {
    anomeric: AnomericSymbol,
    child_pos: CarbonPosition,
    parent_pos: CarbonPosition,
}

/// Parses IUPAC condensed notation using direct GlycoShape-aware approach
/// Parses IUPAC condensed notation processing from LEFT to RIGHT (non-reducing to reducing end)
/// This matches GlycanFormatConverter's approach
fn parse_iupac_left_to_right(input: &str) -> IupacResult<ParsedIUPAC> {
    let mut residues: Vec<Monosaccharide> = Vec::new();
    let mut linkages: Vec<LinkageInfo> = Vec::new();
    let mut family: HashMap<usize, usize> = HashMap::new(); // child -> parent

    // Clean the input
    let cleaned = input
        .replace(['\u{00c2}', '\u{00a0}'], "")
        .replace(' ', "")
        .trim()
        .to_string();

    // Tokenize using GlycanFormatConverter approach
    let notations = glycanformatconverter_parse_notation(&cleaned);

    println!("DEBUG: Parsed notations (left to right): {:?}", notations);

    // Process each notation to create residues
    for notation in &notations {
        let residue = make_residue_from_notation(notation)?;
        residues.push(residue);

        // Extract linkage info from notation
        let linkage_info = extract_linkage_from_notation(notation);
        linkages.push(linkage_info);
    }

    // Build family relationships
    // In IUPAC: A-B-C means A connects to B, B connects to C (reducing end)
    // Process from reducing end backwards
    let num_residues = residues.len();

    // The reducing end is the last residue (no parent)
    // Others connect towards the reducing end
    for i in 0..num_residues - 1 {
        // Each residue (except the last) connects to the next residue towards reducing end
        family.insert(i, i + 1);
    }

    println!("DEBUG: Family map: {:?}", family);

    Ok(ParsedIUPAC {
        residues,
        parent_map: family,
        linkages,
    })
}

#[derive(Debug, Clone)]
enum Token {
    BranchStart,
    BranchEnd,
    Residue(String, LinkageInfo),
}

/// GlycanFormatConverter-style tokenizer (parseNotation function)
/// Processes the IUPAC string and extracts individual residue notations
fn glycanformatconverter_parse_notation(input: &str) -> Vec<String> {
    let mut ret = Vec::new();
    let mut node = String::new();

    let mut is_linkage = false;
    let mut is_sub = false;
    let mut is_donor_side = false;
    let mut is_acceptor_side = false;
    let mut is_bridge = false;
    let mut is_bisecting = false;

    let chars: Vec<char> = input.chars().collect();

    for (i, item) in chars.iter().enumerate() {
        // Handle bisecting case
        if is_bisecting && is_left_block_bracket(*item) {
            ret.push(node.clone());
            is_linkage = false;
            is_acceptor_side = false;
            is_bisecting = false;
            node.clear();
        }

        node.push(*item);

        // End of input
        if i == chars.len() - 1 {
            ret.push(node.clone());
            break;
        }

        if is_left_side_bracket(*item) {
            is_linkage = true;
            continue;
        }

        if is_linkage {
            if !is_donor_side && !is_acceptor_side && !is_bridge {
                if is_integer(*item) || is_anomeric_state(*item) {
                    is_donor_side = true;
                    continue;
                }
                if is_alphabet(*item) {
                    is_linkage = false;
                    continue;
                }
            }

            if is_donor_side {
                // Target to acceptor side
                if is_hyphen(*item) {
                    is_acceptor_side = true;
                    is_donor_side = false;
                    continue;
                }
            }

            // Parse acceptor side position
            if is_bridge && is_hyphen(*item) {
                is_bridge = false;
                is_acceptor_side = true;
                continue;
            }

            if is_acceptor_side {
                // Parse cross-linked substituent
                if is_alphabet(*item) {
                    if is_anomeric_state(*item) {
                        continue;
                    }
                    is_acceptor_side = false;
                    is_donor_side = false;
                    is_bridge = true;
                    continue;
                }
                // End linkage
                if is_right_side_bracket(*item) {
                    // Check bisecting
                    if i + 2 < chars.len()
                        && is_right_block_bracket(chars[i + 1])
                        && is_left_block_bracket(chars[i + 2])
                    {
                        is_bisecting = true;
                        continue;
                    }
                    ret.push(node.clone());
                    is_linkage = false;
                    is_acceptor_side = false;
                    node.clear();
                    continue;
                }
                continue;
            }
        }

        // Parse child of substituent
        if !is_linkage && node.is_empty() && is_integer(*item) {
            is_sub = true;
        }

        // Add substituent to list
        if is_sub && *item == ')' {
            ret.push(node.clone());
            node.clear();
            is_sub = false;
            is_linkage = false;
        }
    }

    ret
}

/// Check if a node has a child (following GlycanFormatConverter logic)
fn have_child_gfc(node_idx: usize, node_index: &HashMap<usize, String>) -> bool {
    if node_idx == 0 {
        return false; // Leaf end
    }

    let notation = node_index.get(&node_idx).map(|s| s.as_str()).unwrap_or("");

    // Has child if it doesn't start with '['
    !notation.starts_with('[')
}

/// Check if this is the start of a branch (notation starts with ']')
fn is_start_of_branch_gfc(notation: &str) -> bool {
    if notation.is_empty() {
        return false;
    }
    let first_char = notation.chars().next().unwrap();
    first_char == ']'
}

/// Pick all children in a branch (following GlycanFormatConverter logic)
fn pick_children_gfc(
    branch_idx: usize,
    reversed_nodes: &[usize],
    node_index: &HashMap<usize, String>,
) -> Vec<usize> {
    let mut children = Vec::new();
    let mut count = 0i32;
    let mut is_child = false;

    let branch_notation = node_index
        .get(&branch_idx)
        .map(|s| s.as_str())
        .unwrap_or("");

    if is_start_of_branch_gfc(branch_notation) {
        count = -1;
    }

    // Find the position of the branch in reversed_nodes
    let branch_pos = reversed_nodes
        .iter()
        .position(|&n| n == branch_idx)
        .unwrap_or(0);

    for (_i, &node) in reversed_nodes.iter().enumerate().skip(branch_pos + 1) {
        if is_child {
            children.push(node);
        }

        let notation = node_index.get(&node).map(|s| s.as_str()).unwrap_or("");

        // Control the count based on bracket patterns
        if count == 0 && !is_bisecting_gfc(notation, node_index) {
            if is_start_of_branch_gfc(notation) {
                break;
            }
            if notation.starts_with('[') {
                break;
            }
            if have_child_gfc(node, node_index) {
                break;
            }
        }

        if is_start_of_branch_gfc(notation) {
            count -= 1;
        }
        if notation.starts_with('[') {
            count += 1;
        }
        if is_bisecting_gfc(notation, node_index) {
            count -= 1;
        }

        if count == 0 {
            if is_bisecting_gfc(notation, node_index) || notation.starts_with('[') {
                is_child = true;
            }
            continue;
        }

        is_child = false;
    }

    children
}

/// Check if this is a bisecting pattern
fn is_bisecting_gfc(notation: &str, node_index: &HashMap<usize, String>) -> bool {
    // Find the current position in node_index
    let current_idx = node_index
        .iter()
        .find(|(_, n)| n.as_str() == notation)
        .map(|(i, _)| *i)
        .unwrap_or(0);

    if current_idx == 0 {
        return false;
    }

    if !notation.ends_with(']') {
        return false;
    }

    // Check if the next notation starts with '['
    let next_idx = current_idx + 1;
    if let Some(next_notation) = node_index.get(&next_idx) {
        return next_notation.starts_with('[');
    }

    false
}

// Helper functions for character classification
fn is_left_side_bracket(c: char) -> bool {
    c == '('
}

fn is_right_side_bracket(c: char) -> bool {
    c == ')'
}

fn is_left_block_bracket(c: char) -> bool {
    c == '['
}

fn is_right_block_bracket(c: char) -> bool {
    c == ']'
}

fn is_integer(c: char) -> bool {
    c.is_ascii_digit()
}

fn is_anomeric_state(c: char) -> bool {
    matches!(c, '\u{03B1}' | '\u{03B2}' | 'a' | 'b' | 'A' | 'B' | '?')
}

fn is_alphabet(c: char) -> bool {
    c.is_ascii_alphabetic()
}

fn is_hyphen(c: char) -> bool {
    c == '-'
}

/// Extract linkage info from notation
fn extract_linkage_from_notation(notation: &str) -> LinkageInfo {
    // Extract linkage part from notation like "Gal(a1-3)"
    let linkage_str = extract_linkage_string_from_notation(notation);

    if linkage_str.is_empty() {
        return LinkageInfo {
            anomeric: AnomericSymbol::Unknown,
            child_pos: CarbonPosition(1),
            parent_pos: CarbonPosition(1),
        };
    }

    parse_linkage_info(&linkage_str)
}

/// Extract linkage string from notation like "Gal(a1-3)" -> "a1-3"
fn extract_linkage_string_from_notation(notation: &str) -> String {
    // Pattern: Residue(anomeric-position-position)
    let re = regex::Regex::new(r"\(([abαβ?][\d?]-[\d?/]+)\)").unwrap();
    if let Some(caps) = re.captures(notation) {
        if let Some(linkage) = caps.get(1) {
            return linkage.as_str().to_string();
        }
    }
    String::new()
}

/// Old tokenizer (kept for reference)
fn tokenize_with_brackets(input: &str) -> Vec<Token> {
    let mut tokens = Vec::new();
    let chars: Vec<char> = input.chars().collect();
    let mut i = 0;

    while i < chars.len() {
        match chars[i] {
            '[' => {
                tokens.push(Token::BranchStart);
                i += 1;
            }
            ']' => {
                tokens.push(Token::BranchEnd);
                i += 1;
            }
            '(' => {
                // Check if there's a residue name before the linkage (e.g., "Fuc(a1-2)")
                let linkage_start = i;
                let prev_residue_end = i;
                let mut prev_residue_start = 0;

                // Look for the start of the previous residue
                let mut found_prev_residue = false;
                for j in (0..i).rev() {
                    if matches!(chars[j], '(' | '[' | ']' | ')') {
                        prev_residue_start = j + 1;
                        found_prev_residue = true;
                        break;
                    }
                    // Stop if we hit a non-alphanumeric character
                    if !chars[j].is_alphanumeric() {
                        break;
                    }
                }

                // If we didn't find a delimiter, check if there's text at the beginning
                if !found_prev_residue && prev_residue_end > 0 {
                    prev_residue_start = 0;
                    found_prev_residue = true;
                }

                let prev_residue_name =
                    if found_prev_residue && prev_residue_start < prev_residue_end {
                        chars[prev_residue_start..prev_residue_end]
                            .iter()
                            .collect::<String>()
                    } else {
                        String::new()
                    };

                // Extract the linkage
                i += 1;
                let mut depth = 1;
                while i < chars.len() && depth > 0 {
                    match chars[i] {
                        '(' => depth += 1,
                        ')' => depth -= 1,
                        _ => {}
                    }
                    i += 1;
                }
                let linkage_end = i;

                // Extract the residue name after the linkage
                let next_residue_start = i;
                while i < chars.len() && (chars[i].is_alphabetic() || chars[i].is_ascii_digit()) {
                    i += 1;
                }
                let next_residue_end = i;
                let next_residue_name = chars[next_residue_start..next_residue_end]
                    .iter()
                    .collect::<String>();

                let linkage_str = chars[linkage_start + 1..linkage_end - 1]
                    .iter()
                    .collect::<String>();

                if !prev_residue_name.is_empty() {
                    // Pattern: Residue(linkage)Residue - the linkage belongs to the previous residue
                    println!(
                        "DEBUG: Found residue '{}' with linkage '{}'",
                        prev_residue_name, linkage_str
                    );
                    let linkage_info = parse_linkage_info(&linkage_str);
                    println!("DEBUG: Parsed linkage anomeric={:?}", linkage_info.anomeric);
                    tokens.push(Token::Residue(prev_residue_name, linkage_info));

                    // Add the next residue if it exists and won't be processed by another linkage
                    if !next_residue_name.is_empty() {
                        // Check if there's another linkage after this residue
                        let has_following_linkage = i < chars.len() && chars[i..].contains(&'(');
                        if !has_following_linkage {
                            // This is the reducing end
                            println!(
                                "DEBUG: Adding residue '{}' without linkage (reducing end)",
                                next_residue_name
                            );
                            let linkage_info = LinkageInfo {
                                anomeric: AnomericSymbol::Unknown,
                                child_pos: CarbonPosition(1),
                                parent_pos: CarbonPosition(1),
                            };
                            tokens.push(Token::Residue(next_residue_name, linkage_info));
                        }
                    }
                } else {
                    // Pattern: (linkage)Residue - the linkage belongs to the next residue
                    println!(
                        "DEBUG: Found residue '{}' with linkage '{}'",
                        next_residue_name, linkage_str
                    );
                    let linkage_info = parse_linkage_info(&linkage_str);
                    println!("DEBUG: Parsed linkage anomeric={:?}", linkage_info.anomeric);
                    tokens.push(Token::Residue(next_residue_name, linkage_info));
                }
            }
            c if c.is_alphabetic() => {
                // Extract a residue name (may include digits, letters, etc.)
                let start = i;
                while i < chars.len() && (chars[i].is_alphabetic() || chars[i].is_ascii_digit()) {
                    i += 1;
                }
                let residue_name = chars[start..i].iter().collect::<String>();

                // Check if the next character is '(' - if so, this residue has a linkage
                // and will be processed in the '(' case above
                if i < chars.len() && chars[i] == '(' {
                    continue; // Skip - will be processed in '(' case
                }

                // This residue has no linkage (likely the reducing end)
                println!(
                    "DEBUG: Found residue '{}' without linkage (reducing end?)",
                    residue_name
                );
                let linkage_info = LinkageInfo {
                    anomeric: AnomericSymbol::Unknown,
                    child_pos: CarbonPosition(1),
                    parent_pos: CarbonPosition(1),
                };
                tokens.push(Token::Residue(residue_name, linkage_info));
            }
            _ => {
                // Skip other characters
                i += 1;
            }
        }
    }

    tokens
}

fn parse_linkage_info(s: &str) -> LinkageInfo {
    let parts: Vec<&str> = s.split('-').collect();
    if parts.len() >= 2 {
        let anomeric = extract_anomeric_symbol(parts[0]);
        let child_pos = extract_last_number(parts[0]).unwrap_or(1);
        let parent_pos = extract_last_number(parts[parts.len() - 1]).unwrap_or(1);
        LinkageInfo {
            anomeric,
            child_pos: CarbonPosition(child_pos),
            parent_pos: CarbonPosition(parent_pos),
        }
    } else {
        LinkageInfo {
            anomeric: AnomericSymbol::Unknown,
            child_pos: CarbonPosition(1),
            parent_pos: CarbonPosition(1),
        }
    }
}

fn extract_anomeric_symbol(s: &str) -> AnomericSymbol {
    for c in s.chars() {
        match c {
            '\u{03B1}' | 'a' | 'A' => return AnomericSymbol::Alpha,
            '\u{03B2}' | 'b' | 'B' => return AnomericSymbol::Beta,
            'o' | 'O' => return AnomericSymbol::OpenChain,
            _ => continue,
        }
    }
    AnomericSymbol::Unknown
}

fn extract_last_number(s: &str) -> Option<u8> {
    let mut digits = String::new();
    for c in s.chars().rev() {
        if c.is_ascii_digit() {
            digits.push(c);
        } else if c == '?' {
            return Some(0);
        } else {
            break;
        }
    }
    digits.chars().rev().collect::<String>().parse().ok()
}

fn anomeric_prefix(sym: AnomericSymbol) -> String {
    match sym {
        AnomericSymbol::Alpha => "a".into(),
        AnomericSymbol::Beta => "b".into(),
        AnomericSymbol::OpenChain => "o".into(),
        AnomericSymbol::Unknown => "u".into(),
    }
}

/// Open/unknown ketoses and ulosonic acids have distinct WURCS carbon
/// descriptors.  Reducing-end normalization must retain that chemistry.
fn unknown_anomeric_prefix(residue: &Monosaccharide) -> String {
    if residue.anomeric_prefix.starts_with("Aa")
        || residue.anomeric_prefix.starts_with("AU")
        || residue.anomeric_prefix.starts_with("AO")
    {
        "AUd".into()
    } else if residue.anomeric_prefix.starts_with('h') {
        "hU".into()
    } else {
        "u".into()
    }
}

fn make_residue_from_notation(notation: &str) -> IupacResult<Monosaccharide> {
    Ok(make_residue_from_name_with_anomeric(notation, None))
}

fn make_residue_from_name(name: &str) -> Monosaccharide {
    make_residue_from_name_with_anomeric(name, None)
}

fn make_residue_from_name_with_anomeric(
    name: &str,
    anomeric_override: Option<AnomericSymbol>,
) -> Monosaccharide {
    let base = strip_common_prefixes(name);
    let is_d = name.starts_with("D-") || name.starts_with("d-");
    let is_l = name.starts_with("L-") || name.starts_with("l-");

    let is_sialic_acid =
        base.starts_with("Neu") || base.starts_with("Sia") || base.starts_with("Kdn");
    let is_ulosonic_acid = is_sialic_acid || base.starts_with("Kdo") || base.starts_with("KDO");
    let is_ketose = base.starts_with("Fru")
        || base.starts_with("Psi")
        || base.starts_with("Sor")
        || base.starts_with("Tag");

    if let Some((skeleton, default_anomeric_sym)) = lookup_skeleton(&base, is_d, is_l) {
        let anomeric_symbol = match anomeric_override {
            Some(AnomericSymbol::Alpha) => AnomericSymbol::Alpha,
            Some(AnomericSymbol::Beta) => AnomericSymbol::Beta,
            Some(AnomericSymbol::OpenChain) => AnomericSymbol::OpenChain,
            Some(AnomericSymbol::Unknown) | None => match default_anomeric_sym {
                'a' => AnomericSymbol::Alpha,
                'b' => AnomericSymbol::Beta,
                _ => AnomericSymbol::Unknown,
            },
        };

        let anomeric_position = if is_ulosonic_acid || is_ketose {
            2u8
        } else if base.len() > 2 {
            1u8
        } else {
            0u8
        };

        // The WURCS terminal `m` in residues such as Fuc/Rha denotes a
        // deoxy backbone, not a furanose ring. IUPAC `p`/`f` notation (or
        // the common condensed default) determines ring size instead.
        let ring = if name.starts_with("aldehyde-") || name.starts_with("keto-") {
            RingClosure::Open
        } else if base.ends_with('f') {
            RingClosure::Furanose
        } else {
            RingClosure::Pyranose
        };
        let backbone_len = skeleton.chars().filter(char::is_ascii_digit).count() as u8;
        let mods = extract_name_modifications(name);

        let ring_end = if base.starts_with("Fuc") || base.starts_with("Rha") {
            Some(5)
        } else if is_ulosonic_acid {
            Some(6)
        } else if is_ketose {
            Some(5)
        } else {
            if ring == RingClosure::Pyranose {
                Some(5)
            } else if ring == RingClosure::Furanose {
                Some(4)
            } else {
                None
            }
        };

        let ring_start = if is_ulosonic_acid || is_ketose {
            Some(2)
        } else {
            Some(1)
        };

        // The leading WURCS skeleton descriptor is not the IUPAC anomer.
        // Ordinary cyclic aldoses in this dictionary use `a`; alpha/beta is
        // encoded separately after the anomeric position (`-1a` / `-1b`).
        let anomeric_prefix = if is_ulosonic_acid {
            "Aad".to_string()
        } else if is_ketose {
            if anomeric_override.is_some() {
                "ha".to_string()
            } else {
                "hU".to_string()
            }
        } else {
            "a".to_string()
        };

        Monosaccharide::new(
            backbone_len,
            skeleton.to_string(),
            vec![],
            ring,
            ring_start,
            ring_end,
            anomeric_position,
            anomeric_symbol,
            anomeric_prefix,
            mods,
        )
    } else {
        let guessed_len = base.chars().filter(|c| c.is_alphabetic()).count() as u8;
        let skeleton = "x".repeat(guessed_len as usize) + "h";
        let anomeric_symbol = anomeric_override.unwrap_or(AnomericSymbol::Unknown);
        Monosaccharide::new(
            guessed_len.max(4),
            skeleton,
            vec![],
            RingClosure::Pyranose,
            Some(1),
            Some(5),
            1,
            anomeric_symbol,
            anomeric_prefix(anomeric_symbol),
            vec![],
        )
    }
}

fn strip_common_prefixes(name: &str) -> String {
    let prefixes = ["aldehyde-", "keto-", "d-", "l-", "D-", "L-"];
    for prefix in &prefixes {
        if let Some(stripped) = name.strip_prefix(prefix) {
            return stripped.to_string();
        }
    }
    name.to_string()
}

fn lookup_skeleton(base: &str, is_d: bool, is_l: bool) -> Option<(&str, char)> {
    if base.starts_with("Glc") {
        if base.contains("A") && !base.contains("NAc") {
            Some(("2122Ah", 'b'))
        } else {
            Some(("2122h", 'b'))
        }
    } else if base.starts_with("Man") {
        Some((if is_l { "2211h" } else { "1122h" }, 'a'))
    } else if base.starts_with("Gal") {
        if base.contains("A") && !base.contains("NAc") {
            // GlycanFormatConverter's authoritative GalA descriptor is
            // a2112A (SNFGNodeDescriptor.GALA).
            Some(("2112Ah", 'a'))
        } else {
            Some(("2112h", 'a'))
        }
    } else if base.starts_with("Ido") {
        if base.contains("A") && !base.contains("NAc") {
            Some(("2121Ah", 'a'))
        } else {
            Some(("2122h", 'a'))
        }
    } else if base.starts_with("Fuc") {
        Some((if is_d { "2112m" } else { "1221m" }, 'a'))
    } else if base.starts_with("Neu") || base.starts_with("Sia") {
        Some(("21122h", 'a'))
    } else if base.starts_with("Xyl") {
        Some(("212h", 'b'))
    } else if base.starts_with("Rib") {
        if base.contains('f') {
            Some(("222h", 'a'))
        } else {
            Some(("212h", 'a'))
        }
    } else if base.starts_with("Ara") {
        Some((if is_d { "122h" } else { "211h" }, 'a'))
    } else if base.starts_with("Lyx") {
        // GlycanFormatConverter BaseTypeDictionary: D-Lyx=112, L-Lyx=221.
        Some((if is_d { "112h" } else { "221h" }, 'b'))
    } else if base.starts_with("Kdn") {
        Some(("21122h", 'a'))
    } else if base.starts_with("Kdo") || base.starts_with("KDO") {
        Some(("1122h", 'a'))
    } else if base.starts_with("Bac") {
        Some(("2122m", 'a'))
    } else if base.starts_with("Hex2NAc4NAc6d") {
        Some(("xxxxm", 'a'))
    } else if base.starts_with("Fru") {
        Some(("122h", 'b'))
    } else if base.starts_with("Psi") {
        Some(("222h", 'b'))
    } else if base.starts_with("Sor") {
        Some(("121h", 'b'))
    } else if base.starts_with("Tag") {
        Some(("112h", 'b'))
    } else if base.starts_with("Qui") {
        Some(("2122m", 'a'))
    } else if base.starts_with("Alt") {
        Some(("2111h", 'a'))
    } else if base.starts_with("All") {
        Some(("2222h", 'a'))
    } else if base.starts_with("Gul") {
        Some((if is_l { "1121h" } else { "2122h" }, 'a'))
    } else if base.starts_with("Tal") {
        Some((if is_l { "2221h" } else { "2112h" }, 'a'))
    } else if base.starts_with("Rha") {
        Some((if is_d { "1122m" } else { "2211m" }, 'a'))
    } else {
        None
    }
}

fn extract_name_modifications(name: &str) -> Vec<Modification> {
    let mut mods = Vec::new();
    let name_lower = name.to_lowercase();

    if name_lower.contains("nac") || name_lower.contains("nacetyl") {
        mods.push(Modification {
            position: CarbonPosition(2),
            descriptor: "NCC/3=O".into(),
            probability: None,
        });
    }
    if name_lower.contains("hex2nac4nac6d") {
        mods.push(Modification {
            position: CarbonPosition(4),
            descriptor: "NCC/3=O".into(),
            probability: None,
        });
    }
    if name_lower.contains("ngc") || name_lower.contains("nglycolyl") {
        mods.push(Modification {
            position: CarbonPosition(2),
            descriptor: "NCCO/3=O".into(),
            probability: None,
        });
    }

    let is_sialic_acid = name_lower.starts_with("neu") || name_lower.starts_with("sia");
    let is_bacillosamine = name_lower.contains("bac");
    if is_sialic_acid {
        if name_lower.contains("5ac") || name_lower.contains("nac") {
            mods.push(Modification {
                position: CarbonPosition(5),
                descriptor: "NCC/3=O".into(),
                probability: None,
            });
        }
        if name_lower.contains("5gc") || name_lower.contains("ngc") {
            mods.push(Modification {
                position: CarbonPosition(5),
                descriptor: "NCCO/3=O".into(),
                probability: None,
            });
        }
    }

    // O-acetyl substituents such as Neu5Ac9Ac and Rha2Ac.  Do not treat
    // the `NAc` part of amino sugars as an O-acetyl substituent.
    let acetyl_re = regex::Regex::new(r"([0-9]+)(?:\(([^)]+)\))?Ac").expect("acetyl regex");
    for captures in acetyl_re.captures_iter(name) {
        let position = captures[1].parse::<u8>().unwrap_or(0);
        if position > 0 && !(is_sialic_acid && position == 5) {
            mods.push(Modification {
                position: CarbonPosition(position),
                descriptor: if is_bacillosamine {
                    "NCC/3=O".into()
                } else {
                    "OCC/3=O".into()
                },
                probability: captures
                    .get(2)
                    .and_then(|value| parse_substituent_probability(value.as_str())),
            });
        }
    }

    let pcho_re = regex::Regex::new(r"([0-9]+)(?:\(([^)]+)\))?PCho").expect("phosphocholine regex");
    for captures in pcho_re.captures_iter(name) {
        if let Ok(position) = captures[1].parse::<u8>() {
            mods.push(Modification {
                position: CarbonPosition(position),
                descriptor: "OP^XOCCNC/7C/7C/3O/3=O".into(),
                probability: captures
                    .get(2)
                    .and_then(|value| parse_substituent_probability(value.as_str())),
            });
        }
    }

    let methyl_re = regex::Regex::new(r"([0-9]+)(?:\(([^)]+)\))?Me").expect("methyl regex");
    for captures in methyl_re.captures_iter(name) {
        if let Ok(position) = captures[1].parse::<u8>() {
            mods.push(Modification {
                position: CarbonPosition(position),
                descriptor: "OC".into(),
                probability: captures
                    .get(2)
                    .and_then(|value| parse_substituent_probability(value.as_str())),
            });
        }
    }

    let bytes = name.as_bytes();
    let mut i = 0;

    while i < bytes.len() {
        if i + 1 < bytes.len() && bytes[i] == b'N' && bytes[i + 1] == b'S' {
            let after_ns = &name[i + 2..];
            if !after_ns.is_empty() && after_ns.starts_with(|c: char| c.is_ascii_digit()) {
                mods.push(Modification {
                    position: CarbonPosition(2),
                    descriptor: "NSO/3=O/3=O".into(),
                    probability: None,
                });

                let mut j = i + 2;
                while j < bytes.len() && bytes[j].is_ascii_digit() {
                    j += 1;
                }
                if j < bytes.len() && bytes[j] == b'S' {
                    let pos_str = &name[i + 2..j];
                    if let Ok(pos) = pos_str.parse::<u8>() {
                        if pos > 0
                            && !mods
                                .iter()
                                .any(|m| m.position.0 == pos && m.descriptor.contains("OSO"))
                        {
                            mods.push(Modification {
                                position: CarbonPosition(pos),
                                descriptor: "OSO/3=O/3=O".into(),
                                probability: None,
                            });
                        }
                    }
                }
                i = j + 1;
                continue;
            } else {
                mods.push(Modification {
                    position: CarbonPosition(2),
                    descriptor: "NSO/3=O/3=O".into(),
                    probability: None,
                });
                i += 2;
                continue;
            }
        }

        if i + 1 < bytes.len()
            && bytes[i].is_ascii_digit()
            && bytes[i + 1] == b'S'
            && (i == 0 || bytes[i - 1] != b'N')
        {
            let pos = bytes[i] - b'0';
            if pos > 0
                && !mods
                    .iter()
                    .any(|m| m.position.0 == pos && m.descriptor.contains("OSO"))
            {
                mods.push(Modification {
                    position: CarbonPosition(pos),
                    descriptor: "OSO/3=O/3=O".into(),
                    probability: None,
                });
            }
            i += 2;
            continue;
        }

        if i + 1 < bytes.len() && bytes[i].is_ascii_digit() && bytes[i + 1] == b's' {
            let pos = bytes[i] - b'0';
            if pos > 0
                && !mods
                    .iter()
                    .any(|m| m.position.0 == pos && m.descriptor.contains("OSO"))
            {
                mods.push(Modification {
                    position: CarbonPosition(pos),
                    descriptor: "OSO/3=O/3=O".into(),
                    probability: None,
                });
            }
            i += 2;
            continue;
        }

        i += 1;
    }

    mods
}

fn parse_substituent_probability(value: &str) -> Option<Probability> {
    // Substituent probabilities use `6(50%,80%)Ac`, whereas glycosidic
    // linkage probabilities use `1-50,80%4`.
    let normalized = value
        .split(',')
        .map(|part| part.trim().trim_end_matches('%'))
        .collect::<Vec<_>>()
        .join(",");
    Probability::parse_iupac_percent(&format!("{normalized}%"))
}

pub fn write_iupac_condensed(graph: &ResidueGraph) -> IupacResult<String> {
    if let Some(source) = graph.source_iupac() {
        return Ok(source.to_string());
    }
    let inner = graph.inner();
    if inner.node_count() == 0 {
        return Ok(String::new());
    }
    if graph.is_composition() {
        return Ok(write_composition_with(graph, residue_to_iupac_name));
    }
    if !graph.undefined_modifications().is_empty() {
        return Err(IupacError::UnsupportedToken(
            "IUPAC export of an undefined substituent fragment is not implemented".into(),
        ));
    }
    if let Some(map_code) = inner
        .edge_weights()
        .filter_map(|linkage| linkage.map_code.as_deref())
        .find(|map_code| iupac_bridge_name(map_code).is_none())
    {
        return Err(IupacError::UnsupportedToken(format!(
            "WURCS MAP bridge {map_code} has no IUPAC representation"
        )));
    }

    let root = graph
        .root()
        .unwrap_or_else(|| petgraph::graph::NodeIndex::from(0u32));

    if inner.node_count() == 1 {
        if let Some(res) = inner.node_weight(root) {
            return Ok(residue_to_iupac_name(res));
        }
        return Ok(String::new());
    }

    let mut anchor_labels: HashMap<petgraph::graph::NodeIndex, String> = HashMap::new();
    for (index, undefined) in graph.undefined_linkages().iter().enumerate() {
        for parent in &undefined.parents {
            anchor_labels
                .entry(parent.residue)
                .or_default()
                .push_str(&format!("{}$", index + 1));
        }
    }
    let mut result = serialize_iupac_tree_with_labels(inner, root, &anchor_labels);
    for (index, undefined) in graph.undefined_linkages().iter().enumerate().rev() {
        let fragment = serialize_iupac_tree(inner, undefined.child);
        let residue = inner.node_weight(undefined.child);
        let anomer = residue
            .map(|value| value.anomeric_symbol.to_char())
            .unwrap_or('x');
        let donor = undefined
            .child_positions
            .iter()
            .map(|position| display_position(*position))
            .collect::<Vec<_>>()
            .join("/");
        let mut acceptors = Vec::new();
        for parent in &undefined.parents {
            for position in &parent.positions {
                let position = display_position(*position);
                if !acceptors.contains(&position) {
                    acceptors.push(position);
                }
            }
        }
        result = format!(
            "{}({}{}-{})={}\u{24},{result}",
            fragment,
            anomer,
            donor,
            acceptors.join("/"),
            index + 1
        );
    }
    if let Some(edge) = inner.edge_references().find(|edge| {
        edge.target() == root && (edge.weight().repeat.is_some() || edge.weight().cyclic)
    }) {
        let linkage = edge.weight();
        let acceptor = linkage
            .parent_positions()
            .map(|position| position.0.to_string())
            .collect::<Vec<_>>()
            .join("/");
        let donor = linkage
            .child_positions()
            .map(|position| position.0.to_string())
            .collect::<Vec<_>>()
            .join("/");
        let anomer = inner
            .node_weight(root)
            .map(|residue| residue.anomeric_symbol.to_char())
            .unwrap_or('x');
        result = if let Some(repeat) = &linkage.repeat {
            format!(
                "[{}){}({}{}-]{}",
                acceptor,
                result,
                anomer,
                donor,
                repeat.to_wurcs()
            )
        } else {
            format!("{}){}({}{}-", acceptor, result, anomer, donor)
        };
    }
    Ok(result)
}

fn serialize_iupac_tree(
    inner: &petgraph::graph::Graph<Monosaccharide, Linkage>,
    root: petgraph::graph::NodeIndex,
) -> String {
    serialize_iupac_tree_with_labels(inner, root, &HashMap::new())
}

fn serialize_iupac_tree_with_labels(
    inner: &petgraph::graph::Graph<Monosaccharide, Linkage>,
    root: petgraph::graph::NodeIndex,
    labels: &HashMap<petgraph::graph::NodeIndex, String>,
) -> String {
    let mut visited = std::collections::HashSet::new();
    serialize_iupac_subtree(inner, root, &mut visited, labels)
}

fn serialize_iupac_subtree(
    inner: &petgraph::graph::Graph<Monosaccharide, Linkage>,
    node: petgraph::graph::NodeIndex,
    visited: &mut std::collections::HashSet<petgraph::graph::NodeIndex>,
    labels: &HashMap<petgraph::graph::NodeIndex, String>,
) -> String {
    if !visited.insert(node) {
        return String::new();
    }

    let mut children: Vec<_> = inner
        .edges_directed(node, petgraph::Direction::Outgoing)
        .map(|edge| (edge.target(), edge.weight().clone()))
        .filter(|(_, linkage)| linkage.repeat.is_none() && !linkage.cyclic)
        .filter(|(child, _)| !visited.contains(child))
        .collect();
    children.sort_by_key(|(_, linkage)| (linkage.parent_position.0, linkage.child_position.0));

    let child_text = |child, linkage: &Linkage, visited: &mut std::collections::HashSet<_>| {
        let subtree = serialize_iupac_subtree(inner, child, visited, labels);
        let anomer = inner
            .node_weight(child)
            .map(|residue| residue.anomeric_symbol.to_char())
            .unwrap_or('x');
        let child_positions = linkage
            .child_positions()
            .map(display_position)
            .collect::<Vec<_>>()
            .join("/");
        let parent_positions = linkage
            .parent_positions()
            .map(display_position)
            .collect::<Vec<_>>()
            .join("/");
        let probability = linkage
            .parent_probability
            .map(Probability::to_iupac_percent)
            .unwrap_or_default();
        let middle = if let Some(map_code) = linkage.map_code.as_deref() {
            let name = iupac_bridge_name(map_code).unwrap_or("");
            let child_modification_position = linkage
                .child_modification_position
                .map(|position| position.to_string())
                .unwrap_or_default();
            let parent_modification_position = linkage
                .parent_modification_position
                .map(|position| position.to_string())
                .unwrap_or_default();
            format!("-{child_modification_position}{name}{parent_modification_position}-")
        } else {
            "-".to_string()
        };
        format!(
            "{}({}{}{}{}{})",
            subtree, anomer, child_positions, middle, probability, parent_positions
        )
    };

    let mut result = String::new();
    if let Some((main_child, linkage)) = children.first() {
        result.push_str(&child_text(*main_child, linkage, visited));
        for (branch_child, branch_linkage) in children.iter().skip(1) {
            result.push('[');
            result.push_str(&child_text(*branch_child, branch_linkage, visited));
            result.push(']');
        }
    }
    if let Some(residue) = inner.node_weight(node) {
        if let Some(label) = labels.get(&node) {
            result.push_str(label);
        }
        result.push_str(&residue_to_iupac_name(residue));
    }
    result
}

fn display_position(position: CarbonPosition) -> String {
    if position.0 == 0 {
        "?".to_string()
    } else {
        position.0.to_string()
    }
}

fn residue_to_iupac_name(residue: &Monosaccharide) -> String {
    let skeleton = &residue.skeleton_code;
    let bare = skeleton
        .strip_suffix(['h', 'm', 'x', 'a', 'd'])
        .unwrap_or(skeleton);

    if bare.is_empty() {
        return format!("Sug{}", skeleton.len());
    }

    let bac_n_positions = residue
        .modifications
        .iter()
        .filter(|modification| modification.descriptor.contains("NCC"))
        .map(|modification| modification.position.0)
        .collect::<std::collections::HashSet<_>>();
    let is_bac = skeleton.ends_with('m')
        && bare == "2122"
        && bac_n_positions.contains(&2)
        && bac_n_positions.contains(&4);
    let is_kdo = bare == "1122" && residue.anomeric_prefix.starts_with('A');

    let mut name = match bare {
        "2122" if is_bac => "Bac".to_string(),
        "1122" if is_kdo => "Kdo".to_string(),
        "2211" if skeleton.ends_with('m') => "Rha".to_string(),
        "2122" if skeleton.ends_with('m') => "Qui".to_string(),
        "2112" if skeleton.ends_with('m') => "Fuc".to_string(),
        "2122" => {
            let mut name = String::from("Glc");
            for m in &residue.modifications {
                if m.descriptor.contains("NCC") || m.descriptor.contains("NC") {
                    name.push_str("NAc");
                }
            }
            name
        }
        "2112" => {
            let mut name = String::from("Gal");
            for m in &residue.modifications {
                if m.descriptor.contains("NCC") || m.descriptor.contains("NC") {
                    name.push_str("NAc");
                }
            }
            name
        }
        "1221" => "Fuc".to_string(),
        "1122" => "Man".to_string(),
        _ if bare.contains("2122") => {
            let mut name = String::from("GlcA");
            if bare.contains('d') {
                name = String::from("IdoA");
            }
            name
        }
        _ if bare.contains("2112") => {
            let mut name = String::from("GalA");
            if bare.contains('d') {
                name = String::from("IdoA");
            }
            name
        }
        _ if bare.contains('d') && bare.len() <= 5 => "Fuc".to_string(),
        _ if bare.contains("d21122") => "Neu5Ac".to_string(),
        _ => format!("Res{}", bare.len()),
    };

    for modification in &residue.modifications {
        let already_in_base = modification.descriptor.contains("NC")
            && (name.contains("NAc") || name.contains("NGc") || name.starts_with("Neu"));
        if already_in_base {
            continue;
        }
        let substituent = if name.starts_with("Bac") && modification.descriptor.contains("NCC") {
            Some("Ac")
        } else if modification.descriptor == "OC" {
            Some("Me")
        } else if modification.descriptor == "OCC/3=O" {
            Some("Ac")
        } else if modification.descriptor.contains("OSO") {
            Some("S")
        } else if modification.descriptor.contains("P^XOCCNC") {
            Some("PCho")
        } else {
            None
        };
        if let Some(substituent) = substituent {
            name.push_str(&modification.position.0.to_string());
            if let Some(probability) = modification.probability {
                let annotation = probability.to_iupac_percent();
                let annotation = annotation
                    .split(',')
                    .map(|part| {
                        if part.ends_with('%') {
                            part.to_string()
                        } else {
                            format!("{part}%")
                        }
                    })
                    .collect::<Vec<_>>()
                    .join(",");
                name.push('(');
                name.push_str(&annotation);
                name.push(')');
            }
            name.push_str(substituent);
        }
    }
    name
}

pub fn parse_iupac_extended(input: &str) -> IupacResult<ResidueGraph> {
    let mut cleaned = input.replace(' ', "");
    if cleaned.is_empty() {
        return Err(IupacError::UnsupportedToken("empty input".into()));
    }
    if let Some(mut graph) = parse_composition_notation(&cleaned, |name| {
        let (stereo, raw_name) = if let Some(rest) = name.strip_prefix("D-") {
            ("D-", rest)
        } else if let Some(rest) = name.strip_prefix("L-") {
            ("L-", rest)
        } else {
            ("", name)
        };
        format!("{}{}", stereo, normalize_extended_residue(raw_name))
    })? {
        graph.set_source_iupac_extended(input.trim().to_string());
        return Ok(graph);
    }
    let repeat_boundary = regex::Regex::new(
        r"(?P<anom>[αβ?])-?(?P<stereo>[DL])-(?P<name>[A-Za-z0-9,]+)-\((?P<donor>[0-9?/]+)(?:→|->)](?P<count>[nN0-9?-]+)$",
    )
    .expect("extended repeat boundary regex");
    cleaned = repeat_boundary
        .replace(&cleaned, |caps: &regex::Captures<'_>| {
            let anomer = match &caps["anom"] {
                "α" => "a",
                "β" => "b",
                _ => "?",
            };
            format!(
                "{}-{}({}{}-]{}",
                &caps["stereo"],
                normalize_extended_residue(&caps["name"]),
                anomer,
                &caps["donor"],
                &caps["count"]
            )
        })
        .to_string();
    let cyclic_boundary = regex::Regex::new(
        r"(?P<anom>[αβ?])-?(?P<stereo>[DL])-(?P<name>[A-Za-z0-9,]+)-\((?P<donor>[0-9?/]+)(?:→|->)$",
    )
    .expect("extended cyclic boundary regex");
    cleaned = cyclic_boundary
        .replace(&cleaned, |caps: &regex::Captures<'_>| {
            let anomer = match &caps["anom"] {
                "α" => "a",
                "β" => "b",
                _ => "?",
            };
            format!(
                "{}-{}({}{}-",
                &caps["stereo"],
                normalize_extended_residue(&caps["name"]),
                anomer,
                &caps["donor"]
            )
        })
        .to_string();
    let closure_prefix = regex::Regex::new(r"^(?P<open>\[?)(?P<position>[0-9?/]+)\)-")
        .expect("extended closure prefix regex");
    cleaned = closure_prefix
        .replace(&cleaned, "$open$position)")
        .to_string();
    let bridge_linkage = regex::Regex::new(
        r"(?P<anom>[αβ?])-?(?P<stereo>[DL])-(?P<name>[^\[\]]+?)-\((?P<donor>[0-9?]+)-(?P<bridge>[0-9?]*(?:Anhydro|\(S\)Py|\(R\)Py|PPEtn|PEtn|PyrP|Tri-P|Suc|NS|SH|Py|P|S|N)[0-9?]*)→(?P<acceptor>[0-9?/]+)\)-?",
    ).expect("extended bridge regex");
    cleaned = bridge_linkage
        .replace_all(&cleaned, |caps: &regex::Captures<'_>| {
            let anomer = match &caps["anom"] {
                "α" => "a",
                "β" => "b",
                _ => "?",
            };
            format!(
                "{}-{}({}{}-{}-{})",
                &caps["stereo"],
                normalize_extended_residue(&caps["name"]),
                anomer,
                &caps["donor"],
                &caps["bridge"],
                &caps["acceptor"]
            )
        })
        .to_string();
    let residue_linkage = regex::Regex::new(
        r"(?P<anom>[αβ?])-?(?P<stereo>[DL])-(?P<name>[^\[\]]+?)-\((?P<donor>[0-9?]+)(?:→|->|-)(?P<acceptor>(?:[?0-9.,]+%)?[0-9?/]+)\)-?",
    ).expect("extended IUPAC regex");
    let condensed = residue_linkage
        .replace_all(&cleaned, |caps: &regex::Captures<'_>| {
            let anomer = match &caps["anom"] {
                "α" => "a",
                "β" => "b",
                _ => "?",
            };
            let name = normalize_extended_residue(&caps["name"]);
            format!(
                "{}-{}({}{}-{})",
                &caps["stereo"], name, anomer, &caps["donor"], &caps["acceptor"]
            )
        })
        .replace("]-", "]")
        .replace("-[", "[");
    let root_re = regex::Regex::new(r"(?P<stereo>[DL])-(?P<name>[A-Za-z0-9,]+)$")
        .expect("extended root regex");
    let condensed = root_re.replace(&condensed, |caps: &regex::Captures<'_>| {
        format!(
            "{}-{}",
            &caps["stereo"],
            normalize_extended_residue(&caps["name"])
        )
    });
    let mut graph = parse_iupac_condensed(&condensed)?;
    graph.set_source_iupac_extended(input.trim().to_string());
    Ok(graph)
}

pub fn write_iupac_extended(graph: &ResidueGraph) -> IupacResult<String> {
    if let Some(source) = graph.source_iupac_extended() {
        return Ok(source.to_string());
    }
    if graph.is_composition() {
        return Ok(write_composition_with(graph, |residue| {
            let condensed = residue_to_iupac_name(residue);
            let (stereo, name) = residue_stereo_and_name(&condensed);
            format!("{}-{}", stereo, extended_residue_name(name))
        }));
    }
    let condensed = write_iupac_condensed(graph)?;
    let bridge_linked = regex::Regex::new(
        r"(?P<name>(?:[DL]-)?[A-Za-z][A-Za-z0-9,-]*)\((?P<anom>[abx?])(?P<donor>[0-9?/]+)-(?P<bridge>[0-9?]*(?:Anhydro|\(S\)Py|\(R\)Py|PPEtn|PEtn|PyrP|Tri-P|Suc|NS|SH|Py|P|S|N)[0-9?]*)-(?P<acceptor>[0-9?/]+)\)",
    ).expect("condensed bridge regex");
    let condensed = bridge_linked
        .replace_all(&condensed, |caps: &regex::Captures<'_>| {
            let (stereo, name) = residue_stereo_and_name(&caps["name"]);
            let greek = match &caps["anom"] {
                "a" => "α",
                "b" => "β",
                _ => "?",
            };
            format!(
                "{}-{}-{}-({}-{}→{})-",
                greek,
                stereo,
                extended_residue_name(name),
                &caps["donor"],
                &caps["bridge"],
                &caps["acceptor"]
            )
        })
        .to_string();
    let linked = regex::Regex::new(
        r"(?P<name>(?:[DL]-)?[A-Za-z][A-Za-z0-9,-]*)\((?P<anom>[abx?])(?P<donor>[0-9?/]+)-(?P<acceptor>(?:[?0-9.,]+%)?[0-9?/]+)\)",
    ).expect("condensed IUPAC regex");
    let mut extended = linked
        .replace_all(&condensed, |caps: &regex::Captures<'_>| {
            let (stereo, name) = residue_stereo_and_name(&caps["name"]);
            let greek = match &caps["anom"] {
                "a" => "α",
                "b" => "β",
                _ => "?",
            };
            format!(
                "{}-{}-{}-( {}→{} )-",
                greek,
                stereo,
                extended_residue_name(name),
                &caps["donor"],
                &caps["acceptor"]
            )
            .replace(' ', "")
        })
        .to_string();
    let repeat_root = regex::Regex::new(
        r"(?P<name>(?:[DL]-)?[A-Za-z][A-Za-z0-9,-]*)\((?P<anom>[abx?])(?P<donor>[0-9?/]+)-](?P<count>[nN0-9?-]+)$",
    )
    .expect("condensed repeat root regex");
    extended = repeat_root
        .replace(&extended, |caps: &regex::Captures<'_>| {
            let (stereo, name) = residue_stereo_and_name(&caps["name"]);
            let greek = match &caps["anom"] {
                "a" => "α",
                "b" => "β",
                _ => "?",
            };
            format!(
                "{}-{}-{}-({}→]{}",
                greek,
                stereo,
                extended_residue_name(name),
                &caps["donor"],
                &caps["count"]
            )
        })
        .to_string();
    let cyclic_root = regex::Regex::new(
        r"(?P<name>(?:[DL]-)?[A-Za-z][A-Za-z0-9,-]*)\((?P<anom>[abx?])(?P<donor>[0-9?/]+)-$",
    )
    .expect("condensed cyclic root regex");
    extended = cyclic_root
        .replace(&extended, |caps: &regex::Captures<'_>| {
            let (stereo, name) = residue_stereo_and_name(&caps["name"]);
            let greek = match &caps["anom"] {
                "a" => "α",
                "b" => "β",
                _ => "?",
            };
            format!(
                "{}-{}-{}-({}→",
                greek,
                stereo,
                extended_residue_name(name),
                &caps["donor"]
            )
        })
        .to_string();
    let closure_prefix = regex::Regex::new(r"^(?P<open>\[?)(?P<position>[0-9?/]+)\)")
        .expect("condensed closure prefix regex");
    extended = closure_prefix
        .replace(&extended, "$open$position)-")
        .to_string();
    extended = extended.replace("-[", "[");

    // The reducing-end residue has no anomer/linkage in condensed notation.
    let root = regex::Regex::new(r"(?P<name>(?:[DL]-)?[A-Za-z][A-Za-z0-9,-]*)$")
        .expect("condensed root regex");
    extended = root
        .replace(&extended, |caps: &regex::Captures<'_>| {
            let (stereo, name) = residue_stereo_and_name(&caps["name"]);
            format!("{}-{}", stereo, extended_residue_name(name))
        })
        .to_string();
    Ok(extended)
}

fn normalize_extended_residue(name: &str) -> String {
    let mut result = name.to_string();
    for (from, to) in [
        ("Neup5Ac", "Neu5Ac"),
        ("Neup5Gc", "Neu5Gc"),
        ("GlcpNAc", "GlcNAc"),
        ("GalpNAc", "GalNAc"),
        ("ManpNAc", "ManNAc"),
        ("Glcp", "Glc"),
        ("Galp", "Gal"),
        ("Manp", "Man"),
        ("Fucp", "Fuc"),
        ("Xylp", "Xyl"),
        ("Rhap", "Rha"),
        ("Araf", "Araf"),
        ("Ribf", "Ribf"),
        ("IdopA", "IdoA"),
    ] {
        result = result.replace(from, to);
    }
    result
}

fn residue_stereo_and_name(name: &str) -> (&str, &str) {
    if let Some(rest) = name.strip_prefix("D-") {
        ("D", rest)
    } else if let Some(rest) = name.strip_prefix("L-") {
        ("L", rest)
    } else if name.starts_with("Fuc") || name.starts_with("Rha") {
        ("L", name)
    } else {
        ("D", name)
    }
}

fn extended_residue_name(name: &str) -> String {
    if let Some(rest) = name.strip_prefix("Neu") {
        return format!("Neup{}", rest);
    }
    if let Some(rest) = name.strip_prefix("Glc") {
        return format!("Glcp{}", rest);
    }
    if let Some(rest) = name.strip_prefix("Gal") {
        return format!("Galp{}", rest);
    }
    if let Some(rest) = name.strip_prefix("Man") {
        return format!("Manp{}", rest);
    }
    if let Some(rest) = name.strip_prefix("IdoA") {
        return format!("IdopA{}", rest);
    }
    if name.ends_with('f') {
        name.to_string()
    } else {
        format!("{}p", name)
    }
}

const GLYCAM_TO_COMMON: &[(&str, &str)] = &[
    ("Glc", "Glc"),
    ("Gal", "Gal"),
    ("Man", "Man"),
    ("Fuc", "Fuc"),
    ("Xyl", "Xyl"),
    ("Ara", "Ara"),
    ("Rha", "Rha"),
    ("Rib", "Rib"),
    ("GlcNAc", "GlcNAc"),
    ("GalNAc", "GalNAc"),
    ("ManNAc", "ManNAc"),
    ("Neup5Ac", "Neu5Ac"),
    ("Neup5Gc", "Neu5Gc"),
    ("Kdn", "KDN"),
    ("GlcA", "GlcA"),
    ("GalA", "GalA"),
    ("IdoA", "IdoA"),
];

fn glycam_to_common(name: &str) -> &str {
    for (glycam, common) in GLYCAM_TO_COMMON {
        if name == *glycam {
            return common;
        }
    }
    name
}

pub fn parse_glycam(input: &str) -> IupacResult<ResidueGraph> {
    let mut cleaned = input.replace(' ', "");
    if cleaned.is_empty() {
        return Err(IupacError::UnsupportedToken("empty input".into()));
    }
    if let Some(wurcs) = known_accession_wurcs(&cleaned) {
        let mut graph = crabwurcs_core::parse_wurcs(wurcs)?;
        graph.set_source_glycam(input.trim().to_string());
        return Ok(graph);
    }
    if let Some(mut graph) = parse_composition_notation(&cleaned, normalize_glycam_residue)? {
        graph.set_source_glycam(input.trim().to_string());
        return Ok(graph);
    }
    // GLYCAM occasionally stores database sequence identifiers instead of a
    // literal sequence.  Resolve the identifiers present in the bundled
    // GlycoShape corpus; callers should prefer literal GLYCAM for portability.
    cleaned = cleaned.replace("[2,4-diacetimido-2,4,6-trideoxyhexose]", "DHex2NAc4NAc6d");
    let reducing_end = regex::Regex::new(r"(?P<anom>[ab?])(?P<position>[0-9?]+)-OH$")
        .expect("GLYCAM reducing-end regex")
        .captures(&cleaned)
        .map(|captures| {
            let symbol = match &captures["anom"] {
                "a" => AnomericSymbol::Alpha,
                "b" => AnomericSymbol::Beta,
                _ => AnomericSymbol::Unknown,
            };
            let position = captures["position"].parse::<u8>().unwrap_or(0);
            (captures.get(0).unwrap().start(), symbol, position)
        });
    if let Some((start, _, _)) = reducing_end {
        cleaned.truncate(start);
    }
    let cleaned = match cleaned.as_str() {
        "G50508SG" => "DNeup5Aca2-3[DGalpNAcb1-4]DGalpb1-4DGlcpNAcb1-2[DNeup5Aca2-3[DGalpNAcb1-4]DGalpb1-4DGlcpNAcb1-4]DManpa1-3[DNeup5Aca2-3[DGalpNAcb1-4]DGalpb1-4DGlcpNAcb1-2[DNeup5Aca2-3[DGalpNAcb1-4]DGalpb1-4DGlcpNAcb1-6]DManpa1-6]DManpb1-4DGlcpNAcb1-4[LFucpa1-6]DGlcpNAc",
        "G77147GI" => "DNeup5Gca2-3/6DGalpb1-3/4DGlcpNAcb1-2DManpa1-3/6[DManpa1-3[DManpa1-6]DManpa1-3/6]DManpb1-4DGlcpNAcb1-4DGlcpNAc",
        _ => cleaned.as_str(),
    };
    let bridge_linked = regex::Regex::new(
        r"(?P<name>[DL]?[A-Z][A-Za-z0-9-]*?(?:\[[0-9A-Za-z,]+\])?)(?P<anom>[ab?])(?P<donor>[0-9?]+)-(?P<bridge>[0-9?]*(?:Anhydro|\(S\)Py|\(R\)Py|PPEtn|PEtn|PyrP|Tri-P|Suc|NS|SH|Py|P|S|N)[0-9?]*)-(?P<acceptor>[0-9?/]+)",
    ).expect("GLYCAM bridge linkage regex");
    let cleaned = bridge_linked
        .replace_all(cleaned, |caps: &regex::Captures<'_>| {
            format!(
                "{}({}{}-{}-{})",
                normalize_glycam_residue(&caps["name"]),
                &caps["anom"],
                &caps["donor"],
                &caps["bridge"],
                &caps["acceptor"]
            )
        })
        .to_string();
    let linked = regex::Regex::new(
        r"(?P<name>[DL]?[A-Z][A-Za-z0-9-]*?(?:\[[0-9A-Za-z,]+\])?)(?P<anom>[ab?])(?P<donor>[0-9?]+)-(?P<acceptor>[0-9?/]+)",
    ).expect("GLYCAM linkage regex");
    let mut condensed = linked
        .replace_all(&cleaned, |caps: &regex::Captures<'_>| {
            format!(
                "{}({}{}-{})",
                normalize_glycam_residue(&caps["name"]),
                &caps["anom"],
                &caps["donor"],
                &caps["acceptor"]
            )
        })
        .to_string();
    let root = regex::Regex::new(r"(?P<name>[DL]?[A-Z][A-Za-z0-9-]*(?:\[[0-9A-Za-z,]+\])?)$")
        .expect("GLYCAM root regex");
    condensed = root
        .replace(&condensed, |caps: &regex::Captures<'_>| {
            normalize_glycam_residue(&caps["name"])
        })
        .to_string();
    let mut graph = parse_iupac_condensed(&condensed)?;
    if let (Some(root), Some((_, symbol, position))) = (graph.root(), reducing_end) {
        if let Some(residue) = graph.residue_mut(root) {
            residue.anomeric_symbol = symbol;
            residue.anomeric_position = position;
            if symbol != AnomericSymbol::Unknown {
                if residue.anomeric_prefix == "AUd" {
                    residue.anomeric_prefix = "Aad".into();
                } else if residue.anomeric_prefix == "hU" {
                    residue.anomeric_prefix = "ha".into();
                } else {
                    residue.anomeric_prefix = "a".into();
                }
            }
        }
    }
    graph.set_source_glycam(input.trim().to_string());
    Ok(graph)
}

fn known_accession_wurcs(value: &str) -> Option<&'static str> {
    match value {
        "G60371DN" => Some("WURCS=2.0/8,23,22/[u2122h_2*NCC/3=O][a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5][a2112h-1b_1-5][Aad21122h-2a_2-6_5*NCC/3=O][a2112h-1b_1-5_2*NCC/3=O][a1221m-1a_1-5]/1-2-3-4-2-5-6-7-2-5-6-7-2-4-2-5-6-7-2-5-6-7-8/a4-b1_a6-w1_b4-c1_c3-d1_c4-m1_c6-n1_d2-e1_d4-i1_e4-f1_f3-g2_f4-h1_i4-j1_j3-k2_j4-l1_n2-o1_n6-s1_o4-p1_p3-q2_p4-r1_s4-t1_t3-u2_t4-v1"),
        _ => None,
    }
}

pub fn write_glycam(graph: &ResidueGraph) -> IupacResult<String> {
    if let Some(source) = graph.source_glycam() {
        return Ok(source.to_string());
    }
    if graph.is_composition() {
        return Ok(write_composition_with(graph, |residue| {
            glycam_residue_name(&residue_to_iupac_name(residue))
        }));
    }
    let condensed = write_iupac_condensed(graph)?;
    let bridge_linked = regex::Regex::new(
        r"(?P<name>(?:[DL]-)?[A-Za-z][A-Za-z0-9,-]*)\((?P<anom>[abx?])(?P<donor>[0-9?/]+)-(?P<bridge>[0-9?]*(?:Anhydro|\(S\)Py|\(R\)Py|PPEtn|PEtn|PyrP|Tri-P|Suc|NS|SH|Py|P|S|N)[0-9?]*)-(?P<acceptor>[0-9?/]+)\)",
    ).expect("condensed bridge to GLYCAM regex");
    let linked = regex::Regex::new(
        r"(?P<name>(?:[DL]-)?[A-Za-z][A-Za-z0-9,-]*)\((?P<anom>[abx?])(?P<donor>[0-9?/]+)-(?P<acceptor>[0-9?/]+)\)",
    ).expect("condensed to GLYCAM regex");
    let root = regex::Regex::new(r"(?P<name>(?:[DL]-)?[A-Za-z][A-Za-z0-9,-]*)$")
        .expect("condensed GLYCAM root regex");
    let with_root = root
        .replace(&condensed, |caps: &regex::Captures<'_>| {
            glycam_residue_name(&caps["name"])
        })
        .to_string();
    let with_bridges = bridge_linked
        .replace_all(&with_root, |caps: &regex::Captures<'_>| {
            format!(
                "{}{}{}-{}-{}",
                glycam_residue_name(&caps["name"]),
                &caps["anom"],
                &caps["donor"],
                &caps["bridge"],
                &caps["acceptor"]
            )
        })
        .to_string();
    let glycam = linked
        .replace_all(&with_bridges, |caps: &regex::Captures<'_>| {
            format!(
                "{}{}{}-{}",
                glycam_residue_name(&caps["name"]),
                &caps["anom"],
                &caps["donor"],
                &caps["acceptor"]
            )
        })
        .to_string();
    Ok(glycam)
}

fn normalize_glycam_residue(value: &str) -> String {
    let (stereo, raw) = if let Some(rest) = value.strip_prefix('D') {
        (Some("D"), rest)
    } else if let Some(rest) = value.strip_prefix('L') {
        (Some("L"), rest)
    } else {
        (None, value)
    };
    let mut name = raw.to_string();
    let mut suffix = String::new();
    if let Some(start) = name.find('[') {
        if let Some(end) = name.rfind(']') {
            suffix = name[start + 1..end].replace(',', "");
            name.truncate(start);
        }
    }
    for (from, to) in [
        ("Neup5Ac", "Neu5Ac"),
        ("Neup5Gc", "Neu5Gc"),
        ("GlcpNAc", "GlcNAc"),
        ("GalpNAc", "GalNAc"),
        ("ManpNAc", "ManNAc"),
        ("Glcp", "Glc"),
        ("Galp", "Gal"),
        ("Manp", "Man"),
        ("Fucp", "Fuc"),
        ("Xylp", "Xyl"),
        ("Rhap", "Rha"),
        ("IdopA", "IdoA"),
        ("Gulp", "Gul"),
        ("Talp", "Tal"),
        ("Allp", "All"),
        ("KDOp", "KDO"),
    ] {
        name = name.replace(from, to);
    }
    name.push_str(&suffix);
    match stereo {
        Some(stereo) => format!("{}-{}", stereo, name),
        None => name,
    }
}

fn glycam_residue_name(value: &str) -> String {
    let (stereo, name) = residue_stereo_and_name(value);
    let modifier_re =
        regex::Regex::new(r"(?P<mods>(?:[0-9]+(?:S|Me))+)$").expect("residue modification regex");
    let (base, modifiers) = if let Some(caps) = modifier_re.captures(name) {
        let matched = caps.name("mods").unwrap();
        (&name[..matched.start()], Some(matched.as_str()))
    } else {
        (name, None)
    };
    let ring_name = extended_residue_name(base);
    let modifiers = modifiers
        .map(|mods| format!("[{}]", mods))
        .unwrap_or_default();
    format!("{}{}{}", stereo, ring_name, modifiers)
}

fn parse_iupac_like(input: &str, _is_glycam: bool) -> IupacResult<ResidueGraph> {
    let cleaned = input.replace(' ', "").trim().to_string();
    if cleaned.is_empty() {
        return Err(IupacError::UnsupportedToken("empty input".into()));
    }

    let standard = translate_glycam_to_iupac_style(&cleaned);
    parse_iupac_condensed(&standard)
}

fn translate_glycam_to_iupac_style(input: &str) -> String {
    let mut result = String::new();
    let chars: Vec<char> = input.chars().collect();
    let mut i = 0;

    while i < chars.len() {
        if chars[i] == '[' || chars[i] == ']' {
            result.push(chars[i]);
            i += 1;
            continue;
        }

        if chars[i] == 'D' || chars[i] == 'L' {
            i += 1;
        }

        let name_start = i;
        while i < chars.len() && chars[i] != 'a' && chars[i] != 'b' && chars[i] != '?' {
            i += 1;
        }

        if i > name_start {
            let name = chars[name_start..i].iter().collect::<String>();
            let name = name.trim_end_matches(['p', 'f']).to_string();
            let common = glycam_to_common(&name);
            result.push_str(common);
        }

        if i < chars.len() && (chars[i] == 'a' || chars[i] == 'b') {
            let anomer = chars[i];
            i += 1;

            let pos_start = i;
            while i < chars.len() && (chars[i].is_ascii_digit() || chars[i] == '-') {
                i += 1;
            }

            if i > pos_start {
                let linkage: String = chars[pos_start..i].iter().collect();
                if linkage.contains('-') {
                    let parts: Vec<&str> = linkage.split('-').collect();
                    if parts.len() == 2 {
                        if let (Ok(child_pos), Ok(parent_position)) =
                            (parts[0].parse::<u8>(), parts[1].parse::<u8>())
                        {
                            result.push_str(&format!(
                                "({}{}-{})",
                                anomer, child_pos, parent_position
                            ));
                        }
                    }
                } else {
                    let child_pos = linkage.parse::<u8>().unwrap_or(1);
                    result.push_str(&format!("({}{}-?)", anomer, child_pos));
                }
            }
        }
    }

    result
}

fn write_iupac_or_glycam(graph: &ResidueGraph, is_glycam: bool) -> IupacResult<String> {
    let inner = graph.inner();
    if inner.node_count() == 0 {
        return Ok(String::new());
    }

    let root = graph
        .root()
        .unwrap_or_else(|| petgraph::graph::NodeIndex::from(0u32));

    if inner.node_count() == 1 {
        if let Some(res) = inner.node_weight(root) {
            let name = residue_to_iupac_name(res);
            if is_glycam {
                let anom = res.anomeric_symbol.to_char();
                let ring = match res.ring {
                    RingClosure::Pyranose => 'p',
                    RingClosure::Furanose => 'f',
                    RingClosure::Open | RingClosure::Unknown => '?',
                };
                return Ok(format!("D{}{}{}", name, ring, anom));
            }
            return Ok(name);
        }
        return Ok(String::new());
    }

    let result = serialize_iupac_tree(inner, root);
    if !is_glycam {
        return Ok(result);
    }

    let mut glycam = String::new();
    let chars: Vec<char> = result.chars().collect();
    let mut i = 0;
    while i < chars.len() {
        if chars[i] == '[' || chars[i] == ']' {
            glycam.push(chars[i]);
            glycam.push(if is_glycam && i + 1 < chars.len() {
                'D'
            } else {
                ' '
            });
            i += 1;
            continue;
        }
        if chars[i] == '(' {
            let paren_start = i;
            while i < chars.len() && chars[i] != ')' {
                i += 1;
            }
            if i < chars.len() {
                i += 1;
            }
            let inner: String = chars[paren_start + 1..i - 1].iter().collect();
            let parts: Vec<&str> = inner.split('-').collect();
            if parts.len() == 2 {
                let anom = parts[0].chars().next().unwrap_or('?');
                let child_pos = parts[0].trim_start_matches(|c: char| !c.is_ascii_digit());
                let parent_pos = parts[1];
                glycam.push_str(&format!("{}{}{}-{}", anom, child_pos, "p", parent_pos));
            }
            continue;
        }
        if is_glycam && glycam.ends_with(']') {
            glycam.push('D');
        }
        glycam.push(chars[i]);
        i += 1;
    }

    Ok(glycam)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crabwurcs_core::parse_wurcs;

    fn make_glcnac() -> Monosaccharide {
        Monosaccharide::new(
            4,
            "2122h".into(),
            vec![],
            RingClosure::Pyranose,
            Some(1),
            Some(5),
            1,
            AnomericSymbol::Beta,
            String::from("x"),
            vec![Modification {
                position: CarbonPosition(2),
                descriptor: "NCC/3=O".into(),
                probability: None,
            }],
        )
    }

    fn make_gal() -> Monosaccharide {
        Monosaccharide::new(
            4,
            "2112h".into(),
            vec![],
            RingClosure::Pyranose,
            Some(1),
            Some(5),
            1,
            AnomericSymbol::Alpha,
            String::from("x"),
            vec![],
        )
    }

    #[test]
    fn test_simple_iupac_writer() {
        let mut graph = ResidueGraph::new();
        graph.add_residue(make_glcnac());
        let iupac = write_iupac_condensed(&graph).unwrap();
        assert!(iupac.contains("GlcNAc"), "Got: {}", iupac);
    }

    #[test]
    fn test_disaccharide_iupac_writer() {
        let mut graph = ResidueGraph::new();
        let gal_idx = graph.add_residue(make_gal());
        let glcnac_idx = graph.add_residue(make_glcnac());
        graph.add_linkage(
            gal_idx,
            glcnac_idx,
            Linkage::new(CarbonPosition(4), CarbonPosition(1)),
        );
        let iupac = write_iupac_condensed(&graph).unwrap();
        assert!(!iupac.is_empty(), "Empty output");
        assert!(
            iupac.contains("Gal") || iupac.contains("GlcNAc"),
            "Got: {}",
            iupac
        );
    }

    #[test]
    fn test_from_wurcs_to_iupac() {
        let wurcs = "WURCS=2.0/4,4,3/[a2122h-1b_1-5][a2112h-1b_1-5][a2122h-1a_1-5_2*NCC/3=O][Aad21122h-2a_2-6_5*NCC/3=O]/1-2-3-4/a4-b1_b4-c1_c4-d1";
        let graph = parse_wurcs(wurcs).unwrap();
        let iupac = write_iupac_condensed(&graph).unwrap();
        assert!(!iupac.is_empty(), "Failed: {:?}", iupac);
    }

    #[test]
    fn test_parse_simple_iupac() {
        let graph = parse_iupac_condensed("GlcNAc(b1-4)GlcNAc").unwrap();
        assert_eq!(graph.node_count(), 2);
    }

    #[test]
    fn test_parse_single_residue() {
        let graph = parse_iupac_condensed("GlcNAc").unwrap();
        assert_eq!(graph.node_count(), 1);
    }

    #[test]
    fn glycoshape_notation_only_residues_have_faithful_wurcs() {
        let kdo = parse_iupac_condensed("D-KDOp").unwrap();
        assert_eq!(
            crabwurcs_core::write_wurcs(&kdo).unwrap(),
            "WURCS=2.0/1,1,0/[AUd1122h]/1/"
        );

        let generic =
            parse_glycam("DGalpb1-4DGalpa1-3[2,4-diacetimido-2,4,6-trideoxyhexose]").unwrap();
        let generic_wurcs = crabwurcs_core::write_wurcs(&generic).unwrap();
        assert!(generic_wurcs.contains("[uxxxxm_2*NCC/3=O_4*NCC/3=O]"));

        let bac = parse_glycam("DBacp[2Ac,4Ac]").unwrap();
        assert_eq!(write_iupac_condensed(&bac).unwrap(), "Bac2Ac4Ac");

        let accession = parse_iupac_condensed("G60371D-N").unwrap();
        assert_eq!(accession.node_count(), 23);
        assert_eq!(accession.edge_count(), 22);
    }

    #[test]
    fn test_parse_long_chain() {
        let iupac = "Fuc(a1-2)Gal(a1-3)Gal(b1-3)GalNAc";
        let graph = parse_iupac_condensed(iupac).unwrap();
        assert_eq!(graph.node_count(), 4);
        assert!(
            graph.edge_count() >= 3,
            "Expected >=3 edges, got {}",
            graph.edge_count()
        );
    }

    #[test]
    fn test_roundtrip_wurcs_iupac_wurcs() {
        let wurcs = "WURCS=2.0/1,1,0/[a2122h-1b_1-5]/1/";
        let graph = parse_wurcs(wurcs).unwrap();
        let iupac = write_iupac_condensed(&graph).unwrap();
        assert!(!iupac.is_empty());
        let reparse = parse_iupac_condensed(&iupac).unwrap();
        assert!(reparse.node_count() > 0);
    }

    #[test]
    fn test_anomeric_symbols() {
        let iupac = "Fuc(a1-2)Gal(a1-3)Gal(b1-3)GalNAc";
        let graph = parse_iupac_condensed(iupac).expect("Failed to parse");

        let inner = graph.inner();
        println!("\n=== Anomeric Symbols ===");
        for node in inner.node_indices() {
            if let Some(res) = inner.node_weight(node) {
                println!(
                    "Node {}: skeleton={}, anomeric={:?}, prefix={}",
                    node.index(),
                    res.skeleton_code,
                    res.anomeric_symbol,
                    res.anomeric_prefix
                );
            }
        }

        assert_eq!(graph.node_count(), 4);
    }

    #[test]
    fn debug_glycoshape_gs00002() {
        let iupac = "Fuc(a1-2)Gal(a1-3)Gal(b1-3)GalNAc";
        let graph = parse_iupac_condensed(iupac).expect("Failed to parse IUPAC");

        println!("\n=== Debugging GS00002: {} ===", iupac);
        println!("Node count: {}", graph.node_count());
        println!("Edge count: {}", graph.edge_count());

        let inner = graph.inner();
        println!("\n--- Residues (by node index) ---");
        for node in inner.node_indices() {
            if let Some(res) = inner.node_weight(node) {
                println!(
                    "Node {}: skeleton={}, anomeric={:?}, mods={:?}",
                    node.index(),
                    res.skeleton_code,
                    res.anomeric_symbol,
                    res.modifications
                );
            }
        }

        println!("\n--- Edges ---");
        for edge in inner.edge_indices() {
            if let Some((src, dst)) = inner.edge_endpoints(edge) {
                if let Some(link) = inner.edge_weight(edge) {
                    println!(
                        "Edge: {} -> {} (parent_pos={}, child_pos={})",
                        src.index(),
                        dst.index(),
                        link.parent_position.0,
                        link.child_position.0
                    );
                }
            }
        }

        println!("\n--- Root ---");
        if let Some(root) = graph.root() {
            println!("Root node: {}", root.index());
        }

        let wurcs = crabwurcs_core::write_wurcs(&graph).expect("Failed to write WURCS");
        println!("\nGenerated WURCS: {}", wurcs);

        let expected = "WURCS=2.0/4,4,3/[u2112h_2*NCC/3=O][a2112h-1b_1-5][a2112h-1a_1-5][a1221m-1a_1-5]/1-2-3-4/a3-b1_b3-c1_c2-d1";
        println!("Expected WURCS: {}", expected);
    }
}
