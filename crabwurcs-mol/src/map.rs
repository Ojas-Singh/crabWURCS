use std::collections::HashMap;

use chematic::core::{Atom, AtomIdx, BondOrder, Chirality, Element, MoleculeBuilder};

use crate::{MolError, MolResult};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MapNode {
    Atom(usize),
    Star(u8),
}

#[derive(Debug, Clone)]
struct MapAtom {
    element: Element,
    aromatic: bool,
    charge: i8,
    chirality: Chirality,
}

#[derive(Debug, Clone, Copy)]
struct MapBond {
    left: usize,
    right: usize,
    order: BondOrder,
}

/// A parsed WURCS MAP graph. Stars are attachment points and are not emitted as
/// atoms; each star maps to the atom immediately bonded to that attachment.
#[derive(Debug, Clone)]
pub(crate) struct MapGraph {
    atoms: Vec<MapAtom>,
    bonds: Vec<MapBond>,
    attachments: HashMap<u8, (usize, BondOrder)>,
}

#[derive(Debug)]
pub(crate) struct BuiltMap {
    pub(crate) attachments: HashMap<u8, (AtomIdx, BondOrder)>,
}

fn error(offset: usize, message: impl Into<String>) -> MolError {
    MolError::InvalidMap {
        offset,
        message: message.into(),
    }
}

fn parse_number(bytes: &[u8], cursor: &mut usize) -> Option<usize> {
    let start = *cursor;
    while bytes.get(*cursor).is_some_and(u8::is_ascii_digit) {
        *cursor += 1;
    }
    (*cursor > start)
        .then(|| {
            std::str::from_utf8(&bytes[start..*cursor])
                .ok()?
                .parse()
                .ok()
        })
        .flatten()
}

fn parse_element(input: &str, cursor: &mut usize) -> MolResult<(Element, bool)> {
    let bytes = input.as_bytes();
    let start = *cursor;
    let first = *bytes
        .get(*cursor)
        .ok_or_else(|| error(*cursor, "expected an atom"))?;
    if !first.is_ascii_alphabetic() {
        return Err(error(*cursor, "expected an atom symbol"));
    }
    *cursor += 1;
    if first.is_ascii_uppercase()
        && bytes
            .get(*cursor)
            .is_some_and(|character| character.is_ascii_lowercase())
    {
        *cursor += 1;
    }
    let raw = &input[start..*cursor];
    let aromatic = first.is_ascii_lowercase();
    let mut symbol = raw.to_string();
    if aromatic {
        let mut chars = symbol.chars();
        symbol = chars
            .next()
            .map(|character| character.to_ascii_uppercase())
            .into_iter()
            .chain(chars)
            .collect();
    }
    let element = Element::from_symbol(&symbol)
        .ok_or_else(|| error(start, format!("unknown element `{raw}`")))?;
    Ok((element, aromatic))
}

fn connect(
    graph: &mut MapGraph,
    left: MapNode,
    right: MapNode,
    order: BondOrder,
    offset: usize,
) -> MolResult<()> {
    match (left, right) {
        (MapNode::Atom(left), MapNode::Atom(right)) => {
            if left == right
                || graph.bonds.iter().any(|bond| {
                    (bond.left == left && bond.right == right)
                        || (bond.left == right && bond.right == left)
                })
            {
                return Err(error(offset, "duplicate or self bond"));
            }
            graph.bonds.push(MapBond { left, right, order });
        }
        (MapNode::Star(star), MapNode::Atom(atom)) | (MapNode::Atom(atom), MapNode::Star(star)) => {
            if graph.attachments.insert(star, (atom, order)).is_some() {
                return Err(error(offset, format!("attachment star {star} is reused")));
            }
        }
        (MapNode::Star(_), MapNode::Star(_)) => {
            return Err(error(
                offset,
                "two attachment stars cannot be directly bonded",
            ));
        }
    }
    Ok(())
}

impl MapGraph {
    /// Parse a complete MAP string. The leading attachment star is required.
    pub(crate) fn parse(input: &str) -> MolResult<Self> {
        let bytes = input.as_bytes();
        let mut cursor = 0usize;
        let mut graph = Self {
            atoms: Vec::new(),
            bonds: Vec::new(),
            attachments: HashMap::new(),
        };
        let mut traversed = Vec::<MapNode>::new();
        let mut current = None::<MapNode>;
        let mut pending_order = BondOrder::Single;
        let mut next_implicit_star = 1u8;

        while cursor < bytes.len() {
            match bytes[cursor] {
                b'/' => {
                    cursor += 1;
                    let branch_offset = cursor;
                    let index = parse_number(bytes, &mut cursor)
                        .ok_or_else(|| error(branch_offset, "branch requires an atom index"))?;
                    current = traversed
                        .get(index.wrapping_sub(1))
                        .copied()
                        .ok_or_else(|| error(branch_offset, "branch atom index is out of range"))?
                        .into();
                }
                b'$' => {
                    cursor += 1;
                    let ring_offset = cursor;
                    let index = parse_number(bytes, &mut cursor)
                        .ok_or_else(|| error(ring_offset, "ring closure requires an atom index"))?;
                    let target = traversed
                        .get(index.wrapping_sub(1))
                        .copied()
                        .ok_or_else(|| error(ring_offset, "ring atom index is out of range"))?;
                    let source = current
                        .ok_or_else(|| error(ring_offset, "ring closure has no source atom"))?;
                    connect(&mut graph, source, target, pending_order, ring_offset)?;
                    pending_order = BondOrder::Single;
                }
                b'=' => {
                    pending_order = BondOrder::Double;
                    cursor += 1;
                }
                b'#' => {
                    pending_order = BondOrder::Triple;
                    cursor += 1;
                }
                b':' => {
                    pending_order = BondOrder::Aromatic;
                    cursor += 1;
                }
                b'^' => {
                    let stereo_offset = cursor;
                    cursor += 1;
                    let descriptor = bytes
                        .get(cursor)
                        .copied()
                        .ok_or_else(|| error(stereo_offset, "missing stereo descriptor"))?;
                    cursor += 1;
                    if let Some(MapNode::Atom(atom)) = current {
                        graph.atoms[atom].chirality = match descriptor {
                            b'R' => Chirality::Clockwise,
                            b'S' => Chirality::CounterClockwise,
                            b'X' | b'E' | b'Z' => Chirality::None,
                            _ => {
                                return Err(error(stereo_offset, "unsupported stereo descriptor"));
                            }
                        };
                    }
                }
                b'+' | b'-' => {
                    let charge_offset = cursor;
                    let sign = if bytes[cursor] == b'+' { 1i8 } else { -1i8 };
                    cursor += 1;
                    let magnitude: i8 = parse_number(bytes, &mut cursor)
                        .unwrap_or(1)
                        .try_into()
                        .map_err(|_| error(charge_offset, "charge magnitude is too large"))?;
                    let Some(MapNode::Atom(atom)) = current else {
                        return Err(error(charge_offset, "charge has no atom"));
                    };
                    graph.atoms[atom].charge = sign * magnitude;
                }
                b'*' => {
                    let node_offset = cursor;
                    cursor += 1;
                    let explicit = parse_number(bytes, &mut cursor);
                    let star = explicit
                        .map(|value| {
                            u8::try_from(value)
                                .map_err(|_| error(node_offset, "star index is too large"))
                        })
                        .transpose()?
                        .unwrap_or_else(|| {
                            let value = next_implicit_star;
                            next_implicit_star = next_implicit_star.saturating_add(1);
                            value
                        });
                    let node = MapNode::Star(star);
                    if let Some(previous) = current {
                        connect(&mut graph, previous, node, pending_order, node_offset)?;
                    }
                    traversed.push(node);
                    current = Some(node);
                    pending_order = BondOrder::Single;
                }
                character if character.is_ascii_alphabetic() => {
                    let node_offset = cursor;
                    let (element, aromatic) = parse_element(input, &mut cursor)?;
                    let atom = graph.atoms.len();
                    graph.atoms.push(MapAtom {
                        element,
                        aromatic,
                        charge: 0,
                        chirality: Chirality::None,
                    });
                    let node = MapNode::Atom(atom);
                    if let Some(previous) = current {
                        connect(&mut graph, previous, node, pending_order, node_offset)?;
                    }
                    traversed.push(node);
                    current = Some(node);
                    pending_order = BondOrder::Single;
                }
                _ => {
                    return Err(error(
                        cursor,
                        format!("unexpected character `{}`", bytes[cursor] as char),
                    ));
                }
            }
        }
        if graph.atoms.is_empty() {
            return Err(error(0, "MAP contains no atoms"));
        }
        if graph.attachments.is_empty() {
            return Err(error(0, "MAP contains no attachment star"));
        }
        Ok(graph)
    }

    pub(crate) fn build(&self, builder: &mut MoleculeBuilder) -> MolResult<BuiltMap> {
        let atoms = self
            .atoms
            .iter()
            .map(|source| {
                let mut atom = Atom::new(source.element);
                atom.aromatic = source.aromatic;
                atom.charge = source.charge;
                atom.chirality = source.chirality;
                builder.add_atom(atom)
            })
            .collect::<Vec<_>>();
        for bond in &self.bonds {
            builder
                .add_bond(atoms[bond.left], atoms[bond.right], bond.order)
                .map_err(|cause| MolError::InvalidValence(cause.to_string()))?;
        }
        Ok(BuiltMap {
            attachments: self
                .attachments
                .iter()
                .map(|(star, (atom, order))| (*star, (atoms[*atom], *order)))
                .collect(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_branched_n_acetyl_map() {
        let map = MapGraph::parse("*NCC/3=O").unwrap();
        assert_eq!(map.atoms.len(), 4);
        assert_eq!(map.bonds.len(), 3);
        assert_eq!(map.attachments.len(), 1);
        assert!(map.bonds.iter().any(|bond| bond.order == BondOrder::Double));
    }

    #[test]
    fn parses_indexed_phosphate_bridge() {
        let map = MapGraph::parse("*1OPO*2/3O/3=O").unwrap();
        assert_eq!(map.attachments.len(), 2);
        assert!(map.attachments.contains_key(&1));
        assert!(map.attachments.contains_key(&2));
    }

    #[test]
    fn reports_map_offsets() {
        assert!(matches!(
            MapGraph::parse("*OC/99N"),
            Err(MolError::InvalidMap { offset: 4, .. })
        ));
    }
}
