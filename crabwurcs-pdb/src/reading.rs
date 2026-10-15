//! Combine what the coordinates say about a monosaccharide with what its
//! residue name says, and record every difference.

use crate::PdbWarningKind;
use crate::geometry::{CarbonReading, DISTORTED_VOLUME};
use crate::names::NameReading;
use crabwurcs_core::{AnomericSymbol, Monosaccharide, ResidueGraph, parse_wurcs};

/// A substituent read from the coordinates.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct CoordinateSubstituent {
    pub position: u8,
    pub descriptor: String,
    /// Its atoms beyond the heteroatom belong to another residue.
    pub external: bool,
    /// A heteroatom with nothing beyond it (`N` for an amino group).
    pub bare: bool,
}

/// Everything the coordinates say about one monosaccharide.
#[derive(Debug, Clone)]
pub(crate) struct CoordinateReading {
    /// Whether the residue's atoms include hydrogens. A model with hydrogens
    /// is taken to be complete: an absent heavy atom is absent.
    pub complete: bool,
    pub anomeric_position: u8,
    pub ring_end: u8,
    pub carbons: Vec<CarbonReading>,
    pub substituents: Vec<CoordinateSubstituent>,
    /// Backbone positions whose heteroatom carries a group this reader does
    /// not recognise.
    pub unrecognized: Vec<u8>,
}

/// Substituents as (position, WURCS MAP descriptor).
type Substituents = Vec<(u8, String)>;

/// One difference between name and coordinates, or one part the
/// coordinates cannot decide.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Note {
    pub kind: PdbWarningKind,
    pub positions: Vec<u8>,
    pub aspect: String,
}

/// The policy-independent result of combining name and coordinates.
pub(crate) struct Combined {
    /// Coordinates decide; the name fills in what they cannot show.
    pub coordinates_first: Monosaccharide,
    /// The name decides where it says something; the coordinates fill in.
    pub names_first: Monosaccharide,
    /// The coordinates alone (`x` where they cannot decide).
    pub coordinates_only: Monosaccharide,
    pub notes: Vec<Note>,
}

fn full_code(residue: &Monosaccharide) -> Vec<char> {
    format!("{}{}", residue.anomeric_prefix, residue.skeleton_code)
        .chars()
        .collect()
}

fn oxygen_count(terminal: char) -> Option<usize> {
    match terminal {
        'm' => Some(0),
        'h' => Some(1),
        'A' => Some(2),
        _ => None,
    }
}

/// The anomeric reference atom: the highest-numbered stereocentre of the
/// first group of (up to four) stereocentres after the anomeric carbon
/// (IUPAC 2-Carb-6.2).
///
/// A centre whose configuration is unknown (`x`) still counts, so that it is
/// never replaced by the centre before it.
fn reference_index(code: &[char], anomeric: usize) -> Option<usize> {
    code.iter()
        .enumerate()
        .skip(anomeric + 1)
        .filter(|(_, character)| matches!(character, '1' | '2' | 'x'))
        .take(4)
        .map(|(index, _)| index)
        .last()
        .filter(|index| matches!(code[*index], '1' | '2'))
}

fn anomer_from_codes(anomeric: char, reference: char) -> AnomericSymbol {
    if anomeric == reference {
        AnomericSymbol::Alpha
    } else {
        AnomericSymbol::Beta
    }
}

/// Build a residue from its parts through the WURCS reader, so that it has
/// exactly the shape of a residue read from WURCS text.
pub(crate) fn residue_from_parts(
    code: &[char],
    anomeric_position: u8,
    anomer: AnomericSymbol,
    ring_end: u8,
    substituents: &[(u8, String)],
) -> Option<Monosaccharide> {
    let mut modifications = substituents.to_vec();
    modifications.sort();
    modifications.dedup();
    let mut unique = format!(
        "{}-{}{}_{}-{}",
        code.iter().collect::<String>(),
        anomeric_position,
        anomer.to_char(),
        anomeric_position,
        ring_end
    );
    for (position, descriptor) in modifications {
        unique.push_str(&format!("_{position}*{descriptor}"));
    }
    let graph = parse_wurcs(&format!("WURCS=2.0/1,1,0/[{unique}]/1/")).ok()?;
    graph.root().and_then(|root| graph.residue(root)).cloned()
}

/// The WURCS UniqueRES of one residue.
pub(crate) fn unique_residue(residue: &Monosaccharide) -> String {
    let mut graph = ResidueGraph::new();
    let node = graph.add_residue(residue.clone());
    graph.set_root(node);
    crabwurcs_core::write_wurcs(&graph)
        .ok()
        .and_then(|wurcs| {
            let start = wurcs.find('[')? + 1;
            let end = wurcs[start..].find(']')? + start;
            Some(wurcs[start..end].to_string())
        })
        .unwrap_or_else(|| residue.skeleton_code.clone())
}

/// A short IUPAC name of one residue (`α-L-IdopA2S`), or its UniqueRES.
pub(crate) fn residue_label(residue: &Monosaccharide) -> String {
    let mut graph = ResidueGraph::new();
    let node = graph.add_residue(residue.clone());
    graph.set_root(node);
    crabwurcs_iupac::write_iupac_extended(&graph)
        .ok()
        .map(|name| {
            name.split("-(")
                .next()
                .unwrap_or(&name)
                .trim_end_matches('-')
                .to_string()
        })
        .filter(|name| !name.is_empty())
        .unwrap_or_else(|| unique_residue(residue))
}

struct Merger<'a> {
    coordinates: &'a CoordinateReading,
    notes: Vec<Note>,
}

impl Merger<'_> {
    fn note(&mut self, kind: PdbWarningKind, positions: Vec<u8>, aspect: impl Into<String>) {
        let note = Note {
            kind,
            positions,
            aspect: aspect.into(),
        };
        if !self.notes.contains(&note) {
            self.notes.push(note);
        }
    }
}

/// Combine name and coordinates.
pub(crate) fn combine(name: Option<&NameReading>, coordinates: &CoordinateReading) -> Combined {
    let mut merger = Merger {
        coordinates,
        notes: Vec::new(),
    };
    let coordinates_only = coordinates_only(coordinates);
    let (coordinates_first, names_first) = match name {
        None => {
            for (index, carbon) in coordinates.carbons.iter().enumerate() {
                let position = index as u8 + 1;
                match carbon {
                    CarbonReading::Stereo { volume, .. } if volume.abs() < DISTORTED_VOLUME => {
                        merger.note(
                            PdbWarningKind::DistortedGeometry,
                            vec![position],
                            format!("C{position}"),
                        );
                    }
                    CarbonReading::Unreadable => merger.note(
                        PdbWarningKind::UnreadableGeometry,
                        vec![position],
                        format!("C{position}"),
                    ),
                    _ => {}
                }
            }
            if let Some(CarbonReading::Anomeric { code: None }) = coordinates
                .carbons
                .get(coordinates.anomeric_position as usize - 1)
            {
                merger.note(
                    PdbWarningKind::MissingAtoms,
                    vec![coordinates.anomeric_position],
                    "anomeric heteroatom",
                );
            }
            (coordinates_only.clone(), coordinates_only.clone())
        }
        Some(name) => merger.merge_with(name),
    };
    for &position in &coordinates.unrecognized {
        merger.note(
            PdbWarningKind::UnrecognizedGroup,
            vec![position],
            format!("C{position}"),
        );
    }
    Combined {
        coordinates_first,
        names_first,
        coordinates_only,
        notes: merger.notes,
    }
}

/// A bare oxygen: the hydroxyl every unsubstituted position has.
fn is_hydroxyl(substituent: &CoordinateSubstituent) -> bool {
    substituent.bare && substituent.descriptor == "O"
}

fn coordinate_code(coordinates: &CoordinateReading) -> Vec<char> {
    coordinates
        .carbons
        .iter()
        .map(|carbon| match carbon {
            CarbonReading::Stereo { volume, .. } if volume.abs() < DISTORTED_VOLUME => 'x',
            other => other.code(),
        })
        .collect()
}

fn coordinate_anomer(coordinates: &CoordinateReading, code: &[char]) -> AnomericSymbol {
    let anomeric = coordinates.anomeric_position as usize - 1;
    let Some(CarbonReading::Anomeric {
        code: Some((side, volume)),
    }) = coordinates.carbons.get(anomeric)
    else {
        return AnomericSymbol::Unknown;
    };
    if volume.abs() < DISTORTED_VOLUME {
        return AnomericSymbol::Unknown;
    }
    match reference_index(code, anomeric) {
        Some(reference) => anomer_from_codes(*side, code[reference]),
        None => AnomericSymbol::Unknown,
    }
}

fn coordinates_only(coordinates: &CoordinateReading) -> Monosaccharide {
    let code = coordinate_code(coordinates);
    let anomer = coordinate_anomer(coordinates, &code);
    let substituents = coordinates
        .substituents
        .iter()
        .filter(|substituent| !is_hydroxyl(substituent))
        .map(|substituent| (substituent.position, substituent.descriptor.clone()))
        .collect::<Vec<_>>();
    residue_from_parts(
        &code,
        coordinates.anomeric_position,
        anomer,
        coordinates.ring_end,
        &substituents,
    )
    .unwrap_or_else(|| {
        crabwurcs_core::residue_from_kind(crabwurcs_core::ResidueKind::Hex)
            .expect("the generic hexose registry entry is valid")
    })
}

impl Merger<'_> {
    fn merge_with(&mut self, name: &NameReading) -> (Monosaccharide, Monosaccharide) {
        let coordinates = self.coordinates;
        let complete = coordinates.complete;
        let name_code = full_code(&name.residue);
        let coordinate_code = coordinate_code(coordinates);
        let label_contradiction = |merger: &mut Self, positions: Vec<u8>, aspect: &str| {
            merger.note(
                PdbWarningKind::NameContradictsCoordinates,
                positions,
                aspect,
            );
        };

        // A different backbone: another monosaccharide altogether.
        let truncated = name_code.len() > coordinate_code.len();
        if name.residue.anomeric_position != coordinates.anomeric_position
            || name_code.len() < coordinate_code.len()
            || (truncated && complete)
        {
            label_contradiction(self, vec![], "backbone");
            let external = coordinates
                .substituents
                .iter()
                .filter(|substituent| {
                    substituent.external && (substituent.position as usize) <= name_code.len()
                })
                .map(|substituent| (substituent.position, substituent.descriptor.clone()))
                .collect::<Vec<_>>();
            let mut names_first = name.residue.clone();
            for (position, descriptor) in external {
                names_first
                    .modifications
                    .push(crabwurcs_core::Modification {
                        position: crabwurcs_core::CarbonPosition(position),
                        descriptor,
                        probability: None,
                    });
            }
            return (coordinates_only(coordinates), names_first);
        }

        // Carbon by carbon.
        let usable = if truncated {
            coordinate_code.len() - 1
        } else {
            coordinate_code.len()
        };
        let mut coordinates_first = Vec::with_capacity(name_code.len());
        let mut names_first = Vec::with_capacity(name_code.len());
        for (index, &named) in name_code.iter().enumerate() {
            let position = index as u8 + 1;
            let aspect = format!("C{position}");
            if index >= usable {
                self.note(PdbWarningKind::MissingAtoms, vec![position], aspect);
                coordinates_first.push(named);
                names_first.push(named);
                continue;
            }
            let reading = coordinates.carbons[index];
            let (choice, kind) = match reading {
                CarbonReading::Anomeric { .. } => (
                    'a',
                    (named != 'a').then_some(PdbWarningKind::NameContradictsCoordinates),
                ),
                CarbonReading::Stereo { code, volume } => {
                    // A centre too flat to read is the name's, whether or not
                    // its sign happens to agree.
                    if matches!(named, '1' | '2') && volume.abs() < DISTORTED_VOLUME {
                        (named, Some(PdbWarningKind::DistortedGeometry))
                    } else if code == named {
                        (code, None)
                    } else {
                        (code, Some(PdbWarningKind::NameContradictsCoordinates))
                    }
                }
                CarbonReading::NoHeteroatom { .. } => match named {
                    'd' => ('d', None),
                    _ if complete => ('d', Some(PdbWarningKind::NameContradictsCoordinates)),
                    _ => (named, Some(PdbWarningKind::MissingAtoms)),
                },
                CarbonReading::Terminal { code, .. } => {
                    match (oxygen_count(code), oxygen_count(named)) {
                        _ if code == named => (code, None),
                        (Some(read), Some(expected)) if read < expected && !complete => {
                            (named, Some(PdbWarningKind::MissingAtoms))
                        }
                        _ => (code, Some(PdbWarningKind::NameContradictsCoordinates)),
                    }
                }
                CarbonReading::Unreadable => (named, Some(PdbWarningKind::UnreadableGeometry)),
            };
            if let Some(kind) = kind {
                self.note(kind, vec![position], aspect);
            }
            coordinates_first.push(choice);
            names_first.push(
                if kind == Some(PdbWarningKind::NameContradictsCoordinates) {
                    named
                } else {
                    choice
                },
            );
        }

        // Ring size.
        let name_ring_end = name.residue.ring_end.unwrap_or(coordinates.ring_end);
        let (coordinate_ring, name_ring) = if name_ring_end != coordinates.ring_end {
            label_contradiction(self, vec![], "ring");
            (coordinates.ring_end, name_ring_end)
        } else {
            (coordinates.ring_end, coordinates.ring_end)
        };

        // Anomer.
        let anomeric = coordinates.anomeric_position as usize - 1;
        let read = coordinate_anomer(coordinates, &coordinate_code);
        let named = name.residue.anomeric_symbol;
        let (coordinate_anomer, name_anomer) = match (read, named) {
            (AnomericSymbol::Unknown, _) => {
                if matches!(named, AnomericSymbol::Alpha | AnomericSymbol::Beta) {
                    let kind = match coordinates.carbons.get(anomeric) {
                        Some(CarbonReading::Anomeric { code: Some(_) }) => {
                            PdbWarningKind::DistortedGeometry
                        }
                        _ => PdbWarningKind::MissingAtoms,
                    };
                    self.note(kind, vec![coordinates.anomeric_position], "anomer");
                }
                // The reference centre may have come from the name.
                let anomer = match coordinates.carbons.get(anomeric) {
                    Some(CarbonReading::Anomeric {
                        code: Some((side, volume)),
                    }) if volume.abs() >= DISTORTED_VOLUME => {
                        reference_index(&coordinates_first, anomeric)
                            .map(|reference| anomer_from_codes(*side, coordinates_first[reference]))
                            .unwrap_or(named)
                    }
                    _ => named,
                };
                (anomer, named)
            }
            (read, AnomericSymbol::Alpha | AnomericSymbol::Beta) if read != named => {
                label_contradiction(self, vec![coordinates.anomeric_position], "anomer");
                (read, named)
            }
            (read, _) => (read, read),
        };

        // Substituents.
        let (coordinate_substituents, name_substituents) = self.merge_substituents(name);

        let coordinates_first = residue_from_parts(
            &coordinates_first,
            coordinates.anomeric_position,
            coordinate_anomer,
            coordinate_ring,
            &coordinate_substituents,
        )
        .unwrap_or_else(|| coordinates_only(coordinates));
        let names_first = residue_from_parts(
            &names_first,
            name.residue.anomeric_position,
            name_anomer,
            name_ring,
            &name_substituents,
        )
        .unwrap_or_else(|| name.residue.clone());
        (
            keep_name_kind(coordinates_first, name),
            keep_name_kind(names_first, name),
        )
    }

    /// Substituents position by position: (coordinates first, names first).
    fn merge_substituents(&mut self, name: &NameReading) -> (Substituents, Substituents) {
        let coordinates = self.coordinates;
        let complete = coordinates.complete;
        let mut coordinates_first = Vec::new();
        let mut names_first = Vec::new();
        let mut positions = name
            .residue
            .modifications
            .iter()
            .map(|modification| modification.position.0)
            .chain(
                coordinates
                    .substituents
                    .iter()
                    .map(|substituent| substituent.position),
            )
            .collect::<Vec<_>>();
        positions.sort_unstable();
        positions.dedup();
        for position in positions {
            let aspect = format!("C{position} substituent");
            let named = name
                .residue
                .modifications
                .iter()
                .filter(|modification| modification.position.0 == position)
                .map(|modification| modification.descriptor.clone())
                .collect::<Vec<_>>();
            let read = coordinates
                .substituents
                .iter()
                .filter(|substituent| substituent.position == position)
                .collect::<Vec<_>>();
            let unrecognized = coordinates.unrecognized.contains(&position);
            let mut matched = vec![false; read.len()];
            for descriptor in &named {
                let same = read
                    .iter()
                    .position(|substituent| &substituent.descriptor == descriptor)
                    .or_else(|| {
                        // An amino group the name lists, extended by a group
                        // in another residue (GLYCAM's SO3 on N2) or by one
                        // the name does not spell out.
                        (descriptor == "N")
                            .then(|| {
                                read.iter().position(|substituent| {
                                    (substituent.external || name.amines_may_carry_groups())
                                        && substituent.descriptor.starts_with('N')
                                })
                            })
                            .flatten()
                    });
                if let Some(index) = same {
                    matched[index] = true;
                    coordinates_first.push((position, read[index].descriptor.clone()));
                    names_first.push((position, read[index].descriptor.clone()));
                    continue;
                }
                // The coordinates show the bare heteroatom (or a group they
                // cannot name): atoms missing unless the model is complete.
                let bare = read.iter().position(|substituent| {
                    substituent.bare && descriptor.starts_with(&substituent.descriptor[..1])
                });
                let nothing = read.iter().all(|substituent| substituent.bare);
                if unrecognized || (nothing && !complete) {
                    if !unrecognized {
                        self.note(PdbWarningKind::MissingAtoms, vec![position], aspect.clone());
                    }
                    if let Some(index) = bare {
                        matched[index] = true;
                    }
                    coordinates_first.push((position, descriptor.clone()));
                    names_first.push((position, descriptor.clone()));
                    continue;
                }
                self.note(
                    PdbWarningKind::NameContradictsCoordinates,
                    vec![position],
                    aspect.clone(),
                );
                names_first.push((position, descriptor.clone()));
                if let Some(index) = bare {
                    matched[index] = true;
                    if !is_hydroxyl(read[index]) {
                        coordinates_first.push((position, read[index].descriptor.clone()));
                    }
                }
            }
            for (index, substituent) in read.iter().enumerate() {
                if matched[index] || is_hydroxyl(substituent) {
                    continue;
                }
                coordinates_first.push((position, substituent.descriptor.clone()));
                if substituent.external || !name.lists_own_substituents() {
                    if named.is_empty() || substituent.external {
                        names_first.push((position, substituent.descriptor.clone()));
                    }
                } else {
                    self.note(
                        PdbWarningKind::NameContradictsCoordinates,
                        vec![position],
                        aspect.clone(),
                    );
                }
            }
        }
        (coordinates_first, names_first)
    }
}

/// Keep the registry identity the name declares when the result is the
/// residue the name describes (with or without extra substituents).
fn keep_name_kind(mut residue: Monosaccharide, name: &NameReading) -> Monosaccharide {
    let same = residue.skeleton_code == name.residue.skeleton_code
        && residue.anomeric_prefix == name.residue.anomeric_prefix
        && residue.anomeric_position == name.residue.anomeric_position
        && residue.ring_end == name.residue.ring_end
        && name
            .residue
            .modifications
            .iter()
            .all(|modification| residue.modifications.contains(modification));
    if same {
        residue.residue_kind = name.residue.residue_kind;
    }
    residue
}
