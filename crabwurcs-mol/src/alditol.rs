//! Reading alditols from a molecular graph.
//!
//! An alditol is an open-chain polyol with no anomeric carbon, most often the
//! reducing end of a reduced glycan (GalNAc-ol). WURCS writes one without an
//! anomeric prefix, as `h2112h`. It has no ring for [`crate::detect_sugar_rings`]
//! to find, so it is found as the unbranched carbon chain that carries a
//! hydroxyl or an amide on every carbon. The caller then has to check that the
//! glycan read accounts for the whole molecule: a fatty ester on a glycerol
//! looks like an acetate to the substituent readers.

use chematic::core::{AtomIdx, Molecule};
use crabwurcs_core::{Monosaccharide, parse_wurcs};

use crate::{
    DetectedSugarRing, atom_is, backbone_atoms, chain_stereo_descriptor, n_acyl_descriptor,
    n_sulfate_descriptor, oxygen_substituent_descriptor,
};

/// An alditol found in a molecule.
#[derive(Debug, Clone)]
pub(crate) struct DetectedAlditol {
    /// Backbone carbons, C1 first.
    backbone: Vec<AtomIdx>,
    pub(crate) residue: Monosaccharide,
}

impl DetectedAlditol {
    /// The WURCS position of a backbone carbon.
    pub(crate) fn position_of(&self, carbon: AtomIdx) -> Option<u8> {
        self.backbone
            .iter()
            .position(|atom| *atom == carbon)
            .map(|index| index as u8 + 1)
    }
}

/// One way of numbering a chain.
struct Numbering {
    backbone: Vec<AtomIdx>,
    skeleton: String,
    substituents: Vec<(u8, &'static str)>,
    /// Positions that carry a sugar.
    attachments: Vec<u8>,
}

impl Numbering {
    /// Lower is preferred. A polyol has no anomeric carbon to number from, and
    /// the two directions are one molecule (D-glucitol is L-gulitol), so this
    /// picks the numbering of the sugar it came from: the deoxy end last, an
    /// amide on the lowest carbon, a D highest stereocentre. What is left is
    /// a real symmetry, and the lowest positions decide.
    fn preference(&self) -> impl Ord + use<> {
        let amide = self
            .substituents
            .iter()
            .filter(|(_, descriptor)| descriptor.starts_with('N'))
            .map(|(position, _)| *position)
            .min()
            .unwrap_or(u8::MAX);
        let highest_stereocentre = self
            .skeleton
            .chars()
            .rev()
            .find(|descriptor| matches!(descriptor, '1' | '2'));
        (
            self.skeleton.starts_with('m'),
            amide,
            highest_stereocentre == Some('1'),
            self.attachments.clone(),
            self.skeleton.clone(),
        )
    }
}

/// The unbranched carbon chain through `start`, end to end, leaving out
/// `sugars`. `None` for a branch or a ring.
fn carbon_chain(molecule: &Molecule, start: AtomIdx, sugars: &[AtomIdx]) -> Option<Vec<AtomIdx>> {
    let carbons_beside = |atom: AtomIdx| {
        molecule
            .neighbors(atom)
            .map(|(neighbor, _)| neighbor)
            .filter(|neighbor| atom_is(molecule, *neighbor, "C") && !sugars.contains(neighbor))
            .collect::<Vec<_>>()
    };
    let mut seen = vec![start];
    let mut stack = vec![start];
    while let Some(atom) = stack.pop() {
        let beside = carbons_beside(atom);
        if beside.len() > 2 {
            return None;
        }
        for neighbor in beside {
            if !seen.contains(&neighbor) {
                seen.push(neighbor);
                stack.push(neighbor);
            }
        }
    }
    let mut end = seen
        .iter()
        .copied()
        .find(|atom| carbons_beside(*atom).len() == 1)?;
    let mut chain = vec![end];
    while let Some(next) = carbons_beside(end)
        .into_iter()
        .find(|neighbor| !chain.contains(neighbor))
    {
        chain.push(next);
        end = next;
    }
    (chain.len() == seen.len()).then_some(chain)
}

/// The heteroatom on a chain carbon, if it is the one single-bonded heteroatom
/// a polyol has there. The outer `None` rejects the carbon.
fn heteroatom(molecule: &Molecule, carbon: AtomIdx, chain: &[AtomIdx]) -> Option<Option<AtomIdx>> {
    let mut found = None;
    for (neighbor, bond) in molecule.neighbors(carbon) {
        if chain.contains(&neighbor) || atom_is(molecule, neighbor, "H") {
            continue;
        }
        let single = molecule.bond(bond).order.order_int() == 1;
        let polyol_atom = atom_is(molecule, neighbor, "O") || atom_is(molecule, neighbor, "N");
        if !single || !polyol_atom || found.is_some() {
            return None;
        }
        found = Some(neighbor);
    }
    Some(found)
}

fn numbering(molecule: &Molecule, backbone: Vec<AtomIdx>, sugars: &[AtomIdx]) -> Option<Numbering> {
    let known = sugars
        .iter()
        .chain(&backbone)
        .copied()
        .collect::<Vec<AtomIdx>>();
    let last = backbone.len() - 1;
    let mut skeleton = String::with_capacity(backbone.len());
    let mut substituents = Vec::new();
    let mut attachments = Vec::new();
    for (index, carbon) in backbone.iter().copied().enumerate() {
        let position = index as u8 + 1;
        let atom = heteroatom(molecule, carbon, &backbone)?;
        if index == 0 || index == last {
            skeleton.push(if atom.is_some() { 'h' } else { 'm' });
        } else {
            atom?;
            skeleton.push(chain_stereo_descriptor(molecule, &backbone, index));
        }
        let Some(atom) = atom else {
            continue;
        };
        let bonded_to_sugar = molecule
            .neighbors(atom)
            .any(|(neighbor, _)| neighbor != carbon && sugars.contains(&neighbor));
        if bonded_to_sugar {
            attachments.push(position);
        }
        let descriptor = n_acyl_descriptor(molecule, carbon)
            .or_else(|| n_sulfate_descriptor(molecule, carbon))
            .or_else(|| oxygen_substituent_descriptor(molecule, carbon, &known));
        if let Some(descriptor) = descriptor {
            substituents.push((position, descriptor));
        }
    }
    Some(Numbering {
        backbone,
        skeleton,
        substituents,
        attachments,
    })
}

/// The alditol whose chain runs through `start`, if there is one.
fn alditol_at(molecule: &Molecule, start: AtomIdx, sugars: &[AtomIdx]) -> Option<DetectedAlditol> {
    let chain = carbon_chain(molecule, start, sugars)?;
    if chain.len() < 3 {
        return None;
    }
    let reversed = chain.iter().rev().copied().collect();
    let forward = numbering(molecule, chain, sugars)?;
    let backward = numbering(molecule, reversed, sugars)?;
    // Three hydroxylated carbons make a polyol; a glycol is not a sugar.
    let hydroxylated = forward.skeleton.chars().filter(|d| *d != 'm').count();
    if hydroxylated < 3 || !forward.skeleton.contains('h') {
        return None;
    }
    let chosen = if backward.preference() < forward.preference() {
        backward
    } else {
        forward
    };

    let mut unique = chosen.skeleton;
    for (position, descriptor) in &chosen.substituents {
        unique.push_str(&format!("_{position}*{descriptor}"));
    }
    let graph = parse_wurcs(&format!("WURCS=2.0/1,1,0/[{unique}]/1/")).ok()?;
    let residue = graph.residue(graph.root()?)?.clone();
    Some(DetectedAlditol {
        backbone: chosen.backbone,
        residue,
    })
}

/// A free alditol: a molecule without a sugar ring.
pub(crate) fn free_alditol(molecule: &Molecule) -> Option<DetectedAlditol> {
    molecule
        .atoms()
        .filter(|(_, atom)| atom.element.symbol() == "C")
        .find_map(|(index, _)| alditol_at(molecule, index, &[]))
}

/// The alditol a sugar is glycosidically bonded to.
pub(crate) fn aglycone_alditol(
    molecule: &Molecule,
    rings: &[DetectedSugarRing],
) -> Option<DetectedAlditol> {
    let sugars = rings.iter().flat_map(backbone_atoms).collect::<Vec<_>>();
    rings.iter().find_map(|ring| {
        let aglycone = molecule
            .neighbors(ring.anomeric_oxygen)
            .map(|(neighbor, _)| neighbor)
            .find(|neighbor| *neighbor != ring.carbons[0])?;
        if sugars.contains(&aglycone) || !atom_is(molecule, aglycone, "C") {
            return None;
        }
        alditol_at(molecule, aglycone, &sugars)
    })
}
