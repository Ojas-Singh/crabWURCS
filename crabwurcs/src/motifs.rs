//! Named glycan motifs, usable wherever a motif string is accepted.
//!
//! Each definition is an IUPAC pattern with the position its residues must
//! take. An epitope that is a cap by definition requires its non-reducing
//! residues to be terminal, so `LewisX` is not found inside sialyl-Lewis X; a
//! core requires its last residue to be the reducing end, so `CoreFucose` is
//! not an antennal fucose.

use std::sync::OnceLock;

use crabwurcs_core::{MotifPosition, MotifQuery};

use crate::{ConversionResult, Format, iupac, parse_notation};

const F: MotifPosition = MotifPosition::Flexible;
const T: MotifPosition = MotifPosition::Terminal;
const R: MotifPosition = MotifPosition::Reducing;

/// A motif of the bundled library.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NamedMotif {
    pub name: &'static str,
    /// Other names it is looked up by. Lookups ignore case, spaces and
    /// hyphens, so `sialyl-Lewis X` finds `SialylLewisX`.
    pub aliases: &'static [&'static str],
    pub iupac: &'static str,
    /// One position per residue of `iupac`, in reading order.
    pub positions: &'static [MotifPosition],
}

const fn motif(
    name: &'static str,
    aliases: &'static [&'static str],
    iupac: &'static str,
    positions: &'static [MotifPosition],
) -> NamedMotif {
    NamedMotif {
        name,
        aliases,
        iupac,
        positions,
    }
}

const LIBRARY: &[NamedMotif] = &[
    // Lewis and blood-group epitopes.
    motif(
        "LewisX",
        &["LeX", "CD15"],
        "Gal(b1-4)[Fuc(a1-3)]GlcNAc",
        &[T, T, F],
    ),
    motif("LewisA", &["LeA"], "Gal(b1-3)[Fuc(a1-4)]GlcNAc", &[T, T, F]),
    motif(
        "LewisY",
        &["LeY"],
        "Fuc(a1-2)Gal(b1-4)[Fuc(a1-3)]GlcNAc",
        &[T, F, T, F],
    ),
    motif(
        "LewisB",
        &["LeB"],
        "Fuc(a1-2)Gal(b1-3)[Fuc(a1-4)]GlcNAc",
        &[T, F, T, F],
    ),
    motif(
        "SialylLewisX",
        &["sLeX", "CD15s"],
        "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc",
        &[T, F, T, F],
    ),
    motif(
        "SialylLewisA",
        &["sLeA", "CA19-9"],
        "Neu5Ac(a2-3)Gal(b1-3)[Fuc(a1-4)]GlcNAc",
        &[T, F, T, F],
    ),
    motif("HType1", &["H1"], "Fuc(a1-2)Gal(b1-3)GlcNAc", &[T, F, F]),
    motif("HType2", &["H2"], "Fuc(a1-2)Gal(b1-4)GlcNAc", &[T, F, F]),
    motif(
        "BloodGroupA",
        &["AAntigen"],
        "GalNAc(a1-3)[Fuc(a1-2)]Gal",
        &[T, T, F],
    ),
    motif(
        "BloodGroupB",
        &["BAntigen"],
        "Gal(a1-3)[Fuc(a1-2)]Gal",
        &[T, T, F],
    ),
    motif(
        "AlphaGal",
        &["aGal"],
        "Gal(a1-3)Gal(b1-4)GlcNAc",
        &[T, F, F],
    ),
    motif(
        "Sda",
        &["SdaAntigen"],
        "GalNAc(b1-4)[Neu5Ac(a2-3)]Gal",
        &[T, T, F],
    ),
    // Antennal disaccharides and their caps.
    motif("LacNAc", &["LacNAcType2"], "Gal(b1-4)GlcNAc", &[F, F]),
    motif("LacNAcType1", &["Type1LacNAc"], "Gal(b1-3)GlcNAc", &[F, F]),
    motif("LacdiNAc", &["LDN"], "GalNAc(b1-4)GlcNAc", &[F, F]),
    motif(
        "PolyLacNAc",
        &[],
        "Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)GlcNAc",
        &[F, F, F, F],
    ),
    motif(
        "Sialyl3LacNAc",
        &["a23SialylLacNAc"],
        "Neu5Ac(a2-3)Gal(b1-4)GlcNAc",
        &[T, F, F],
    ),
    motif(
        "Sialyl6LacNAc",
        &["a26SialylLacNAc"],
        "Neu5Ac(a2-6)Gal(b1-4)GlcNAc",
        &[T, F, F],
    ),
    // N-glycan core features.
    motif(
        "NGlycanCore",
        &["ChitobioseCore", "TrimannosylCore"],
        "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc",
        &[F, F, F, F, R],
    ),
    motif("CoreFucose", &["CoreFuc"], "Fuc(a1-6)GlcNAc", &[T, R]),
    motif(
        "CoreFucoseAlpha3",
        &["CoreFucA3"],
        "Fuc(a1-3)GlcNAc",
        &[T, R],
    ),
    motif(
        "BisectingGlcNAc",
        &["Bisecting"],
        "Man(a1-3)[Man(a1-6)][GlcNAc(b1-4)]Man(b1-4)GlcNAc",
        &[F, F, F, F, F],
    ),
    // Mucin-type O-glycan cores on their reducing-end GalNAc.
    motif("Core1", &["TAntigen"], "Gal(b1-3)GalNAc", &[F, R]),
    motif("Core2", &[], "Gal(b1-3)[GlcNAc(b1-6)]GalNAc", &[F, F, R]),
    motif("Core3", &[], "GlcNAc(b1-3)GalNAc", &[F, R]),
    motif("Core4", &[], "GlcNAc(b1-3)[GlcNAc(b1-6)]GalNAc", &[F, F, R]),
    motif("SialylTn", &["sTn"], "Neu5Ac(a2-6)GalNAc", &[T, R]),
];

fn normalized(name: &str) -> String {
    name.chars()
        .filter(|character| !matches!(character, ' ' | '-' | '_'))
        .flat_map(char::to_lowercase)
        .collect()
}

/// Every motif of the library.
pub fn named_motifs() -> &'static [NamedMotif] {
    LIBRARY
}

/// The library motif called `name`.
pub fn named_motif(name: &str) -> Option<&'static NamedMotif> {
    let wanted = normalized(name);
    LIBRARY.iter().find(|motif| {
        normalized(motif.name) == wanted
            || motif
                .aliases
                .iter()
                .any(|alias| normalized(alias) == wanted)
    })
}

impl NamedMotif {
    /// This motif as a query.
    pub fn query(&self) -> MotifQuery {
        queries()[LIBRARY
            .iter()
            .position(|motif| motif.name == self.name)
            .expect("a named motif comes from the library")]
        .clone()
    }
}

fn queries() -> &'static [MotifQuery] {
    static QUERIES: OnceLock<Vec<MotifQuery>> = OnceLock::new();
    QUERIES.get_or_init(|| {
        LIBRARY
            .iter()
            .map(|motif| {
                let graph = iupac::parse_iupac_condensed(motif.iupac)
                    .expect("every library motif is valid IUPAC");
                MotifQuery::with_positions(graph, motif.positions)
                    .expect("every library motif states all its positions")
            })
            .collect()
    })
}

/// The motif `text` stands for. With [`Format::Auto`] a name of the library
/// wins; anything else is parsed as a notation of that format.
pub fn resolve_motif(text: &str, format: Format) -> ConversionResult<MotifQuery> {
    if format == Format::Auto
        && let Some(motif) = named_motif(text)
    {
        return Ok(motif.query());
    }
    Ok(MotifQuery::new(parse_notation(text, format)?))
}

/// The names of the library motifs found in `graph`.
pub fn annotate_motifs(graph: &crabwurcs_core::ResidueGraph) -> Vec<&'static str> {
    LIBRARY
        .iter()
        .zip(queries())
        .filter(|(_, query)| query.contained_in(graph).unwrap_or(false))
        .map(|(motif, _)| motif.name)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ResidueKind;

    fn glycan(iupac: &str) -> crate::ResidueGraph {
        parse_notation(iupac, Format::IupacCondensed).unwrap()
    }

    fn found_in(iupac: &str) -> Vec<&'static str> {
        annotate_motifs(&glycan(iupac))
    }

    #[test]
    fn every_library_motif_parses_with_a_position_per_residue() {
        for motif in named_motifs() {
            let query = motif.query();
            assert_eq!(
                query.motif().node_count(),
                motif.positions.len(),
                "{}",
                motif.name
            );
            // A motif finds itself unless its positions need more around it.
            if motif.positions.iter().all(|position| *position == F) {
                assert!(
                    query.contained_in(&glycan(motif.iupac)).unwrap(),
                    "{}",
                    motif.name
                );
            }
        }
    }

    #[test]
    fn no_motif_name_shadows_a_residue() {
        let residues = ResidueKind::ALL
            .iter()
            .flat_map(|kind| {
                std::iter::once(kind.canonical_name()).chain(kind.aliases().iter().copied())
            })
            .map(normalized)
            .collect::<std::collections::HashSet<_>>();
        let mut names = std::collections::HashSet::new();
        for motif in named_motifs() {
            for name in std::iter::once(motif.name).chain(motif.aliases.iter().copied()) {
                let name = normalized(name);
                assert!(!residues.contains(&name), "{name} is a residue");
                assert!(names.insert(name.clone()), "{name} is used twice");
            }
        }
    }

    #[test]
    fn lookups_ignore_case_spaces_and_hyphens() {
        for name in ["LewisX", "lewis x", "Le-X", "LEX", "CD15"] {
            assert_eq!(named_motif(name).unwrap().name, "LewisX", "{name}");
        }
        assert_eq!(named_motif("sialyl-Lewis X").unwrap().name, "SialylLewisX");
        assert_eq!(named_motif("core fucose").unwrap().name, "CoreFucose");
        assert_eq!(named_motif("Core 1").unwrap().name, "Core1");
        assert_eq!(
            named_motif("bisecting GlcNAc").unwrap().name,
            "BisectingGlcNAc"
        );
        assert!(named_motif("Gal(b1-4)GlcNAc").is_none());
    }

    #[test]
    fn a_name_wins_only_when_the_format_is_left_open() {
        let target = glycan("Gal(b1-4)[Fuc(a1-3)]GlcNAc");
        let named = resolve_motif("LewisX", Format::Auto).unwrap();
        assert_eq!(named.count_in(&target).unwrap(), 1);
        let notation = resolve_motif("Gal(b1-4)GlcNAc", Format::Auto).unwrap();
        assert_eq!(notation.count_in(&target).unwrap(), 1);
        // Asked for as IUPAC, the same word is an unknown residue name.
        let literal = resolve_motif("LewisX", Format::IupacCondensed).unwrap();
        assert_eq!(literal.count_in(&target).unwrap(), 0);
    }

    #[test]
    fn caps_are_terminal_so_an_extended_epitope_is_not_the_shorter_one() {
        let lewis_x = found_in("Gal(b1-4)[Fuc(a1-3)]GlcNAc");
        assert!(lewis_x.contains(&"LewisX"));
        assert!(!lewis_x.contains(&"SialylLewisX"));

        let sialyl = found_in("Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc");
        assert!(sialyl.contains(&"SialylLewisX"));
        assert!(sialyl.contains(&"Sialyl3LacNAc"), "{sialyl:?}");
        assert!(!sialyl.contains(&"LewisX"), "{sialyl:?}");
    }

    #[test]
    fn core_fucose_is_on_the_reducing_end_only() {
        let core = found_in("Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc");
        assert!(core.contains(&"CoreFucose"));
        let antennal = found_in("Fuc(a1-6)GlcNAc(b1-2)Man(a1-3)Man");
        assert!(!antennal.contains(&"CoreFucose"), "{antennal:?}");
    }

    #[test]
    fn an_n_glycan_is_annotated_with_what_it_carries() {
        let complex = found_in(
            "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-4)][Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc",
        );
        for expected in [
            "NGlycanCore",
            "CoreFucose",
            "BisectingGlcNAc",
            "LacNAc",
            "Sialyl3LacNAc",
        ] {
            assert!(complex.contains(&expected), "{expected}: {complex:?}");
        }

        let unbisected =
            found_in("Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc");
        assert!(unbisected.contains(&"NGlycanCore"));
        assert!(!unbisected.contains(&"BisectingGlcNAc"));
        assert!(!unbisected.contains(&"CoreFucose"));
    }

    #[test]
    fn o_glycan_cores_sit_on_the_reducing_galnac() {
        assert!(found_in("Gal(b1-3)GalNAc").contains(&"Core1"));
        let core2 = found_in("Gal(b1-3)[GlcNAc(b1-6)]GalNAc");
        assert!(core2.contains(&"Core2"));
        assert!(found_in("Neu5Ac(a2-6)GalNAc").contains(&"SialylTn"));
        // The same disaccharide away from the reducing end is no core.
        assert!(!found_in("Gal(b1-3)GalNAc(b1-4)Gal").contains(&"Core1"));
    }
}
