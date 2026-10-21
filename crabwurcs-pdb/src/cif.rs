//! The few mmCIF items the extractor reads itself, and the PDB `LINK`
//! records: covalent connections between residues, and the crystal items
//! some archive files state in a way pdbtbx refuses.
//!
//! pdbtbx reads the coordinates (`_atom_site`) but neither `_struct_conn`
//! nor `LINK`. It refuses a file whose space-group items disagree
//! (`_symmetry.space_group_name_H-M 'I 1 2 1'` with
//! `_symmetry.space_group_name_Hall 'C 2y (x,y,-x+z)'`) and panics on an
//! unknown (`?`) space group or cell; the extractor needs neither.

use std::collections::HashMap;

/// One CIF token: a tag, a reserved word or a value, with the byte span it
/// takes in the text (quotes and text-field delimiters included).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Token<'a> {
    text: &'a str,
    /// A quoted string or a text field: always a value, never a tag or a
    /// reserved word.
    quoted: bool,
    start: usize,
    end: usize,
}

impl Token<'_> {
    fn is_tag(&self) -> bool {
        !self.quoted && self.text.starts_with('_')
    }

    fn is_reserved(&self) -> bool {
        if self.quoted {
            return false;
        }
        let lower = self.text.to_ascii_lowercase();
        lower == "loop_"
            || lower == "stop_"
            || lower == "global_"
            || lower.starts_with("data_")
            || lower.starts_with("save_")
    }

    /// The value, `None` for the CIF null values `?` and `.`.
    fn value(&self) -> Option<&str> {
        (self.quoted || !matches!(self.text, "?" | ".")).then_some(self.text)
    }
}

/// CIF 1.1 tokens from byte `start` on: `#` comments are skipped; a quoted
/// string ends at its quote followed by whitespace; a text field runs from a
/// `;` at the start of a line to the next line starting with `;`.
struct Tokens<'a> {
    text: &'a str,
    cursor: usize,
}

impl<'a> Tokens<'a> {
    fn new(text: &'a str, start: usize) -> Self {
        Self {
            text,
            cursor: start,
        }
    }
}

impl<'a> Iterator for Tokens<'a> {
    type Item = Token<'a>;

    fn next(&mut self) -> Option<Token<'a>> {
        let bytes = self.text.as_bytes();
        loop {
            while self.cursor < bytes.len() && bytes[self.cursor].is_ascii_whitespace() {
                self.cursor += 1;
            }
            if self.cursor >= bytes.len() {
                return None;
            }
            if bytes[self.cursor] != b'#' {
                break;
            }
            while self.cursor < bytes.len() && bytes[self.cursor] != b'\n' {
                self.cursor += 1;
            }
        }
        let start = self.cursor;
        let line_start = start == 0 || bytes[start - 1] == b'\n' || bytes[start - 1] == b'\r';
        if line_start && bytes[start] == b';' {
            let mut end = start + 1;
            let close = loop {
                match self.text[end..].find('\n') {
                    Some(offset) => {
                        end += offset + 1;
                        if bytes.get(end) == Some(&b';') {
                            break Some(end);
                        }
                    }
                    None => break None,
                }
            };
            let (value_end, token_end) = match close {
                Some(close) => (close - 1, close + 1),
                None => (bytes.len(), bytes.len()),
            };
            self.cursor = token_end;
            return Some(Token {
                text: self.text[start + 1..value_end].trim_end_matches('\r'),
                quoted: true,
                start,
                end: token_end,
            });
        }
        if matches!(bytes[start], b'\'' | b'"') {
            let quote = bytes[start];
            let mut end = start + 1;
            while end < bytes.len()
                && !(bytes[end] == quote
                    && bytes
                        .get(end + 1)
                        .is_none_or(|next| next.is_ascii_whitespace()))
            {
                end += 1;
            }
            self.cursor = (end + 1).min(bytes.len());
            return Some(Token {
                text: &self.text[start + 1..end.min(bytes.len())],
                quoted: true,
                start,
                end: self.cursor,
            });
        }
        let mut end = start;
        while end < bytes.len() && !bytes[end].is_ascii_whitespace() {
            end += 1;
        }
        self.cursor = end;
        Some(Token {
            text: &self.text[start..end],
            quoted: false,
            start,
            end,
        })
    }
}

/// Byte offset of the first line that starts with `prefix`, after leading
/// blanks.
fn line_starting_with(text: &str, prefix: &str) -> Option<usize> {
    let mut offset = 0;
    for line in text.split_inclusive('\n') {
        let indent = line.len() - line.trim_start_matches([' ', '\t']).len();
        if line[indent..].starts_with(prefix) {
            return Some(offset + indent);
        }
        offset += line.len();
    }
    None
}

/// The rows of one category (`_struct_conn`), looped or as single items, as
/// maps from the item name (`ptnr1_label_atom_id`) to its value; `?` and `.`
/// are left out.
fn category_rows(text: &str, category: &str) -> Vec<HashMap<String, String>> {
    let prefix = format!("{category}.");
    let Some(first) = line_starting_with(text, &prefix) else {
        return Vec::new();
    };
    // A loop starts with `loop_` on an earlier line (comments and blank
    // lines may come between).
    let mut start = first;
    let mut line_end = first;
    while line_end > 0 {
        let line_start = text[..line_end - 1]
            .rfind('\n')
            .map_or(0, |index| index + 1);
        let trimmed = text[line_start..line_end].trim();
        line_end = line_start;
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        if trimmed.eq_ignore_ascii_case("loop_") {
            start = line_start;
        }
        break;
    }
    let item = |tag: &str| {
        tag.get(..prefix.len())
            .filter(|head| head.eq_ignore_ascii_case(&prefix))
            .map(|_| tag[prefix.len()..].to_ascii_lowercase())
    };
    let mut tokens = Tokens::new(text, start).peekable();
    let mut rows = Vec::new();
    if tokens
        .peek()
        .is_some_and(|token| !token.quoted && token.text.eq_ignore_ascii_case("loop_"))
    {
        tokens.next();
        let mut headers = Vec::new();
        while let Some(token) = tokens.next_if(Token::is_tag) {
            headers.push(item(token.text));
        }
        if headers.is_empty() || headers.iter().all(Option::is_none) {
            return rows;
        }
        let mut row = HashMap::new();
        let mut column = 0;
        while let Some(token) = tokens.next_if(|token| !token.is_tag() && !token.is_reserved()) {
            if let (Some(Some(name)), Some(value)) = (headers.get(column), token.value()) {
                row.insert(name.clone(), value.to_string());
            }
            column += 1;
            if column == headers.len() {
                rows.push(std::mem::take(&mut row));
                column = 0;
            }
        }
        return rows;
    }
    let mut row = HashMap::new();
    while let Some(tag) = tokens.next_if(Token::is_tag) {
        let Some(name) = item(tag.text) else {
            break;
        };
        let Some(value) = tokens.next_if(|token| !token.is_tag() && !token.is_reserved()) else {
            break;
        };
        if let Some(value) = value.value() {
            row.insert(name, value.to_string());
        }
    }
    if !row.is_empty() {
        rows.push(row);
    }
    rows
}

/// One end of a covalent connection, identified as the coordinate reader
/// identifies residues: author chain and sequence number (label ones when
/// the author ones are absent), insertion code and residue name.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub(crate) struct AtomReference {
    pub chain: String,
    pub sequence: isize,
    pub insertion_code: Option<String>,
    pub residue: String,
    pub atom: String,
}

/// The covalent connections a file records between residues.
#[derive(Debug, Clone, Default)]
pub(crate) struct ConnectionRecords {
    /// Whether the file records any covalent connection between residues
    /// (`_struct_conn` of type `covale*`, `LINK`), whatever their partners.
    pub any: bool,
    /// Those within the asymmetric unit whose partners could be read.
    pub covalent: Vec<(AtomReference, AtomReference)>,
}

fn same_asymmetric_unit(symmetry: Option<&str>) -> bool {
    symmetry.is_none_or(|operator| matches!(operator.trim(), "1_555" | "1555" | ""))
}

/// The covalent connections of `_struct_conn` (`conn_type_id` `covale`,
/// `covale_base`, `covale_phosphate`, `covale_sugar`); metal coordination,
/// hydrogen bonds and disulfides are left out.
pub(crate) fn mmcif_connections(text: &str) -> ConnectionRecords {
    let mut records = ConnectionRecords::default();
    for row in category_rows(text, "_struct_conn") {
        let covalent = row
            .get("conn_type_id")
            .is_some_and(|kind| kind.to_ascii_lowercase().starts_with("covale"));
        if !covalent {
            continue;
        }
        records.any = true;
        let partner = |number: &str| {
            let get = |auth: &str, label: &str| {
                row.get(&format!("ptnr{number}_{auth}"))
                    .or_else(|| row.get(&format!("ptnr{number}_{label}")))
            };
            let symmetry = row.get(&format!("ptnr{number}_symmetry"));
            if !same_asymmetric_unit(symmetry.map(String::as_str)) {
                return None;
            }
            Some(AtomReference {
                chain: get("auth_asym_id", "label_asym_id")?.clone(),
                sequence: get("auth_seq_id", "label_seq_id")?.parse().ok()?,
                insertion_code: row.get(&format!("pdbx_ptnr{number}_pdb_ins_code")).cloned(),
                residue: get("label_comp_id", "auth_comp_id")?.to_ascii_uppercase(),
                atom: get("label_atom_id", "auth_atom_id")?.clone(),
            })
        };
        if let (Some(first), Some(second)) = (partner("1"), partner("2")) {
            records.covalent.push((first, second));
        }
    }
    records
}

/// The connections of PDB `LINK` records. The records do not say which are
/// covalent; a metal coordination is told apart by the elements later.
pub(crate) fn pdb_connections(text: &str) -> ConnectionRecords {
    let mut records = ConnectionRecords::default();
    for line in text.lines().filter(|line| line.starts_with("LINK  ")) {
        records.any = true;
        let field = |range: std::ops::Range<usize>| {
            line.get(range.start..range.end.min(line.len()))
                .map(str::trim)
                .unwrap_or("")
        };
        let partner = |offset: usize, symmetry: std::ops::Range<usize>| {
            let symmetry = field(symmetry);
            if !same_asymmetric_unit(Some(symmetry)) {
                return None;
            }
            let insertion = field(offset + 14..offset + 15);
            Some(AtomReference {
                atom: field(offset..offset + 4).to_string(),
                residue: field(offset + 5..offset + 8).to_ascii_uppercase(),
                chain: field(offset + 9..offset + 10).to_string(),
                sequence: field(offset + 10..offset + 14).parse().ok()?,
                insertion_code: (!insertion.is_empty()).then(|| insertion.to_string()),
            })
        };
        if let (Some(first), Some(second)) = (partner(12, 59..65), partner(42, 66..72)) {
            records.covalent.push((first, second));
        }
    }
    records
}

/// Crystal items pdbtbx reads and the extractor does not need: the
/// space-group items, which pdbtbx compares with each other, and the cell,
/// whose unknown values (`?`) stop pdbtbx.
const CRYSTAL_ITEMS: [&str; 12] = [
    "_symmetry.space_group_name_h-m",
    "_symmetry.space_group_name_hall",
    "_symmetry.int_tables_number",
    "_space_group.name_h-m_alt",
    "_space_group.name_hall",
    "_space_group.it_number",
    "_cell.length_a",
    "_cell.length_b",
    "_cell.length_c",
    "_cell.angle_alpha",
    "_cell.angle_beta",
    "_cell.angle_gamma",
];

/// Whether a crystal item has the value `?` or `.` on its own line (pdbtbx
/// panics on one).
pub(crate) fn has_unknown_crystal_item(text: &str) -> bool {
    text.lines().any(|line| {
        let mut words = line.split_ascii_whitespace();
        let (Some(tag), Some(value), None) = (words.next(), words.next(), words.next()) else {
            return false;
        };
        matches!(value, "?" | ".") && CRYSTAL_ITEMS.contains(&tag.to_ascii_lowercase().as_str())
    })
}

/// `text` without its crystal items (single items, as pdbtbx reads them),
/// or `None` when it has none.
pub(crate) fn without_crystal_items(text: &str) -> Option<String> {
    let mut removed: Vec<(usize, usize)> = Vec::new();
    let mut tokens = Tokens::new(text, 0).peekable();
    let mut in_loop_header = false;
    while let Some(token) = tokens.next() {
        if !token.quoted && token.text.eq_ignore_ascii_case("loop_") {
            in_loop_header = true;
            continue;
        }
        if !token.is_tag() {
            in_loop_header = false;
            continue;
        }
        if in_loop_header {
            continue;
        }
        let Some(value) = tokens.next_if(|next| !next.is_tag() && !next.is_reserved()) else {
            continue;
        };
        if CRYSTAL_ITEMS.contains(&token.text.to_ascii_lowercase().as_str()) {
            removed.push((token.start, value.end));
        }
    }
    if removed.is_empty() {
        return None;
    }
    let mut kept = String::with_capacity(text.len());
    let mut cursor = 0;
    for (start, end) in removed {
        kept.push_str(&text[cursor..start]);
        cursor = end;
    }
    kept.push_str(&text[cursor..]);
    Some(kept)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tokens_follow_cif_quoting() {
        let text = "data_x\n_a.b 'it''s' # note\n_a.c \"O3'\"\n;\nline one\nline two\n;\n_a.d ?\n";
        let tokens = Tokens::new(text, 0)
            .map(|token| (token.text, token.quoted))
            .collect::<Vec<_>>();
        assert_eq!(
            tokens,
            vec![
                ("data_x", false),
                ("_a.b", false),
                ("it''s", true),
                ("_a.c", false),
                ("O3'", true),
                ("\nline one\nline two", true),
                ("_a.d", false),
                ("?", false),
            ]
        );
    }

    const STRUCT_CONN: &str = "data_1EO8
#
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.pdbx_ptnr1_PDB_ins_code
_struct_conn.ptnr1_symmetry
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.pdbx_ptnr2_PDB_ins_code
_struct_conn.ptnr1_auth_asym_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.ptnr2_auth_asym_id
_struct_conn.ptnr2_auth_seq_id
_struct_conn.ptnr2_symmetry
covale1 covale A ASN 81 ND2 ? 1_555 G NAG . C1 ? A 81 A 440 1_555
metalc1 metalc A ASP 12 OD1 ? 1_555 J CA . CA ? A 12 A 501 1_555
covale2 covale E NAG . O4 ? 1_555 E NAG . C1 ? C 1 C 2 1_555
covale3 covale A SER 5 OG ? 1_555 K MAN . C1 ? A 5 A 601 2_655
#
_struct_conn_type.id covale
";

    #[test]
    fn struct_conn_rows_give_covalent_connections_by_author_identifiers() {
        let records = mmcif_connections(STRUCT_CONN);
        assert!(records.any);
        let pairs = records
            .covalent
            .iter()
            .map(|(first, second)| {
                (
                    format!(
                        "{}/{}/{}/{}",
                        first.chain, first.residue, first.sequence, first.atom
                    ),
                    format!(
                        "{}/{}/{}/{}",
                        second.chain, second.residue, second.sequence, second.atom
                    ),
                )
            })
            .collect::<Vec<_>>();
        // The metal coordination and the bond to a symmetry mate are left out.
        assert_eq!(
            pairs,
            vec![
                ("A/ASN/81/ND2".to_string(), "A/NAG/440/C1".to_string()),
                ("C/NAG/1/O4".to_string(), "C/NAG/2/C1".to_string()),
            ]
        );
    }

    #[test]
    fn a_single_struct_conn_may_be_written_as_items() {
        let text = "data_x\n_struct_conn.id covale1\n_struct_conn.conn_type_id covale\n_struct_conn.ptnr1_auth_asym_id A\n_struct_conn.ptnr1_label_comp_id TRP\n_struct_conn.ptnr1_auth_seq_id 7\n_struct_conn.ptnr1_label_atom_id CD1\n_struct_conn.ptnr2_auth_asym_id A\n_struct_conn.ptnr2_label_comp_id MAN\n_struct_conn.ptnr2_auth_seq_id 101\n_struct_conn.ptnr2_label_atom_id C1\n#\n_cell.length_a 1\n";
        let records = mmcif_connections(text);
        assert_eq!(records.covalent.len(), 1);
        assert_eq!(records.covalent[0].0.atom, "CD1");
        assert_eq!(records.covalent[0].1.residue, "MAN");
        assert_eq!(records.covalent[0].1.sequence, 101);
    }

    #[test]
    fn link_records_name_both_partners() {
        let text = "LINK         ND2 ASN A  81                 C1  NAG A 440     1555   1555  1.46  \nLINK         O4  NAG C   1                 C1  NAG C   2     1555   1555  1.42  \nLINK         OG  SER A   5                 C1  MAN A 601     1555   2655  1.42  \n";
        let records = pdb_connections(text);
        assert!(records.any);
        assert_eq!(records.covalent.len(), 2);
        assert_eq!(
            records.covalent[0],
            (
                AtomReference {
                    chain: "A".into(),
                    sequence: 81,
                    insertion_code: None,
                    residue: "ASN".into(),
                    atom: "ND2".into(),
                },
                AtomReference {
                    chain: "A".into(),
                    sequence: 440,
                    insertion_code: None,
                    residue: "NAG".into(),
                    atom: "C1".into(),
                }
            )
        );
    }

    #[test]
    fn crystal_items_are_removed_and_nothing_else() {
        let text = "data_31NP\n_cell.length_a 10.0\n_cell.volume 1.0\n_symmetry.entry_id 31NP\n_symmetry.space_group_name_H-M 'I 1 2 1'\n_symmetry.space_group_name_Hall\n'C 2y (x,y,-x+z)'\n_symmetry.Int_Tables_number 5\nloop_\n_space_group_symop.id\n_space_group_symop.operation_xyz\n1 x,y,z\n";
        assert!(!has_unknown_crystal_item(text));
        assert_eq!(
            without_crystal_items(text).unwrap(),
            "data_31NP\n\n_cell.volume 1.0\n_symmetry.entry_id 31NP\n\n\n\nloop_\n_space_group_symop.id\n_space_group_symop.operation_xyz\n1 x,y,z\n"
        );
        assert_eq!(without_crystal_items("data_x\n_cell.volume 1\n"), None);
        assert!(has_unknown_crystal_item(
            "data_x\n_symmetry.Int_Tables_number   ?\n"
        ));
    }
}
