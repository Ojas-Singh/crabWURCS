use std::collections::BTreeMap;
use std::fs;
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let workspace = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..");
    let data_dir = workspace.join("crabwurcs/data");
    let mut index = BTreeMap::<String, String>::new();
    let mut records = 0usize;

    for filename in [
        "glycoshape_notations.tsv",
        "glycoshape_derived_notations.tsv",
    ] {
        let input = fs::read_to_string(data_dir.join(filename))?;
        for line in input.lines() {
            let fields = line.split('\t').collect::<Vec<_>>();
            if fields.len() != 5 {
                return Err(
                    format!("{filename} contains a record with {} fields", fields.len()).into(),
                );
            }
            let molecule = chematic::smiles::parse(fields[4])?;
            let canonical = chematic::smiles::canonical_smiles(&molecule);
            if let Some(previous) = index.insert(canonical.clone(), fields[0].to_owned()) {
                if previous != fields[0] {
                    return Err(format!(
                        "canonical molecular collision maps to two WURCS:\n{previous}\n{}\n{canonical}",
                        fields[0]
                    )
                    .into());
                }
            }
            records += 1;
        }
    }

    if records != 938 {
        return Err(format!("expected 938 molecular records, found {records}").into());
    }

    let output = index
        .into_iter()
        .map(|(canonical, wurcs)| format!("{wurcs}\t{canonical}\n"))
        .collect::<String>();
    let path = data_dir.join("glycoshape_canonical_smiles.tsv");
    fs::write(&path, output)?;
    println!("generated {}", path.display());
    Ok(())
}
