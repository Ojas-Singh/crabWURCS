#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_dir=$(cd "${script_dir}/.." && pwd)
output_path="${repo_dir}/crabwurcs/data/glycoshape_derived_notations.tsv"
temp_path=$(mktemp "${TMPDIR:-/tmp}/crabwurcs-derived.XXXXXX")
trap 'rm -f "${temp_path}"' EXIT

cd "${repo_dir}"
cargo build --quiet -p crabwurcs-cli

jq -c '
  to_entries[].value.archetype
  | select((.smiles // "") != "" and .wurcs == null)
  | {iupac, iupac_extended, glycam, smiles}
' GLYCOSHAPE.json | while IFS= read -r record; do
    iupac=$(jq -r '.iupac // ""' <<<"${record}")
    extended=$(jq -r '.iupac_extended // ""' <<<"${record}")
    glycam=$(jq -r '.glycam // ""' <<<"${record}")
    smiles=$(jq -r '.smiles' <<<"${record}")
    wurcs=$(printf '%s' "${iupac}" \
      | target/debug/crabwurcs convert --from iupac-condensed --to wurcs)
    printf '%s\t%s\t%s\t%s\t%s\n' \
      "${wurcs}" "${iupac}" "${extended}" "${glycam}" "${smiles}" \
      >>"${temp_path}"
done

record_count=$(wc -l <"${temp_path}" | tr -d ' ')
if [[ "${record_count}" != "100" ]]; then
  printf 'expected 100 derived molecular records, generated %s\n' "${record_count}" >&2
  exit 1
fi

mv "${temp_path}" "${output_path}"
trap - EXIT
printf 'generated %s\n' "${output_path}"
cargo run --quiet -p crabwurcs-mol --example regenerate_canonical_index
