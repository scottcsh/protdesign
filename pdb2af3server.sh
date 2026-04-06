#!/usr/bin/env bash

set -euo pipefail

DIR=""
TARGET=""
CHAIN=""

usage() {
    echo "Usage: ./pdb2af3server.sh --dir <pdb directory> --target <target protein pdb> --chain <chain id>"
    exit 1
}

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --dir)
            DIR="${2:-}"
            shift 2
            ;;
        --target)
            TARGET="${2:-}"
            shift 2
            ;;
        --chain)
            CHAIN="${2:-}"
            shift 2
            ;;
        *)
            echo "Unknown parameter passed: $1"
            usage
            ;;
    esac
done

if [[ -z "$DIR" || -z "$TARGET" || -z "$CHAIN" ]]; then
    usage
fi

if [[ ! -d "$DIR" ]]; then
    echo "Error: directory not found: $DIR"
    exit 1
fi

if [[ ! -f "$TARGET" ]]; then
    echo "Error: target pdb not found: $TARGET"
    exit 1
fi

if [[ ${#CHAIN} -ne 1 ]]; then
    echo "Error: --chain must be a single chain identifier."
    exit 1
fi

aa3to1_awk='
function aa3to1(res) {
    if (res == "ALA") return "A"
    else if (res == "ARG") return "R"
    else if (res == "ASN") return "N"
    else if (res == "ASP") return "D"
    else if (res == "CYS") return "C"
    else if (res == "GLU") return "E"
    else if (res == "GLN") return "Q"
    else if (res == "GLY") return "G"
    else if (res == "HIS") return "H"
    else if (res == "ILE") return "I"
    else if (res == "LEU") return "L"
    else if (res == "LYS") return "K"
    else if (res == "MET") return "M"
    else if (res == "PHE") return "F"
    else if (res == "PRO") return "P"
    else if (res == "SER") return "S"
    else if (res == "THR") return "T"
    else if (res == "TRP") return "W"
    else if (res == "TYR") return "Y"
    else if (res == "VAL") return "V"
    else return "X"
}
'

pdb_to_seq_all_chains() {
    local pdb_file="$1"
    awk '
    '"$aa3to1_awk"'
    /^ATOM/ {
        atom_name = substr($0, 13, 4)
        gsub(/ /, "", atom_name)
        if (atom_name != "CA") next

        chain = substr($0, 22, 1)
        resname = substr($0, 18, 3)
        resseq = substr($0, 23, 4)
        icode = substr($0, 27, 1)
        resid = chain "_" resseq icode

        if (!(resid in seen)) {
            if (!(chain in chain_seen)) {
                chain_order[++nchains] = chain
                chain_seen[chain] = 1
            }
            seqs[chain] = seqs[chain] aa3to1(resname)
            seen[resid] = 1
        }
    }
    END {
        out = ""
        for (i = 1; i <= nchains; i++) {
            c = chain_order[i]
            if (i > 1) out = out ":"
            out = out seqs[c]
        }
        print out
    }' "$pdb_file"
}

pdb_to_seq_chain() {
    local pdb_file="$1"
    local chain_id="$2"
    awk -v chain="$chain_id" '
    '"$aa3to1_awk"'
    /^ATOM/ {
        atom_name = substr($0, 13, 4)
        gsub(/ /, "", atom_name)
        pdb_chain = substr($0, 22, 1)
        if (atom_name != "CA" || pdb_chain != chain) next

        resname = substr($0, 18, 3)
        resseq = substr($0, 23, 4)
        icode = substr($0, 27, 1)
        resid = resseq icode

        if (!(resid in seen)) {
            seq = seq aa3to1(resname)
            seen[resid] = 1
        }
    }
    END { print seq }' "$pdb_file"
}

TARGET_ABS="$(realpath "$TARGET")"
TARGET_SEQ="$(pdb_to_seq_all_chains "$TARGET")"

if [[ -z "$TARGET_SEQ" ]]; then
    echo "Error: failed to extract sequence from target PDB: $TARGET"
    exit 1
fi

output_dir="$DIR/af3_jsons"
mkdir -p "$output_dir"

job_count=0

for pdb in "$DIR"/*.pdb; do
    [[ -e "$pdb" ]] || continue

    if [[ "$(realpath "$pdb")" == "$TARGET_ABS" ]]; then
        continue
    fi

    candidate_seq="$(pdb_to_seq_chain "$pdb" "$CHAIN")"
    if [[ -z "$candidate_seq" ]]; then
        continue
    fi

    job_name="$(basename "$pdb" .pdb)"
    output_json="$output_dir/${job_name}.json"

    cat > "$output_json" <<EOF
[
  {
    "name": "$job_name",
    "modelSeeds": [],
    "sequences": [
      {
        "proteinChain": {
          "sequence": "$TARGET_SEQ",
          "count": 1
        }
      },
      {
        "proteinChain": {
          "sequence": "$candidate_seq",
          "count": 1
        }
      }
    ]
  }
]
EOF

    job_count=$((job_count + 1))
done

if [[ $job_count -eq 0 ]]; then
    echo "Error: no candidate PDB with chain '$CHAIN' was found in $DIR"
    exit 1
fi

echo "AF3 JSON files generated in: $output_dir"
echo "Target: all chains kept"
echo "Candidate chain: $CHAIN"
echo "JSON files written: $job_count (1 per PDB)"
