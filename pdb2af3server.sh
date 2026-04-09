#!/usr/bin/env bash

set -euo pipefail

DIR=""
TARGET=""
TARGET_CHAIN=""
CANDIDATE_CHAIN=""
RANGE=""
COPY_COUNT=""

usage() {
    echo "Usage: ./pdb2af3server.sh --dir <pdb directory> --target <target protein pdb> --target-chain <chain id> --candidate-chain <chain id> --range <residue range> --copy <count>"
    echo
    echo "Examples:"
    echo "  ./pdb2af3server.sh --dir ./pdbs --target target.pdb --target-chain A --candidate-chain B --range 10-120 --copy 2"
    echo "  ./pdb2af3server.sh --dir ./pdbs --target target.pdb --target-chain A --candidate-chain C --range 10-50,80-120,150 --copy 3"
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
        --target-chain)
            TARGET_CHAIN="${2:-}"
            shift 2
            ;;
        --candidate-chain)
            CANDIDATE_CHAIN="${2:-}"
            shift 2
            ;;
        --range)
            RANGE="${2:-}"
            shift 2
            ;;
        --copy)
            COPY_COUNT="${2:-}"
            shift 2
            ;;
        *)
            echo "Unknown parameter passed: $1"
            usage
            ;;
    esac
done

if [[ -z "$DIR" || -z "$TARGET" || -z "$TARGET_CHAIN" || -z "$CANDIDATE_CHAIN" || -z "$RANGE" || -z "$COPY_COUNT" ]]; then
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

if [[ ${#TARGET_CHAIN} -ne 1 ]]; then
    echo "Error: --target-chain must be a single chain identifier."
    exit 1
fi

if [[ ${#CANDIDATE_CHAIN} -ne 1 ]]; then
    echo "Error: --candidate-chain must be a single chain identifier."
    exit 1
fi

if ! [[ "$COPY_COUNT" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: --copy must be a positive integer."
    exit 1
fi

if ! [[ "$RANGE" =~ ^[0-9,[:space:]-]+$ ]]; then
    echo "Error: --range must contain only digits, commas, hyphens, and spaces."
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

function trim(s) {
    gsub(/^[ \t]+|[ \t]+$/, "", s)
    return s
}

function build_allowed_ranges(range_spec,    i, n, token, bounds, start, end, tmp, j) {
    n = split(range_spec, parts, ",")
    for (i = 1; i <= n; i++) {
        token = trim(parts[i])
        if (token == "") continue

        if (token ~ /^[0-9]+-[0-9]+$/) {
            split(token, bounds, "-")
            start = bounds[1] + 0
            end = bounds[2] + 0
            if (start > end) {
                tmp = start
                start = end
                end = tmp
            }
            for (j = start; j <= end; j++) {
                allowed[j] = 1
            }
        } else if (token ~ /^[0-9]+$/) {
            allowed[token + 0] = 1
        } else {
            invalid_range = token
            return 0
        }
    }
    return 1
}
'

pdb_to_seq_chain_range() {
    local pdb_file="$1"
    local chain_id="$2"
    local range_spec="$3"
    awk -v chain="$chain_id" -v range_spec="$range_spec" '
    '"$aa3to1_awk"'
    BEGIN {
        ok = build_allowed_ranges(range_spec)
        if (!ok) {
            printf("Error: invalid --range token: %s\n", invalid_range) > "/dev/stderr"
            exit 2
        }
    }
    /^ATOM/ {
        atom_name = substr($0, 13, 4)
        gsub(/ /, "", atom_name)
        pdb_chain = substr($0, 22, 1)
        if (atom_name != "CA" || pdb_chain != chain) next

        resname = substr($0, 18, 3)
        resseq_raw = substr($0, 23, 4)
        icode = substr($0, 27, 1)
        resseq = trim(resseq_raw) + 0
        resid = resseq "_" icode

        if (!(resseq in allowed)) next

        if (!(resid in seen)) {
            seq = seq aa3to1(resname)
            seen[resid] = 1
        }
    }
    END { print seq }' "$pdb_file"
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
TARGET_SEQ="$(pdb_to_seq_chain_range "$TARGET" "$TARGET_CHAIN" "$RANGE")"

if [[ -z "$TARGET_SEQ" ]]; then
    echo "Error: failed to extract sequence from target PDB for chain '$TARGET_CHAIN' and range '$RANGE': $TARGET"
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

    candidate_seq="$(pdb_to_seq_chain "$pdb" "$CANDIDATE_CHAIN")"
    if [[ -z "$candidate_seq" ]]; then
        continue
    fi

    job_name="$(basename "$pdb" .pdb)"
    output_json="$output_dir/${job_name}.json"

    target_block_seq="$TARGET_SEQ"
    candidate_block_seq="$candidate_seq"

    cat > "$output_json" <<EOF_JSON
[
  {
    "name": "$job_name",
    "modelSeeds": [],
    "sequences": [
      {
        "proteinChain": {
          "sequence": "$target_block_seq",
          "count": $COPY_COUNT
        }
      },
      {
        "proteinChain": {
          "sequence": "$candidate_block_seq",
          "count": 1
        }
      }
    ]
  }
]
EOF_JSON

    job_count=$((job_count + 1))
done

if [[ $job_count -eq 0 ]]; then
    echo "Error: no candidate PDB with chain '$CANDIDATE_CHAIN' was found in $DIR"
    exit 1
fi

echo "AF3 JSON files generated in: $output_dir"
echo "Target chain: $TARGET_CHAIN"
echo "Target residue range: $RANGE"
echo "Target copy count: $COPY_COUNT"
echo "Candidate chain: $CANDIDATE_CHAIN"
echo "JSON files written: $job_count (1 per PDB)"
