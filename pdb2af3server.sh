#!/usr/bin/env bash

set -euo pipefail

DIR=""
TARGET=""
CHAIN=""
RANGE=""
COPY_COUNT=""

usage() {
    echo "Usage: ./pdb2af3server.sh --dir <pdb directory> --target <target protein pdb> --chain <chain id> --range <residue range> --copy <count>"
    exit 1
}

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --dir) DIR="${2:-}"; shift 2 ;;
        --target) TARGET="${2:-}"; shift 2 ;;
        --chain) CHAIN="${2:-}"; shift 2 ;;
        --range) RANGE="${2:-}"; shift 2 ;;
        --copy) COPY_COUNT="${2:-}"; shift 2 ;;
        *) echo "Unknown parameter: $1"; usage ;;
    esac
done

if [[ -z "$DIR" || -z "$TARGET" || -z "$CHAIN" || -z "$RANGE" || -z "$COPY_COUNT" ]]; then
    usage
fi

aa3to1_awk='
function aa3to1(res) {
    if (res=="ALA") return "A"; else if (res=="ARG") return "R";
    else if (res=="ASN") return "N"; else if (res=="ASP") return "D";
    else if (res=="CYS") return "C"; else if (res=="GLU") return "E";
    else if (res=="GLN") return "Q"; else if (res=="GLY") return "G";
    else if (res=="HIS") return "H"; else if (res=="ILE") return "I";
    else if (res=="LEU") return "L"; else if (res=="LYS") return "K";
    else if (res=="MET") return "M"; else if (res=="PHE") return "F";
    else if (res=="PRO") return "P"; else if (res=="SER") return "S";
    else if (res=="THR") return "T"; else if (res=="TRP") return "W";
    else if (res=="TYR") return "Y"; else if (res=="VAL") return "V";
    else return "X";
}
function trim(s){gsub(/^[ \t]+|[ \t]+$/,"",s); return s}
function build_allowed(range_spec){
    n=split(range_spec,parts,",");
    for(i=1;i<=n;i++){
        token=trim(parts[i]);
        if(token~/^[0-9]+-[0-9]+$/){
            split(token,b,"-");
            for(j=b[1];j<=b[2];j++) allowed[j]=1;
        } else if(token~/^[0-9]+$/){
            allowed[token]=1;
        }
    }
}
'

pdb_seq_range() {
    awk -v chain="$2" -v range="$3" '
    '"$aa3to1_awk"'
    BEGIN{build_allowed(range)}
    /^ATOM/{
        if(substr($0,13,4)!~/CA/) next
        if(substr($0,22,1)!=chain) next
        res=substr($0,23,4)+0
        if(!(res in allowed)) next
        if(!(res in seen)){
            seq=seq aa3to1(substr($0,18,3))
            seen[res]=1
        }
    }
    END{print seq}' "$1"
}

pdb_seq_chain() {
    awk -v chain="$2" '
    '"$aa3to1_awk"'
    /^ATOM/{
        if(substr($0,13,4)!~/CA/) next
        if(substr($0,22,1)!=chain) next
        res=substr($0,23,4)
        if(!(res in seen)){
            seq=seq aa3to1(substr($0,18,3))
            seen[res]=1
        }
    }
    END{print seq}' "$1"
}

TARGET_SEQ=$(pdb_seq_range "$TARGET" "$CHAIN" "$RANGE")

mkdir -p "$DIR/af3_jsons"

for pdb in "$DIR"/*.pdb; do
    [[ "$(realpath "$pdb")" == "$(realpath "$TARGET")" ]] && continue

    CAND_SEQ=$(pdb_seq_chain "$pdb" "$CHAIN")
    [[ -z "$CAND_SEQ" ]] && continue

    name=$(basename "$pdb" .pdb)
    out="$DIR/af3_jsons/${name}.json"

    # *** 핵심 수정: 변수 명시적으로 분리 ***
    TARGET_BLOCK_SEQ="$TARGET_SEQ"
    CANDIDATE_BLOCK_SEQ="$CAND_SEQ"

    cat > "$out" <<EOF
[
  {
    "name": "$name",
    "modelSeeds": [],
    "sequences": [
      {
        "proteinChain": {
          "sequence": "$TARGET_BLOCK_SEQ",
          "count": $COPY_COUNT
        }
      },
      {
        "proteinChain": {
          "sequence": "$CANDIDATE_BLOCK_SEQ",
          "count": 1
        }
      }
    ]
  }
]
EOF

done
