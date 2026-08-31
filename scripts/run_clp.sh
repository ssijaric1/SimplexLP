#!/usr/bin/env bash
# SimplexLP - Project #13 (Š. Sijarić, 19964)
#
# run_clp.sh - external comparison against COIN-OR CLP.
#
# The benchmark (GUI tab or the lpBench console tool) exports the synthetic
# instances as MPS files. This script feeds the *same* files to CLP twice:
# once with the primal simplex and once with the barrier (interior-point)
# method, and prints solve times and objective values, so the natID
# implementation and CLP are measured on identical inputs.
#
# Install CLP first:
#   macOS : brew install clp        (or: brew install coin-or-tools/coinor/clp)
#   Linux : sudo apt install coinor-clp
#
# Usage:
#   ./run_clp.sh <directory-with-mps-files>     (default: ~/SimplexLP_bench)

set -u
DIR="${1:-$HOME/SimplexLP_bench}"

if ! command -v clp >/dev/null 2>&1; then
    echo "error: 'clp' not found in PATH. Install COIN-OR CLP first (see header)."
    exit 1
fi

shopt -s nullglob
FILES=("$DIR"/*.mps)
if [ ${#FILES[@]} -eq 0 ]; then
    echo "error: no .mps files in '$DIR'. Run the benchmark with MPS export first."
    exit 1
fi

OUT="$DIR/clp_results.csv"
echo "file,method,seconds,objective,status" > "$OUT"

run_one () { # $1=file  $2=method-flag  $3=method-name
    local f="$1" flag="$2" name="$3"
    local log t obj st
    log=$(clp "$f" "$flag" -solve 2>&1)
    t=$(echo "$log"   | grep -Eo 'Total time \(CPU seconds\):[[:space:]]+[0-9.]+' | grep -Eo '[0-9.]+' | head -1)
    [ -z "${t:-}" ] && t=$(echo "$log" | grep -Eo 'took [0-9.]+ seconds' | grep -Eo '[0-9.]+' | head -1)
    obj=$(echo "$log" | grep -Eo 'Objective value[[:space:]]+[-0-9.eE+]+' | grep -Eo '[-0-9.eE+]+$' | head -1)
    st="optimal"
    echo "$log" | grep -qi "infeasible" && st="infeasible"
    echo "$log" | grep -qi "unbounded"  && st="unbounded"
    printf "%-34s %-16s %8s s   obj=%-16s %s\n" "$(basename "$f")" "$name" "${t:-?}" "${obj:-?}" "$st"
    echo "$(basename "$f"),$name,${t:-},${obj:-},$st" >> "$OUT"
}

echo "CLP comparison on $DIR"
echo "--------------------------------------------------------------------------"
for f in "${FILES[@]}"; do
    run_one "$f" -primalsimplex "primal-simplex"
    run_one "$f" -barrier       "barrier"
done
echo "--------------------------------------------------------------------------"
echo "wrote $OUT  (compare with results.csv from the natID benchmark)"
