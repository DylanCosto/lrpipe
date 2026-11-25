#!/usr/bin/env bash
# Regenerates crates/snv-select/tests/data/*/expected with Clair3's own SelectQual,
# SelectHetSnp and SelectCandidates, run the way the pipeline runs them (input through `pigz -fdc`, under PyPy).
#
# Each case directory holds `input.vcf` and `runs.txt`, one run per line:
#   phase <var_pct_phasing>                     SelectQual --phase
#   full <var_pct_full> <ref_pct_full> <vcf_fn>  SelectQual (vcf_fn relative to the case dir)
#   het <ctgName|->                             SelectHetSnp, with the phase_qual of the latest
#                                               phase run (none before one: the fallback path)
#   cand <ctgName|-> <var_pct_full> <ref_pct_full> <platform> <split_bed_size> <qual>
#                                               SelectCandidates; <qual> is `latest` (the qual
#                                               of the latest full run), `none` (the fallback
#                                               path) or a file in the case dir
# Run N's outputs, stdout and stderr go to expected/N-<kind>/. In FULL_ALN_FILE_*, the split
# folder is replaced by @SPLIT_FOLDER@.
set -euo pipefail
ENV=${TORCH_ENV:-/data2/bioinformatics/envs/clair3-t213}/bin
DATA=$(cd "$(dirname "$0")/../crates/snv-select/tests/data" && pwd)
clair3() { PATH="$ENV:$PATH" pypy3 "$ENV/clair3.py" "$@"; }

for case in "$DATA"/*/; do
    case=${case%/}
    rm -rf "$case/expected"
    mkdir -p "$case/expected"
    cd "$case"
    n=0
    phase_qual=
    full_qual=
    while read -r kind a b c d e f; do
        n=$((n + 1))
        out="$case/expected/$n-$kind"
        mkdir -p "$out"
        case $kind in
        phase)
            pigz -fdc -p 2 input.vcf | clair3 SelectQual --phase --output_fn "$out" --var_pct_phasing "$a" \
                >"$out/stdout" 2>"$out/stderr"
            phase_qual="$out/phase_qual"
            ;;
        full)
            pigz -fdc -p 2 input.vcf | clair3 SelectQual --output_fn "$out" --var_pct_full "$a" \
                --ref_pct_full "$b" --platform ont --vcf_fn "$c" >"$out/stdout" 2>"$out/stderr"
            full_qual="$out/qual"
            ;;
        het)
            [[ -n $phase_qual ]] && cp "$phase_qual" "$out/phase_qual"
            ctg=()
            [[ $a != - ]] && ctg=(--ctgName "$a")
            clair3 SelectHetSnp --vcf_fn input.vcf --split_folder "$out" "${ctg[@]}" \
                >"$out/stdout" 2>"$out/stderr"
            rm -f "$out/phase_qual"
            ;;
        cand)
            case $f in
            latest) cp "$full_qual" "$out/qual" ;;
            none) ;;
            *) cp "$f" "$out/qual" ;;
            esac
            ctg=()
            [[ $a != - ]] && ctg=(--ctgName "$a")
            clair3 SelectCandidates --pileup_vcf_fn input.vcf --split_folder "$out" "${ctg[@]}" \
                --var_pct_full "$b" --ref_pct_full "$c" --platform "$d" --split_bed_size "$e" \
                >"$out/stdout" 2>"$out/stderr"
            rm -f "$out/qual"
            for list in "$out"/FULL_ALN_FILE_*; do
                [[ -e $list ]] && sed -i.bak "s#^$out/#@SPLIT_FOLDER@/#" "$list" && rm "$list.bak"
            done
            ;;
        *)
            echo "unknown run kind $kind in $case/runs.txt" >&2
            exit 1
            ;;
        esac
    done <runs.txt
    echo "$(basename "$case"): $n runs"
done
