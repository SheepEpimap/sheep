#!/usr/bin/env bash
set -euo pipefail

: "${GWAS_BASE:?Set GWAS_BASE in the paths configuration}"
: "${OLD_GWAS_BASE:?Set OLD_GWAS_BASE in the paths configuration}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p \
    "$GWAS_BASE/beds" \
    "$GWAS_BASE/control_beds" \
    "$GWAS_BASE/ldcts" \
    "$GWAS_BASE/logs" \
    "$GWAS_BASE/manifests" \
    "$GWAS_BASE/outputs" \
    "$GWAS_BASE/summary" \
    "$GWAS_BASE/figures" \
    "$GWAS_BASE/scripts"

echo "[INFO] GWAS_BASE=$GWAS_BASE"
echo "[INFO] OLD_GWAS_BASE=$OLD_GWAS_BASE"
echo "[INFO] SCRIPT_DIR=$SCRIPT_DIR"

build_control_beds() {
    local bed_dir="$GWAS_BASE/beds"
    local control_dir="$GWAS_BASE/control_beds"
    local out
    local inputs=()
    local t
    local c
    local missing=()

    local tissues_with_muscle=(Adipose Colon Cortex Heart Liver Lung Muscle Ovary Sintest Spleen Stomach)
    local tissues_no_muscle=(Adipose Colon Cortex Heart Liver Lung Ovary Sintest Spleen Stomach)

    if ! compgen -G "$bed_dir/*.bed" > /dev/null; then
        echo "[WARN] cannot build control BEDs because no BED files exist under $bed_dir"
        return 0
    fi

    command -v bedtools >/dev/null 2>&1 || {
        echo "[ERROR] bedtools is required to build control BEDs" >&2
        exit 1
    }

    mkdir -p "$control_dir"

    write_union() {
        out="$1"
        shift
        inputs=("$@")
        missing=()

        for f in "${inputs[@]}"; do
            [[ -s "$f" ]] || missing+=("$f")
        done

        if [[ "${#missing[@]}" -gt 0 ]]; then
            echo "[ERROR] missing input BEDs for $out:" >&2
            printf '  %s\n' "${missing[@]}" >&2
            exit 1
        fi

        cat "${inputs[@]}" \
            | LC_ALL=C sort -k1,1V -k2,2n -k3,3n \
            | bedtools merge -i - > "$out"
        echo "[INFO] wrote control BED: $out"
    }

    inputs=()
    for t in "${tissues_no_muscle[@]}"; do
        inputs+=("$bed_dir/sfCRE_ts_${t}.bed")
    done
    write_union "$control_dir/sfCRE_ts_union.bed" "${inputs[@]}"

    for c in sdCRE soCRE ssCRE; do
        inputs=()
        for t in "${tissues_with_muscle[@]}"; do
            inputs+=("$bed_dir/${c}_ts_${t}.bed")
        done
        write_union "$control_dir/${c}_ts_union.bed" "${inputs[@]}"
    done

    inputs=()
    for t in "${tissues_no_muscle[@]}"; do
        inputs+=("$bed_dir/sf_proj_ts_${t}.bed")
    done
    write_union "$control_dir/sf_proj_ts_union.bed" "${inputs[@]}"
}

validate_bed_manifest() {
    local mf="$1"
    local tmp="$GWAS_BASE/manifests/.$(basename "$mf").bed_paths.tmp"
    local missing=0
    local bed
    local resolved

    [[ -s "$mf" ]] || return 0

    awk -v base="$GWAS_BASE" 'BEGIN{FS="\t"} NF >= 3 {
        bed=$3
        sub(/\r$/, "", bed)
        if (tolower(bed) == "bed") next
        if (bed !~ /^\//) bed=base "/" bed
        print bed
    }' "$mf" | sort -u > "$tmp"

    while IFS= read -r bed; do
        [[ -n "$bed" ]] || continue
        resolved="$bed"
        if [[ ! -s "$resolved" ]]; then
            echo "[ERROR] manifest references missing BED: $mf -> $resolved" >&2
            missing=$((missing + 1))
        fi
    done < "$tmp"

    rm -f "$tmp"

    if [[ "$missing" -gt 0 ]]; then
        echo "[ERROR] $mf has $missing missing BED reference(s)" >&2
        exit 1
    fi
}

ldscore_prefix_from_manifest() {
    local mf="$1"
    local annot_id="$2"
    local row
    local chrom
    local out_prefix

    row=$(awk -v id="$annot_id" 'BEGIN{FS="\t"; OFS="\t"}
        $1 == id {
            gsub(/\r$/, "", $2)
            gsub(/\r$/, "", $4)
            print $2, $4
            exit
        }' "$mf")

    if [[ -z "$row" ]]; then
        echo "[ERROR] cannot find annotation_id=$annot_id in $mf" >&2
        exit 1
    fi

    IFS=$'\t' read -r chrom out_prefix <<< "$row"

    if [[ "$out_prefix" == *".${chrom}" ]]; then
        printf '%s\n' "${out_prefix%"${chrom}"}"
    elif [[ "$out_prefix" == *"${chrom}" ]]; then
        printf '%s\n' "${out_prefix%"${chrom}"}"
    else
        printf '%s\n' "$out_prefix"
    fi
}

build_ldcts_files() {
    local ldscore_manifest="$GWAS_BASE/manifests/ldscore_jobs.tsv"
    local control_manifest="$GWAS_BASE/manifests/control_ldscore_jobs.tsv"
    local ldcts_dir="$GWAS_BASE/ldcts"
    local tissues_with_muscle=(Adipose Colon Cortex Heart Liver Lung Muscle Ovary Sintest Spleen Stomach)
    local tissues_no_muscle=(Adipose Colon Cortex Heart Liver Lung Ovary Sintest Spleen Stomach)

    [[ -s "$ldscore_manifest" ]] || {
        echo "[WARN] cannot build LDCTS files; missing $ldscore_manifest"
        return 0
    }
    [[ -s "$control_manifest" ]] || {
        echo "[WARN] cannot build LDCTS files; missing $control_manifest"
        return 0
    }

    mkdir -p "$ldcts_dir"

    write_ldcts() {
        local class="$1"
        local control_id="$2"
        local out="$ldcts_dir/${class}.ldcts"
        shift 2
        local tissues=("$@")
        local control_prefix
        local annot_id
        local target_prefix
        local t

        control_prefix=$(ldscore_prefix_from_manifest "$control_manifest" "$control_id")
        : > "$out"

        for t in "${tissues[@]}"; do
            annot_id="${class}_ts_${t}"
            target_prefix=$(ldscore_prefix_from_manifest "$ldscore_manifest" "$annot_id")
            printf '%s\t%s,%s\n' "$t" "$target_prefix" "$control_prefix" >> "$out"
        done

        echo "[INFO] wrote LDCTS file: $out"
    }

    write_ldcts sfCRE sfCRE_ts_union "${tissues_no_muscle[@]}"
    write_ldcts sdCRE sdCRE_ts_union "${tissues_with_muscle[@]}"
    write_ldcts soCRE soCRE_ts_union "${tissues_with_muscle[@]}"
    write_ldcts ssCRE ssCRE_ts_union "${tissues_with_muscle[@]}"
    write_ldcts sf_proj sf_proj_ts_union "${tissues_no_muscle[@]}"
}

validate_cts_manifest() {
    local mf="$1"
    local missing=0
    local path

    [[ -s "$mf" ]] || return 0

    while IFS= read -r path; do
        [[ -n "$path" ]] || continue
        if [[ ! -s "$path" ]]; then
            echo "[ERROR] CTS manifest references missing LDCTS file: $mf -> $path" >&2
            missing=$((missing + 1))
        fi
    done < <(awk 'BEGIN{FS="\t"} NF >= 4 {
        path=$4
        sub(/\r$/, "", path)
        if (tolower(path) == "ldcts_file") next
        print path
    }' "$mf" | sort -u)

    if [[ "$missing" -gt 0 ]]; then
        echo "[ERROR] $mf has $missing missing LDCTS reference(s)" >&2
        exit 1
    fi
}

manifest_header() {
    case "$1" in
        make_annot_jobs.tsv|control_make_annot_jobs.tsv)
            printf 'annotation_id\tgroup\tbed\tchr\tout\tld_prefix\n' ;;
        ldscore_jobs.tsv|control_ldscore_jobs.tsv)
            printf 'annotation_id\tchrom\tannot\tout_prefix\n' ;;
        A1_h2_manifest.tsv|A1_h2_manifest_152.tsv|B1_h2_manifest.tsv|B1_h2_manifest_152.tsv)
            printf 'annotation_id\tcre_class\ttrait_id\tsumstats\tannot_prefix\toutdir\n' ;;
        A2_B2_cts_manifest.tsv|A2_B2_cts_manifest_152.tsv)
            printf 'ldcts_name\ttrait_id\tsumstats\tldcts_file\toutdir\n' ;;
        *) return 1 ;;
    esac
}

copy_manifest_with_header() {
    local src="$1"
    local dst="$2"
    local tmp="${dst}.tmp"
    local header
    local first

    sed "s#${OLD_GWAS_BASE}#${GWAS_BASE}#g" "$src" > "$tmp"
    header=$(manifest_header "$(basename "$src")") || {
        echo "[WARN] skipping ungoverned manifest schema: $src" >&2
        rm -f "$tmp"
        return 0
    }
    first=$(head -n 1 "$tmp" | cut -f1 | tr -d '\r')
    {
        printf '%s\n' "$header"
        case "$first" in
            annotation_id|ldcts_name) tail -n +2 "$tmp" ;;
            *) cat "$tmp" ;;
        esac
    } > "$dst"
    rm -f "$tmp"
}

for f in "$SCRIPT_DIR"/*; do
    [[ -f "$f" ]] || continue
    dst="$GWAS_BASE/scripts/$(basename "$f")"
    if [[ ! -e "$dst" || "$(readlink -f "$f")" != "$(readlink -f "$dst")" ]]; then
        cp -f "$f" "$dst"
    fi
done

if [[ "${RUN_MOVE_BED:-0}" == "1" ]]; then
    echo "[INFO] collecting final BED files into $GWAS_BASE/beds"
    GWAS_BASE="$GWAS_BASE" bash "$GWAS_BASE/scripts/00_move_bed.sh"
fi

if compgen -G "$GWAS_BASE/beds/*.bed" > /dev/null; then
    bed_count=$(find "$GWAS_BASE/beds" -maxdepth 1 -name '*.bed' | wc -l)
    echo "[INFO] BED count: $bed_count"
else
    echo "[WARN] no BED files found under $GWAS_BASE/beds"
fi

if [[ "${BUILD_CONTROL_BEDS:-1}" == "1" ]]; then
    echo "[INFO] building CTS control BEDs under $GWAS_BASE/control_beds"
    build_control_beds
fi

if [[ -d "$OLD_GWAS_BASE/manifests" ]]; then
    for mf in "$OLD_GWAS_BASE"/manifests/*.tsv; do
        [[ -f "$mf" ]] || continue
        out="$GWAS_BASE/manifests/$(basename "$mf")"
        copy_manifest_with_header "$mf" "$out"
        echo "[INFO] wrote manifest: $out"
    done
else
    echo "[WARN] old manifest directory not found: $OLD_GWAS_BASE/manifests"
fi

if [[ "${VALIDATE_MANIFEST_BEDS:-1}" == "1" ]]; then
    echo "[INFO] validating BED paths in make-annot manifests"
    validate_bed_manifest "$GWAS_BASE/manifests/make_annot_jobs.tsv"
    validate_bed_manifest "$GWAS_BASE/manifests/control_make_annot_jobs.tsv"
fi

if [[ "${BUILD_LDCTS:-1}" == "1" ]]; then
    echo "[INFO] building CTS LDCTS files under $GWAS_BASE/ldcts"
    build_ldcts_files
fi

if [[ "${VALIDATE_LDCTS:-1}" == "1" ]]; then
    echo "[INFO] validating LDCTS paths in CTS manifests"
    validate_cts_manifest "$GWAS_BASE/manifests/A2_B2_cts_manifest.tsv"
    validate_cts_manifest "$GWAS_BASE/manifests/A2_B2_cts_manifest_152.tsv"
fi

if [[ -s "$OLD_GWAS_BASE/all_152_traits.tsv" ]]; then
    if [[ -e "$GWAS_BASE/all_152_traits.tsv" ]]; then
        cmp -s "$OLD_GWAS_BASE/all_152_traits.tsv" "$GWAS_BASE/all_152_traits.tsv" || {
            echo "[ERROR] existing all_152_traits.tsv differs from historical source; refusing overwrite" >&2
            exit 1
        }
    else
        cp "$OLD_GWAS_BASE/all_152_traits.tsv" "$GWAS_BASE/all_152_traits.tsv"
    fi
else
    echo "[WARN] old all_152_traits.tsv not found; run expand_to_152.py after setup"
fi

if [[ -s "$OLD_GWAS_BASE/selected_41_traits.tsv" ]]; then
    if [[ -e "$GWAS_BASE/selected_41_traits.tsv" ]]; then
        cmp -s "$OLD_GWAS_BASE/selected_41_traits.tsv" "$GWAS_BASE/selected_41_traits.tsv" || {
            echo "[ERROR] existing selected_41_traits.tsv differs from historical source; refusing overwrite" >&2
            exit 1
        }
    else
        cp "$OLD_GWAS_BASE/selected_41_traits.tsv" "$GWAS_BASE/selected_41_traits.tsv"
    fi
fi

echo "[DONE] GWAS recalculation workspace is prepared."
echo "[NEXT] cd $GWAS_BASE && bash scripts/submit_E1E9_reuse_existing_152.sh"
