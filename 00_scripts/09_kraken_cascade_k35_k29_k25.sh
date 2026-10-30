#!/usr/bin/env bash

#SBATCH --job-name=sed_kraken_cascade
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --partition=smp
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=/home/plstenge/seda_DNA_Corsican_wreck_comparison/00_scripts/09_kraken_cascade_k35_k29_k25_%j.out
#SBATCH --error=/home/plstenge/seda_DNA_Corsican_wreck_comparison/00_scripts/09_kraken_cascade_k35_k29_k25_%j.err

# =============================================================================
# Title
#   Sequential Kraken2 classification and cumulative MPA profiling of sediment
#   samples sed6 and sed8.
#
# Analytical design
#   1. Classify all merged reads against the k35 database.
#   2. Classify only reads remaining unclassified after k35 against k29.
#   3. Classify only reads remaining unclassified after k29 against k25.
#
# Primary outputs
#   For each sample:
#       <sample>_k35.mpa
#       <sample>_k35_k29.mpa
#       <sample>_k35_k29_k25.mpa
#
#   Across samples:
#       combined_sed6_sed8_k35.tsv
#       combined_sed6_sed8_k35_k29.tsv
#       combined_sed6_sed8_k35_k29_k25.tsv
#
# Methodological considerations
#   - Merged reads are analysed as single-end sequences.
#   - Original input files are accessed in place and are never modified.
#   - MPA profiles contain integer clade counts, not relative abundances.
#   - Each read can contribute a classified assignment at only one stage.
#   - Cumulative profiles preserve earlier assignments without reclassification.
#   - Database taxonomies must be compatible for cumulative aggregation.
#
# Storage considerations
#   Intermediate unclassified FASTQ files are uncompressed.
#   These derived files may be substantially larger than the compressed inputs.
#   Temporary FASTQ files are deleted as soon as they are no longer required.
#
# Execution policy
#   Samples and stages are processed sequentially.
#   Each execution creates a separate run directory.
#   Failed runs are retained for diagnosis; automatic resumption is not provided.
# =============================================================================


# =============================================================================
# 1. Initialise the software environment
# =============================================================================

# Environment activation is performed before enabling nounset because some
# Conda activation scripts reference unset shell variables.

#set -eo pipefail

module load conda/4.12.0
source ~/.bashrc
conda activate metagenomics

#set -Eeuo pipefail
#export LC_ALL=C
#umask 0022

trap 'printf "[ERROR] Failure at line %s: %s\n" \
    "${LINENO}" "${BASH_COMMAND}" >&2' ERR


# =============================================================================
# 2. Define paths and analytical parameters
# =============================================================================

PROJECT="/home/plstenge/seda_DNA_Corsican_wreck_comparison"
ANALYSIS="${PROJECT}/09_kraken_cascade_k35_k29_k25"

INPUT_DIR="/storage/groups/gdec/shared_paleo/CapoSagro/fastp"
KRAKENTOOLS_DIR="${ANALYSIS}/tools/KrakenTools"

DB_K35="/storage/groups/gdec/shared/Kraken_database/k2_core_nt_20250609"
DB_K29="/storage/groups/gdec/shared/Kraken_database/core_nt_k29"
DB_K25="/storage/groups/gdec/shared/Kraken_database/core_nt_k25"

THREADS="${SLURM_CPUS_PER_TASK:-36}"
CONFIDENCE="0.2"

SAMPLES=("sed6" "sed8")

# A unique directory prevents accidental replacement of previous results.
RUN_ID="run_$(date +%Y%m%d_%H%M%S)_${SLURM_JOB_ID:-manual_$$}"
OUT="${ANALYSIS}/runs/${RUN_ID}"

REPORT_DIR="${OUT}/kraken"
STAGE_MPA_DIR="${OUT}/mpa/stage_specific"
CUMULATIVE_MPA_DIR="${OUT}/mpa/cumulative"
MATRIX_DIR="${OUT}/mpa/matrices"
QC_DIR="${OUT}/qc"
TMP_DIR="${OUT}/tmp_unclassified"

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

die() {
    printf '[ERROR] %s\n' "$*" >&2
    exit 1
}


# =============================================================================
# 3. Validate software, input files, and database components
# =============================================================================

for executable in kraken2 python3 awk; do
    command -v "${executable}" >/dev/null 2>&1 \
        || die "Required executable not found: ${executable}"
done

[[ -r "${KRAKENTOOLS_DIR}/kreport2mpa.py" ]] \
    || die "Missing KrakenTools script: ${KRAKENTOOLS_DIR}/kreport2mpa.py"

for sample in "${SAMPLES[@]}"; do
    input="${INPUT_DIR}/clean_${sample}_concat_dedup_fastp_merged.fastq.gz"

    [[ -r "${input}" && -s "${input}" ]] \
        || die "Input file is missing, unreadable, or empty: ${input}"
done

for database in "${DB_K35}" "${DB_K29}" "${DB_K25}"; do
    for component in hash.k2d opts.k2d taxo.k2d; do
        [[ -r "${database}/${component}" && -s "${database}/${component}" ]] \
            || die "Missing or unreadable database component: ${database}/${component}"
    done
done

mkdir -p "${REPORT_DIR}" \
         "${STAGE_MPA_DIR}" \
         "${CUMULATIVE_MPA_DIR}" \
         "${MATRIX_DIR}" \
         "${QC_DIR}" \
         "${TMP_DIR}"

log "Run directory: ${OUT}"
log "Samples: ${SAMPLES[*]}"
log "Threads: ${THREADS}"
log "Kraken2 confidence threshold: ${CONFIDENCE}"


# =============================================================================
# 4. Record provenance
# =============================================================================

{
    printf 'parameter\tvalue\n'
    printf 'run_id\t%s\n' "${RUN_ID}"
    printf 'date\t%s\n' "$(date -Is)"
    printf 'hostname\t%s\n' "$(hostname)"
    printf 'slurm_job_id\t%s\n' "${SLURM_JOB_ID:-not_available}"
    printf 'conda_environment\t%s\n' "${CONDA_DEFAULT_ENV:-unknown}"
    printf 'threads\t%s\n' "${THREADS}"
    printf 'confidence\t%s\n' "${CONFIDENCE}"
    printf 'input_directory\t%s\n' "${INPUT_DIR}"
    printf 'database_k35\t%s\n' "${DB_K35}"
    printf 'database_k29\t%s\n' "${DB_K29}"
    printf 'database_k25\t%s\n' "${DB_K25}"

    for sample in "${SAMPLES[@]}"; do
        printf 'input_%s\t%s\n' "${sample}" \
            "${INPUT_DIR}/clean_${sample}_concat_dedup_fastp_merged.fastq.gz"
    done
} > "${QC_DIR}/run_parameters.tsv"

kraken2 --version > "${QC_DIR}/kraken2_version.txt" 2>&1
python3 --version > "${QC_DIR}/python_version.txt" 2>&1
conda list > "${QC_DIR}/conda_packages.txt"

if command -v git >/dev/null 2>&1 \
    && git -C "${KRAKENTOOLS_DIR}" rev-parse HEAD \
        > "${QC_DIR}/krakentools_commit.txt" 2>/dev/null; then
    :
else
    printf 'Commit information unavailable.\n' \
        > "${QC_DIR}/krakentools_commit.txt"
fi

printf 'sample\tstage\tinput_reads\tclassified_reads\tunclassified_reads\n' \
    > "${QC_DIR}/classification_by_stage.tsv"


# =============================================================================
# 5. Define a single classification stage
# =============================================================================

run_stage() {
    local sample="$1"
    local stage="$2"
    local database="$3"
    local input="$4"
    local unclassified_output="$5"
    local expected_reads="$6"

    local stage_dir="${REPORT_DIR}/${stage}"
    local raw_output="${stage_dir}/${sample}_${stage}.kraken"
    local report="${stage_dir}/${sample}_${stage}.report"
    local stage_log="${stage_dir}/${sample}_${stage}.log"
    local mpa="${STAGE_MPA_DIR}/${sample}_${stage}.mpa"

    local total classified unclassified
    local -a options=()

    mkdir -p "${stage_dir}"

    # An empty residual population requires no further classification.
    # Empty stage-specific files represent a genuine zero contribution.
    if [[ "${expected_reads}" == "0" ]]; then
        log "${sample}: skipping ${stage}; no residual reads."

        : > "${raw_output}"
        : > "${report}"
        : > "${mpa}"

        printf 'Stage skipped: no residual reads.\n' > "${stage_log}"

        if [[ "${unclassified_output}" != "-" ]]; then
            : > "${unclassified_output}"
        fi

        printf '%s\t%s\t0\t0\t0\n' "${sample}" "${stage}" \
            >> "${QC_DIR}/classification_by_stage.tsv"

        STAGE_TOTAL=0
        STAGE_CLASSIFIED=0
        STAGE_UNCLASSIFIED=0
        return
    fi

    [[ -r "${input}" && -s "${input}" ]] \
        || die "Expected non-empty input for ${sample}/${stage}: ${input}"

    # Compression is explicitly specified only for the original gzipped input.
    if [[ "${input}" == *.gz ]]; then
        options+=(--gzip-compressed)
    fi

    # Residual FASTQ output is required at k35 and k29, but not at k25.
    if [[ "${unclassified_output}" != "-" ]]; then
        options+=(--unclassified-out "${unclassified_output}")
    fi

    log "${sample}: starting ${stage}."

    kraken2 \
        --db "${database}" \
        --threads "${THREADS}" \
        --confidence "${CONFIDENCE}" \
        --output "${raw_output}" \
        --report "${report}" \
        "${options[@]}" \
        "${input}" \
        > "${stage_log}" 2>&1

    # Every input sequence must generate exactly one C or U output record.
    # These counts also provide stage-to-stage consistency checks.
    read -r total classified unclassified < <(
        awk -F '\t' '
            $1 == "C" { c++; next }
            $1 == "U" { u++; next }
            { bad++ }
            END {
                if (bad > 0) {
                    print "Invalid Kraken output records." > "/dev/stderr"
                    exit 1
                }
                printf "%.0f %.0f %.0f\n", c + u, c + 0, u + 0
            }
        ' "${raw_output}"
    )

    [[ "${total}" -gt 0 ]] \
        || die "No classification records produced for ${sample}/${stage}."

    if [[ "${expected_reads}" != "unknown" ]]; then
        [[ "${total}" -eq "${expected_reads}" ]] \
            || die "Read-count mismatch for ${sample}/${stage}: expected ${expected_reads}, observed ${total}."
    fi

    if [[ "${unclassified_output}" != "-" ]]; then
        if [[ "${unclassified}" -gt 0 ]]; then
            [[ -s "${unclassified_output}" ]] \
                || die "Residual reads were reported, but residual FASTQ is empty."
        else
            # Ensure that the next stage receives an explicit empty file.
            : > "${unclassified_output}"
        fi
    fi

    [[ -s "${report}" ]] \
        || die "Kraken report is empty for ${sample}/${stage}."

    # Explicitly request clade read counts and standard taxonomic ranks.
    # The resulting two-column files have no header.
    python3 "${KRAKENTOOLS_DIR}/kreport2mpa.py" \
        -r "${report}" \
        -o "${mpa}" \
        --read-count \
        --no-intermediate-ranks \
        >> "${stage_log}" 2>&1

    printf '%s\t%s\t%s\t%s\t%s\n' \
        "${sample}" "${stage}" "${total}" "${classified}" "${unclassified}" \
        >> "${QC_DIR}/classification_by_stage.tsv"

    STAGE_TOTAL="${total}"
    STAGE_CLASSIFIED="${classified}"
    STAGE_UNCLASSIFIED="${unclassified}"

    log "${sample}: ${stage} completed; classified=${classified}; unclassified=${unclassified}."
}


# =============================================================================
# 6. Execute the cascade independently for each sample
# =============================================================================

for sample in "${SAMPLES[@]}"; do
    log "============================================================"
    log "Processing sample: ${sample}"

    original="${INPUT_DIR}/clean_${sample}_concat_dedup_fastp_merged.fastq.gz"

    residual_k35="${TMP_DIR}/${sample}_unclassified_after_k35.fastq"
    residual_k29="${TMP_DIR}/${sample}_unclassified_after_k29.fastq"

    # Stage 1: classify the complete original merged-read dataset.
    run_stage \
        "${sample}" "k35" "${DB_K35}" \
        "${original}" "${residual_k35}" "unknown"

    n_original="${STAGE_TOTAL}"
    n_c35="${STAGE_CLASSIFIED}"
    n_u35="${STAGE_UNCLASSIFIED}"

    # Stage 2: classify only reads left unclassified by k35.
    run_stage \
        "${sample}" "k29" "${DB_K29}" \
        "${residual_k35}" "${residual_k29}" "${n_u35}"

    n_c29="${STAGE_CLASSIFIED}"
    n_u29="${STAGE_UNCLASSIFIED}"

    # The k35 residual FASTQ is no longer needed once k29 has completed.
    rm -f -- "${residual_k35}"

    # Stage 3: classify only reads left unclassified by k29.
    # Final unclassified sequences are not written as a separate FASTQ.
    run_stage \
        "${sample}" "k25" "${DB_K25}" \
        "${residual_k29}" "-" "${n_u29}"

    n_c25="${STAGE_CLASSIFIED}"
    n_u25="${STAGE_UNCLASSIFIED}"

    rm -f -- "${residual_k29}"

    # Classified populations are disjoint by construction.
    # Their sum plus the final unclassified population must equal the input.
    [[ $((n_c35 + n_c29 + n_c25 + n_u25)) -eq "${n_original}" ]] \
        || die "Global read-count conservation failed for ${sample}."

    log "${sample}: read-count conservation verified."
done


# =============================================================================
# 7. Generate cumulative MPA profiles, cross-sample matrices, and QC summaries
# =============================================================================

export OUT

python3 <<'PYTHON'
import csv
import os
from collections import Counter
from pathlib import Path

root = Path(os.environ["OUT"])
samples = ("sed6", "sed8")

stage_dir = root / "mpa" / "stage_specific"
cumulative_dir = root / "mpa" / "cumulative"
matrix_dir = root / "mpa" / "matrices"
qc_dir = root / "qc"

profiles = {
    "k35": ("k35",),
    "k35_k29": ("k35", "k29"),
    "k35_k29_k25": ("k35", "k29", "k25"),
}


def load_mpa(path):
    """Read an MPA profile containing integer clade counts."""
    counts = Counter()

    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")

            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")

            if len(fields) != 2:
                raise ValueError(
                    f"Invalid MPA record at {path}:{line_number}"
                )

            lineage, raw_count = fields
            count = int(raw_count)

            if count < 0:
                raise ValueError(
                    f"Negative count at {path}:{line_number}"
                )

            counts[lineage] += count

    return counts


# Load each stage once. Empty files correctly represent zero contribution.
stage_counts = {
    sample: {
        stage: load_mpa(stage_dir / f"{sample}_{stage}.mpa")
        for stage in ("k35", "k29", "k25")
    }
    for sample in samples
}

cumulative_counts = {}

for profile_name, stages in profiles.items():
    cumulative_counts[profile_name] = {}

    for sample in samples:
        combined = Counter()

        for stage in stages:
            combined.update(stage_counts[sample][stage])

        cumulative_counts[profile_name][sample] = combined

        output = cumulative_dir / f"{sample}_{profile_name}.mpa"

        with output.open("w") as handle:
            for lineage in sorted(combined):
                handle.write(f"{lineage}\t{combined[lineage]}\n")

    # Use the union of observed lineages and fill absent taxa with zero.
    # Explicit sample names facilitate subsequent import into R.
    all_lineages = sorted(
        set().union(
            *(set(cumulative_counts[profile_name][s]) for s in samples)
        )
    )

    output = matrix_dir / f"combined_sed6_sed8_{profile_name}.tsv"

    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["#Classification", *samples])

        for lineage in all_lineages:
            writer.writerow([
                lineage,
                *[
                    cumulative_counts[profile_name][sample].get(lineage, 0)
                    for sample in samples
                ],
            ])


# Derive cumulative classification statistics from per-stage read counts.
statistics = {}

with (qc_dir / "classification_by_stage.tsv").open() as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        statistics[(row["sample"], row["stage"])] = {
            key: int(row[key])
            for key in (
                "input_reads",
                "classified_reads",
                "unclassified_reads",
            )
        }

with (qc_dir / "classification_cumulative.tsv").open(
    "w", newline=""
) as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow([
        "sample",
        "profile",
        "total_input_reads",
        "cumulative_classified_reads",
        "remaining_unclassified_reads",
        "classified_percent",
    ])

    for sample in samples:
        total = statistics[(sample, "k35")]["input_reads"]

        for profile_name, stages in profiles.items():
            classified = sum(
                statistics[(sample, stage)]["classified_reads"]
                for stage in stages
            )

            remaining = statistics[
                (sample, stages[-1])
            ]["unclassified_reads"]

            if classified + remaining != total:
                raise RuntimeError(
                    f"Read-count conservation failed: "
                    f"{sample}/{profile_name}"
                )

            writer.writerow([
                sample,
                profile_name,
                total,
                classified,
                remaining,
                f"{100 * classified / total:.6f}",
            ])
PYTHON


# =============================================================================
# 8. Finalise the analysis
# =============================================================================

# Remove the temporary directory only if it is empty.
rmdir "${TMP_DIR}"

# This marker is created only after all preceding steps have succeeded.
date -Is > "${OUT}/SUCCESS"

log "============================================================"
log "All analyses completed successfully."
log "Run directory: ${OUT}"
log "Cumulative sample profiles: ${CUMULATIVE_MPA_DIR}"
log "Cross-sample matrices: ${MATRIX_DIR}"
log "Quality-control summaries: ${QC_DIR}"
log "Original input FASTQ files were neither copied nor modified."
log "============================================================"
