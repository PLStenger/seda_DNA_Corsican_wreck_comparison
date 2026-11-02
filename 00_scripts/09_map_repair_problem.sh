cat > repair_mpa_9307771.sh <<'BASH'
#!/usr/bin/env bash
set -Eeuo pipefail

trap 'printf "[ERROR] Ligne %s : %s\n" \
    "${LINENO}" "${BASH_COMMAND}" >&2' ERR

export OUT="/home/plstenge/seda_DNA_Corsican_wreck_comparison/09_kraken_cascade_k35_k29_k25/runs/run_20261030_162223_9307771"

CONVERTER="/home/plstenge/coprolites_comparison/08_krakentools/KrakenTools/kreport2mpa.py"

[[ -d "${OUT}" ]] || {
    echo "[ERROR] Run introuvable : ${OUT}" >&2
    exit 1
}

# Le marqueur existant ne prouve pas la réussite du post-traitement.
if [[ -f "${OUT}/SUCCESS" ]]; then
    mv -- "${OUT}/SUCCESS" \
        "${OUT}/SUCCESS.before_mpa_repair.$(date +%Y%m%d_%H%M%S)"
fi

[[ -r "${CONVERTER}" ]] || {
    echo "[ERROR] Convertisseur introuvable : ${CONVERTER}" >&2
    exit 1
}

# Vérifier d'abord tous les rapports avant de lancer les conversions.
for sample in sed6 sed8; do
    for stage in k35 k29 k25; do
        report="${OUT}/kraken/${stage}/${sample}_${stage}.report"

        [[ -r "${report}" ]] || {
            echo "[ERROR] Rapport absent ou illisible : ${report}" >&2
            exit 1
        }
    done
done

mkdir -p \
    "${OUT}/mpa/stage_specific" \
    "${OUT}/mpa/cumulative" \
    "${OUT}/mpa/matrices" \
    "${OUT}/qc"

for sample in sed6 sed8; do
    for stage in k35 k29 k25; do
        report="${OUT}/kraken/${stage}/${sample}_${stage}.report"
        mpa="${OUT}/mpa/stage_specific/${sample}_${stage}.mpa"
        log="${OUT}/kraken/${stage}/${sample}_${stage}.mpa_repair.log"

        echo "[INFO] Conversion ${sample}/${stage}"

        if [[ -s "${report}" ]]; then
            python3 "${CONVERTER}" \
                -r "${report}" \
                -o "${mpa}.tmp" \
                --read_count \
                --no-intermediate-ranks \
                > "${log}" 2>&1

            mv -- "${mpa}.tmp" "${mpa}"
        else
            # Le contrôle Python vérifiera que cette étape a bien zéro reads.
            : > "${mpa}"
            printf 'Rapport vide : vérification requise dans le QC.\n' > "${log}"
        fi

        [[ -f "${mpa}" ]] || {
            echo "[ERROR] MPA non créé : ${mpa}" >&2
            exit 1
        }
    done
done

python3 <<'PYTHON'
import csv
import os
from collections import Counter
from pathlib import Path

root = Path(os.environ["OUT"])
samples = ("sed6", "sed8")
stages = ("k35", "k29", "k25")

profiles = (
    ("k35", ("k35",)),
    ("k35_k29", ("k35", "k29")),
    ("k35_k29_k25", ("k35", "k29", "k25")),
)

stage_dir = root / "mpa" / "stage_specific"
cumulative_dir = root / "mpa" / "cumulative"
matrix_dir = root / "mpa" / "matrices"
qc_dir = root / "qc"

statistics = {}
qc_file = qc_dir / "classification_by_stage.tsv"

with qc_file.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    required = {
        "sample", "stage", "input_reads",
        "classified_reads", "unclassified_reads",
    }

    if not required.issubset(set(reader.fieldnames or [])):
        raise ValueError("Colonnes QC manquantes dans {}".format(qc_file))

    for row in reader:
        key = (row["sample"], row["stage"])

        if key in statistics:
            raise ValueError("Ligne QC dupliquée : {}".format(key))

        values = {
            field: int(row[field])
            for field in (
                "input_reads", "classified_reads", "unclassified_reads"
            )
        }

        if any(value < 0 for value in values.values()):
            raise ValueError("Valeur QC négative : {}".format(key))

        if (
            values["classified_reads"] + values["unclassified_reads"]
            != values["input_reads"]
        ):
            raise ValueError("Conservation des reads invalide : {}".format(key))

        statistics[key] = values

for sample in samples:
    for stage in stages:
        if (sample, stage) not in statistics:
            raise ValueError("QC manquant : {}/{}".format(sample, stage))

    for previous, current in zip(stages, stages[1:]):
        if (
            statistics[(sample, current)]["input_reads"]
            != statistics[(sample, previous)]["unclassified_reads"]
        ):
            raise ValueError(
                "Cascade incohérente : {}/{} -> {}".format(
                    sample, previous, current
                )
            )

def load_mpa(path):
    counts = Counter()

    with path.open() as handle:
        for number, line in enumerate(handle, 1):
            line = line.rstrip("\r\n")

            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) != 2 or not fields[0]:
                raise ValueError("MPA invalide : {}:{}".format(path, number))

            lineage, raw_count = fields
            count = int(raw_count)

            if count < 0:
                raise ValueError("Compte négatif : {}:{}".format(path, number))

            counts[lineage] += count

    return counts

stage_counts = {}

for sample in samples:
    stage_counts[sample] = {}

    for stage in stages:
        report = root / "kraken" / stage / (
            "{}_{}.report".format(sample, stage)
        )
        values = statistics[(sample, stage)]

        if report.stat().st_size == 0 and values["input_reads"] != 0:
            raise ValueError(
                "Rapport vide malgré des reads : {}/{}".format(sample, stage)
            )

        path = stage_dir / "{}_{}.mpa".format(sample, stage)
        counts = load_mpa(path)

        if values["classified_reads"] > 0 and not counts:
            raise ValueError(
                "MPA sans taxons malgré des reads classifiés : {}".format(path)
            )

        stage_counts[sample][stage] = counts
        print("[OK] {} : {} lignées".format(path.name, len(counts)))

qc_rows = []

for profile_name, selected_stages in profiles:
    combined_by_sample = {}

    for sample in samples:
        combined = Counter()
        for stage in selected_stages:
            combined.update(stage_counts[sample][stage])

        combined_by_sample[sample] = combined

        output = cumulative_dir / "{}_{}.mpa".format(sample, profile_name)
        with output.open("w") as handle:
            for lineage in sorted(combined):
                handle.write("{}\t{}\n".format(lineage, combined[lineage]))

        total = statistics[(sample, "k35")]["input_reads"]
        classified = sum(
            statistics[(sample, stage)]["classified_reads"]
            for stage in selected_stages
        )
        remaining = statistics[
            (sample, selected_stages[-1])
        ]["unclassified_reads"]

        if classified + remaining != total:
            raise ValueError(
                "Conservation cumulative invalide : {}/{}".format(
                    sample, profile_name
                )
            )

        percent = (
            "{:.6f}".format(100 * classified / total)
            if total else "NA"
        )
        qc_rows.append([
            sample, profile_name, total, classified, remaining, percent
        ])

    lineages = sorted(
        set().union(*(set(combined_by_sample[s]) for s in samples))
    )

    output = matrix_dir / "combined_sed6_sed8_{}.tsv".format(profile_name)
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["#Classification", *samples])
        for lineage in lineages:
            writer.writerow([
                lineage,
                *[combined_by_sample[s].get(lineage, 0) for s in samples],
            ])

    print("[OK] {} : {} lignées".format(output.name, len(lineages)))

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
    writer.writerows(qc_rows)
PYTHON

date -Is > "${OUT}/SUCCESS"

echo "[OK] MPA par étape, profils cumulatifs, matrices et QC reconstruits."

find "${OUT}/mpa" -type f -name '*.tmp' -prune -o \
    -type f \( -name '*.mpa' -o -name '*.tsv' \) -print

cat "${OUT}/qc/classification_cumulative.tsv"
BASH

bash repair_mpa_9307771.sh
