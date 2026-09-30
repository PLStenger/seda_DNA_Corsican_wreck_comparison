#!/usr/bin/env python3

import csv
import gzip
import sys
from pathlib import Path

RUNS = {
    "RUN1": Path("/storage/groups/gdec/shared_paleo/Illumina/01_raw_data"),
    "RUN2": Path("/storage/groups/gdec/shared_paleo/E1531_final/run1_20250320_AV241601_E1531_Ps5Lane1_Ps6Lane2"),
    "RUN3": Path("/storage/groups/gdec/shared_paleo/E1531_final/run2_20250414_AV241601_E1531_Ps5_Ps6_11022026_CORRECTED"),
    "RUN4": Path("/storage/groups/gdec/shared_paleo/E1531_final/run3_20251008_AV241601_E1531_Ps7_Ps8"),
    "RUN5": Path("/storage/groups/gdec/shared_paleo/E1531_final/run4_20251104_AV241601_E1531_Ps7_Ps8_04112025"),
    "RUN6": Path("/storage/groups/gdec/shared_paleo/E1672"),
    "RUN7": Path("/storage/groups/gdec/shared_paleo/E1672/Ps10_Ps11"),
}

SAMPLES = (
    "1120_sed6_rep3",
    "1129_sed6_rep1",
    "1130_sed6_rep2",
    "1121_sed8_rep1",
    "1122_sed8_rep2",
    "1131_sed8_rep3",
)

# Échantillons explicitement présents dans la liste fournie.
PRESENT = {
    "RUN1": set(SAMPLES),
    "RUN2": set(SAMPLES),
    "RUN3": set(SAMPLES),
    "RUN4": {
        "1129_sed6_rep1", "1130_sed6_rep2",
        "1121_sed8_rep1", "1122_sed8_rep2",
    },
    "RUN5": {
        "1129_sed6_rep1", "1130_sed6_rep2",
        "1121_sed8_rep1", "1122_sed8_rep2",
    },
    "RUN6": {
        "1129_sed6_rep1", "1130_sed6_rep2",
        "1121_sed8_rep1", "1122_sed8_rep2",
    },
    "RUN7": {
        "1129_sed6_rep1", "1130_sed6_rep2",
        "1121_sed8_rep1", "1122_sed8_rep2",
    },
}

OUTPUT = Path("nombre_sequences_par_run.csv")


def fastq_path(run, sample, read):
    root = RUNS[run]

    if run == "RUN1":
        return root / f"{sample}_{read}.fastq.gz"

    if run in {"RUN6", "RUN7"}:
        file_sample = sample.replace("_rep", "-rep")
        return root / file_sample / f"{file_sample}_{read}.fastq.gz"

    return root / sample / f"{sample}_{read}.fastq.gz"


def count_sequences(path):
    """Compte les lignes décompressées, sans charger tout le FASTQ en mémoire."""
    if not path.is_file():
        return None, "FICHIER_MANQUANT"

    lines = 0
    try:
        with gzip.open(path, "rb") as handle:
            while chunk := handle.read(8 * 1024 * 1024):
                lines += chunk.count(b"\n")
    except (OSError, EOFError) as exc:
        return None, f"ERREUR_LECTURE: {exc}"

    if lines % 4 != 0:
        return None, f"FASTQ_INCOMPLET: {lines} lignes"

    return lines // 4, "OK"


def main():
    problems = 0

    with OUTPUT.open("w", newline="", encoding="utf-8") as output:
        writer = csv.writer(output)
        writer.writerow([
            "run", "echantillon", "sediment", "replicat",
            "sequences_R1", "sequences_R2", "paires",
            "statut", "fichier_R1", "fichier_R2",
        ])

        for run in RUNS:
            for sample in SAMPLES:
                sediment, replicat = sample.split("_")[1:]
                r1_path = fastq_path(run, sample, "R1")
                r2_path = fastq_path(run, sample, "R2")

                if sample not in PRESENT[run]:
                    writer.writerow([
                        run, sample, sediment, replicat,
                        "", "", "", "ABSENT", "", "",
                    ])
                    continue

                n1, status1 = count_sequences(r1_path)
                n2, status2 = count_sequences(r2_path)

                if status1 != "OK" or status2 != "OK":
                    statut = f"R1={status1};R2={status2}"
                    problems += 1
                elif n1 != n2:
                    statut = "R1_R2_DIFFERENTS"
                    problems += 1
                else:
                    statut = "OK"

                writer.writerow([
                    run, sample, sediment, replicat,
                    n1 if n1 is not None else "",
                    n2 if n2 is not None else "",
                    n1 if statut == "OK" else "",
                    statut, str(r1_path), str(r2_path),
                ])

    print(f"Tableau créé : {OUTPUT.resolve()}")
    print(f"Anomalies détectées : {problems}")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
