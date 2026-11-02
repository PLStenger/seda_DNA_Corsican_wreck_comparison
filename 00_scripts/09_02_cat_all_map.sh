python3 <<'PYTHON'
import csv
from collections import Counter
from pathlib import Path

input_dir = Path.cwd()
output = input_dir.parent / "matrices" / "combined_sed6_sed8_k35_k29_k25.tsv"

samples = ("sed6", "sed8")
stages = ("k35", "k29", "k25")

counts = {sample: Counter() for sample in samples}

# Vérifier la présence des six fichiers avant de commencer.
for sample in samples:
    for stage in stages:
        path = input_dir / "{}_{}.mpa".format(sample, stage)
        if not path.is_file():
            raise FileNotFoundError("Fichier manquant : {}".format(path))

# Additionner les trois étapes séparément pour chaque échantillon.
for sample in samples:
    for stage in stages:
        path = input_dir / "{}_{}.mpa".format(sample, stage)

        with path.open() as handle:
            for line_number, line in enumerate(handle, start=1):
                line = line.rstrip("\r\n")

                if not line or line.startswith("#"):
                    continue

                fields = line.split("\t")
                if len(fields) != 2 or not fields[0]:
                    raise ValueError(
                        "Format invalide : {}:{}".format(path, line_number)
                    )

                lineage, raw_count = fields
                count = int(raw_count)

                if count < 0:
                    raise ValueError(
                        "Compte négatif : {}:{}".format(path, line_number)
                    )

                counts[sample][lineage] += count

# Conserver toutes les lignées ; attribuer zéro si absentes d'un échantillon.
lineages = sorted(set(counts["sed6"]) | set(counts["sed8"]))

output.parent.mkdir(parents=True, exist_ok=True)

with output.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["#Classification", *samples])

    for lineage in lineages:
        writer.writerow([
            lineage,
            counts["sed6"].get(lineage, 0),
            counts["sed8"].get(lineage, 0),
        ])

print("Fichier créé : {}".format(output))
print("Nombre de lignées : {}".format(len(lineages)))
PYTHON
