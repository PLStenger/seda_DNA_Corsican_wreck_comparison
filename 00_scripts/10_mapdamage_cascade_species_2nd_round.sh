#!/usr/bin/env bash
#SBATCH --job-name=sed_mapdamage_cascade
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --partition=smp
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=/home/plstenge/seda_DNA_Corsican_wreck_comparison/00_scripts/10_mapdamage_cascade_%j.out
#SBATCH --error=/home/plstenge/seda_DNA_Corsican_wreck_comparison/00_scripts/10_mapdamage_cascade_%j.err

# Usage: sbatch 10_mapdamage_cascade_species.sh /chemin/09_kraken_cascade_k35_k29_k25/runs/run_...
# Aucun lancement de Kraken2. Utilise exclusivement une cascade terminee.
# k29 et k25 = contributions residuelles, PAS classifications independantes.
# Les references et leurs index ne sont jamais modifies : index locaux au run.
# Pas de deduplication supplementaire, pas de filtre MAPQ par defaut.
# Environnement unique : mapdamage_py39, comme dans le script 10 fourni.
#set -eo pipefail
module load conda/4.12.0
source ~/.bashrc
conda activate mapdamage_py39
#set -Eeuo pipefail
export LC_ALL=C
trap 'printf "ERREUR ligne %s : %s\n" "$LINENO" "$BASH_COMMAND" >&2' ERR

PROJECT=/home/plstenge/seda_DNA_Corsican_wreck_comparison
INPUT_DIR=/storage/groups/gdec/shared_paleo/CapoSagro/fastp
KRAKEN_RUN="${1:-/home/plstenge/seda_DNA_Corsican_wreck_comparison/09_kraken_cascade_k35_k29_k25/runs/run_20261030_162223_9307771}"
THREADS="${SLURM_CPUS_PER_TASK:-36}"
SORT_THREADS=4
SORT_MEM=2G
# Les taux tous MAPQ et MAPQ >= seuil sont rapportes separement.
# mapDamage utilise les alignements primaires mappes MAPQ >= seuil.
MIN_MAPQ="${MIN_MAPQ:-0}"
KEEP_EXTRACTED="${KEEP_EXTRACTED:-0}"
SAMPLES=(sed6 sed8)
PROFILES=(k35 k29 k25 k35_k29_k25)
RUN_ID="run_$(date +%Y%m%d_%H%M%S)_${SLURM_JOB_ID:-manual_$$}"
OUT="${PROJECT}/10_mapdamage_cascade/runs/${RUN_ID}"
[[ -d "$KRAKEN_RUN/kraken" ]] || {
    printf 'ERREUR : répertoire Kraken absent : %s\n' \
        "$KRAKEN_RUN/kraken" >&2
    exit 1
}

printf 'Run Kraken utilisé : %s\n' "$KRAKEN_RUN"
[[ "$MIN_MAPQ" =~ ^[0-9]+$ && "$MIN_MAPQ" -le 255 ]] || exit 1
[[ "$KEEP_EXTRACTED" == 0 || "$KEEP_EXTRACTED" == 1 ]] || exit 1
for cmd in python3 bwa samtools mapDamage; do command -v "$cmd" >/dev/null; done
mkdir -p "$OUT"/{qc,references,extracted,results}
log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
exec > >(tee -a "$OUT/pipeline.log") 2>&1
log "Sortie : $OUT"

# References cibles : uniquement les 22 taxons demandes.
cat > "$OUT/qc/targets.tsv" <<'TARGETS'
reference_id\ttaxon\ttaxid\treference_path
Homo_sapiens\tHomo_sapiens\t9606\t/home/plstenge/genomes/Homo_sapiens/GCF_000001405.40_GRCh38.p14_genomic.fna
Gobiusculus_flavescens\tGobiusculus_flavescens\t257540\t/home/plstenge/genomes/Gobiusculus_flavescens/Gobiusculus_flavescens_fGobFla1.fasta
Mus_musculus\tMus_musculus\t10090\t/home/plstenge/genomes/Mus_musculus/Mus_musculus.GRCm39.dna.toplevel.fa
Rattus_norvegicus\tRattus_norvegicus\t10116\t/home/plstenge/genomes/Rattus_norvegicus/GCF_036323735.1_GRCr8_genomic.fna
Quercus\tQuercus\t3511\t/home/plstenge/genomes/Quercus_robur/Qrob_PM1N.fa
Calluna_vulgaris\tCalluna_vulgaris\t13385\t/home/plstenge/genomes/Calluna_vulgaris/GCA_964145215.1_ddCalVulg4.hap2.1_genomic.fna
Cannabis_sativa\tCannabis_sativa\t3483\t/home/plstenge/genomes/Cannabis_sativa/GCF_029168945.1_ASM2916894v1_genomic.fna
Canis_lupus_familiaris\tCanis_lupus_familiaris\t9615\t/home/plstenge/genomes/Canis_lupus_familiaris/GCF_011100685.1_UU_Cfam_GSD_1.0_genomic.fna
Canis_lupus\tCanis_lupus\t9612\t/home/plstenge/genomes/Canis_lupus/GCA_905319855.2_mCanLor1.2_genomic.fna
Alnus_glutinosa\tAlnus_glutinosa\t3517\t/home/plstenge/genomes/Alnus_glutinosa/Alnus_glutinosa_genome_assembly_dhAlnGlut1.fa
Bertholletia_excelsa\tBertholletia_excelsa\t3645\t/home/plstenge/genomes/Bertholletia_excelsa/GCA_039639745.1_Bnut_1.0_genomic.fna
Brassica_oleracea\tBrassica_oleracea\t3712\t/home/plstenge/genomes/Brassica_oleracea/GCF_000695525.1_BOL_genomic.fna
Eragrostis_tef\tEragrostis_tef\t110835\t/home/plstenge/genomes/Eragrostis_tef/GCA_000970635.1_ASM97063v1_genomic.fna
Humulus_lupulus\tHumulus_lupulus\t3486\t/home/plstenge/genomes/Humulus_lupulus/GCF_963169125.1_drHumLupu1.1_genomic.fna
Malus_domestica\tMalus_domestica\t3750\t/home/plstenge/genomes/Malus_domestica/GCF_042453785.1_GDT2T_hap1_genomic.fna
Medicago_truncatula\tMedicago_truncatula\t3880\t/home/plstenge/genomes/Medicago_truncatula/GCF_003473485.1_MtrunA17r5.0-ANR_genomic.fna
Olea_europaea\tOlea_europaea\t4146\t/home/plstenge/genomes/Olea_europaea/CACTIH01.fasta
Oryza_sativa\tOryza_sativa\t4530\t/home/plstenge/genomes/Oryza_sativa/GCF_034140825.1_ASM3414082v1_genomic.fna
Populus_trichocarpa\tPopulus_trichocarpa\t3694\t/home/plstenge/genomes/Populus_trichocarpa/GCF_000002775.5_P.trichocarpa_v4.1_genomic.fna
Prunus_dulcis\tPrunus_dulcis\t3755\t/home/plstenge/genomes/Prunus_dulcis/CABIKO01.fasta
Sus_scrofa\tSus_scrofa\t9823\t/home/plstenge/genomes/Sus_scrofa/GCF_000003025.6_Sscrofa11.1_genomic.fna
Orcinus_orca\tOrcinus_orca\t9733\t/home/plstenge/genomes/Orcinus_orca/GCF_937001465.1_mOrcOrc1.1_genomic.fna
TARGETS
# Convertit les separateurs litteraux du manifeste en tabulations.
sed -i 's/\\t/\t/g' "$OUT/qc/targets.tsv"
while IFS=$'\t' read -r rid taxon taxid ref; do
    [[ "$rid" == reference_id ]] && continue
    [[ -r "$ref" && -s "$ref" ]] || { echo "Reference absente : $ref"; exit 1; }
done < "$OUT/qc/targets.tsv"
for sample in "${SAMPLES[@]}"; do
    [[ -s "$INPUT_DIR/clean_${sample}_concat_dedup_fastp_merged.fastq.gz" ]] || exit 1
    for stage in k35 k29 k25; do
        for ext in kraken report; do
            [[ -r "$KRAKEN_RUN/kraken/$stage/${sample}_${stage}.$ext" ]] || {
                echo "Entree absente : $sample $stage $ext"; exit 1;
            }
        done
    done
done
printf 'parameter\tvalue\nkraken_run\t%s\ninput_dir\t%s\nmin_mapq\t%s\nthreads\t%s\nbwa_aln\t-n 0.08 -l 24 -k 2 -q 20\n' \
    "$KRAKEN_RUN" "$INPUT_DIR" "$MIN_MAPQ" "$THREADS" > "$OUT/qc/parameters.tsv"
conda list > "$OUT/qc/conda_packages.txt"
cp -- "$0" "$OUT/qc/script.sh"

# Liste derivee du manifeste pour ne conserver aucune ancienne cible.
mapfile -t TARGET_TAXA < <(awk -F '\t' 'NR>1 {print $2}' "$OUT/qc/targets.tsv")

# Index Orcinus_orca cree au debut, avant extraction et alignement.
# Index locaux au run : le repertoire du genome source reste inchange.
ORCA_REF=/home/plstenge/genomes/Orcinus_orca/GCF_937001465.1_mOrcOrc1.1_genomic.fna
ORCA_DIR="$OUT/references/Orcinus_orca"
mkdir -p "$ORCA_DIR"
ln -s "$ORCA_REF" "$ORCA_DIR/reference.fa"
log "Creation initiale de l'index BWA : Orcinus_orca"
bwa index "$ORCA_DIR/reference.fa" > "$ORCA_DIR/bwa_index.log" 2>&1 || exit 1
samtools faidx "$ORCA_DIR/reference.fa" || exit 1

# Extraction synchronisee en une lecture du FASTQ original par echantillon.
# Pour chaque read : consulter k35, puis k29 seulement si U, puis k25 si U.
# Les descendants sont definis par l'indentation du report DU niveau concerne.
# Controle strict des identifiants, de la structure FASTQ et des fins de fichiers.
cat > "$OUT/qc/extract_cascade.py" <<'PY'
import csv, gzip, sys
from pathlib import Path
from contextlib import ExitStack
run, fastp, output, sample = map(str, sys.argv[1:])
root = Path(output)
stages = ('k35', 'k29', 'k25')
with (root/'qc/targets.tsv').open() as h:
    targets = {r['taxon']: r['taxid'] for r in csv.DictReader(h, delimiter='\t')}

def descendants(path):
    result = {t: {i} for t, i in targets.items()}
    stack = []
    with path.open() as h:
        for line in h:
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 6:
                raise ValueError(f'Report non standard : {path}')
            taxid, name = fields[4].strip(), fields[5]
            depth = len(name) - len(name.lstrip(' '))
            while stack and stack[-1][0] >= depth:
                stack.pop()
            ancestors = {x[1] for x in stack}
            for taxon, parent in targets.items():
                if taxid == parent or parent in ancestors:
                    result[taxon].add(taxid)
            stack.append((depth, taxid))
    return result

allowed = {s: descendants(Path(run)/'kraken'/s/f'{sample}_{s}.report') for s in stages}
counts = {(t,s): 0 for t in targets for s in stages}
seen = dict.fromkeys(stages, 0)
classified = dict.fromkeys(stages, 0)
n_input = 0
folder = root/'extracted'/sample
folder.mkdir(parents=True)
with ExitStack() as ctx:
    kh = {s: ctx.enter_context((Path(run)/'kraken'/s/f'{sample}_{s}.kraken').open()) for s in stages}
    fh = {(t,s): ctx.enter_context((folder/f'{t}_{s}.fastq').open('w')) for t in targets for s in stages}
    fq = ctx.enter_context(gzip.open(Path(fastp)/f'clean_{sample}_concat_dedup_fastp_merged.fastq.gz', 'rt'))
    while True:
        header = fq.readline()
        if not header:
            break
        seq, plus, qual = fq.readline(), fq.readline(), fq.readline()
        if not header.startswith('@') or not plus.startswith('+') or not qual or len(seq.rstrip('\r\n')) != len(qual.rstrip('\r\n')):
            raise ValueError(f'FASTQ invalide au read {n_input+1}')
        name = header[1:].split()[0]
        n_input += 1
        for stage in stages:
            fields = kh[stage].readline().rstrip('\n').split('\t')
            if len(fields) < 5 or fields[0] not in ('C','U') or fields[1] != name:
                raise ValueError(f'Desynchronisation {sample}/{stage}, read {n_input}, ID {name}')
            seen[stage] += 1
            if fields[0] == 'C':
                classified[stage] += 1
                for taxon in targets:
                    if fields[2] in allowed[stage][taxon]:
                        fh[taxon,stage].write(header+seq+plus+qual)
                        counts[taxon,stage] += 1
                break
    for stage in stages:
        if kh[stage].readline():
            raise ValueError(f'Records Kraken surnumeraires : {sample}/{stage}')
if seen['k29'] != seen['k35']-classified['k35'] or seen['k25'] != seen['k29']-classified['k29']:
    raise ValueError('Incoherence des effectifs de cascade')
with (root/'qc'/f'{sample}_extraction.tsv').open('w') as h:
    h.write('sample\ttaxon\tprofile\textracted_reads\n')
    for taxon in targets:
        for stage in stages:
            h.write(f'{sample}\t{taxon}\t{stage}\t{counts[taxon,stage]}\n')
        h.write(f'{sample}\t{taxon}\tk35_k29_k25\t{sum(counts[taxon,s] for s in stages)}\n')
with (root/'qc'/f'{sample}_cascade_counts.tsv').open('w') as h:
    h.write('stage\tinput_reads\tclassified_reads\n')
    for stage in stages:
        h.write(f'{stage}\t{seen[stage]}\t{classified[stage]}\n')
print(f'{sample}: {n_input} reads verifies et extraction terminee', flush=True)
PY

for sample in "${SAMPLES[@]}"; do
    python3 "$OUT/qc/extract_cascade.py" "$KRAKEN_RUN" "$INPUT_DIR" "$OUT" "$sample"
    for taxon in "${TARGET_TAXA[@]}"; do
        cat "$OUT/extracted/$sample/${taxon}_k35.fastq" \
            "$OUT/extracted/$sample/${taxon}_k29.fastq" \
            "$OUT/extracted/$sample/${taxon}_k25.fastq" \
            > "$OUT/extracted/$sample/${taxon}_k35_k29_k25.fastq"
    done
done

SUMMARY="$OUT/mapping_summary.tsv"
printf 'sample\ttaxon\ttaxid\treference_id\treference_path\tprofile\textracted_reads\tprimary_records\tmapped_reads\tmapping_rate_percent\tmapped_reads_mapq_ge_threshold\tmapping_rate_mapq_percent\tmin_mapq\tmapdamage_status\tbam_path\tmapdamage_dir\n' > "$SUMMARY"
failures=0
while IFS=$'\t' read -r rid taxon taxid ref; do
    [[ "$rid" == reference_id ]] && continue
    # Un seul index par genome. Les index existants sont reutilises via liens.
    # Les index manquants sont crees localement, jamais a cote de la source.
    refdir="$OUT/references/$rid"
    mkdir -p "$refdir"
    localref="$refdir/reference.fa"
    if [[ "$rid" != Orcinus_orca ]]; then
ln -s "$ref" "$localref"
    complete=1
    for ext in amb ann bwt pac sa; do [[ -s "$ref.$ext" ]] || complete=0; done
    if [[ "$complete" == 1 ]]; then
        for ext in amb ann bwt pac sa; do ln -s "$ref.$ext" "$localref.$ext"; done
    else
        log "Index BWA local : $rid"
        bwa index "$localref" > "$refdir/bwa_index.log" 2>&1
    fi
    if [[ -s "$ref.fai" ]]; then ln -s "$ref.fai" "$localref.fai"
    else samtools faidx "$localref"; fi
fi

for sample in "${SAMPLES[@]}"; do
        for profile in "${PROFILES[@]}"; do
            d="$OUT/results/$rid/$sample/$profile"
            mkdir -p "$d"
            fq="$OUT/extracted/$sample/${taxon}_${profile}.fastq"
            n=$(awk -F '\t' -v t="$taxon" -v p="$profile" 'NR>1 && $2==t && $3==p {print $4}' "$OUT/qc/${sample}_extraction.tsv")
            [[ "$n" =~ ^[0-9]+$ ]] || exit 1
            bam="$d/primary.sorted.bam"
            md="$d/mapDamage"
            if [[ "$n" == 0 ]]; then
                printf '%s\t%s\t%s\t%s\t%s\t%s\t0\t0\t0\tNA\t0\tNA\t%s\tNO_READS\tNA\tNA\n' \
                    "$sample" "$taxon" "$taxid" "$rid" "$ref" "$profile" "$MIN_MAPQ" >> "$SUMMARY"
                continue
            fi
            log "$sample / $rid / $profile : $n reads"
            bwa aln -n 0.08 -l 24 -k 2 -q 20 -t "$THREADS" "$localref" "$fq" > "$d/reads.sai" 2> "$d/bwa_aln.log"
            # Exclut secondaires (256) et supplementaires (2048), conserve unmapped.
            bwa samse "$localref" "$d/reads.sai" "$fq" 2> "$d/bwa_samse.log" \
                | samtools view -u -F 2304 - \
                | samtools sort -@ "$SORT_THREADS" -m "$SORT_MEM" -T "$d/sorttmp" -o "$bam" -
            samtools quickcheck "$bam"
            samtools index "$bam"
            samtools flagstat "$bam" > "$d/flagstat.txt"
            samtools stats "$bam" > "$d/samtools_stats.txt"
            total=$(samtools view -c "$bam")
            mapped=$(samtools view -c -F 4 "$bam")
            filtered=$(samtools view -c -F 4 -q "$MIN_MAPQ" "$bam")
            [[ "$total" -eq "$n" ]] || { echo 'Effectif BAM different de FASTQ'; exit 1; }
            rate=$(awk -v m="$mapped" -v n="$n" 'BEGIN {printf "%.6f",100*m/n}')
            qrate=$(awk -v m="$filtered" -v n="$n" 'BEGIN {printf "%.6f",100*m/n}')
            if [[ "$filtered" -eq 0 ]]; then
                status=NO_MAPPED_READS
                md=NA
            else
                mdbam="$d/mapdamage_input.bam"
                samtools view -b -F 4 -q "$MIN_MAPQ" -o "$mdbam" "$bam"
                samtools index "$mdbam"
                if mapDamage -i "$mdbam" -r "$localref" --folder "$md" --no-stats > "$d/mapDamage.log" 2>&1; then
                    status=OK
                else
                    status=FAILED
                    failures=$((failures+1))
                    log "ECHEC mapDamage : $sample / $rid / $profile (voir $d/mapDamage.log)"
                fi
            fi
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "$sample" "$taxon" "$taxid" "$rid" "$ref" "$profile" "$n" "$total" "$mapped" "$rate" \
                "$filtered" "$qrate" "$MIN_MAPQ" "$status" "$bam" "$md" >> "$SUMMARY"
            rm -f "$d/reads.sai"
        done
    done
done < "$OUT/qc/targets.tsv"
if [[ "$KEEP_EXTRACTED" == 0 ]]; then
    find "$OUT/extracted" -type f -name '*.fastq' -delete
fi
if [[ "$failures" -gt 0 ]]; then
    printf '%s echecs mapDamage\n' "$failures" > "$OUT/COMPLETED_WITH_ERRORS"
    log "Termine avec $failures echecs mapDamage. Table : $SUMMARY"
    exit 1
fi
date -Is > "$OUT/SUCCESS"
log "Termine. Table : $SUMMARY"
