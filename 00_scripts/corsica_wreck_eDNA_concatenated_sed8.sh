#!/bin/bash

#SBATCH --job-name=corsica_wreck_eDNA_concatenated_sed8
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/seda_DNA_Corsican_wreck_comparison/00_scripts/corsica_wreck_eDNA_concatenated_sed8.err"
#SBATCH --output="/home/plstenge/seda_DNA_Corsican_wreck_comparison/00_scripts/corsica_wreck_eDNA_concatenated_sed8.out"

################################################################################
# Pipeline aDNA - Projet Corsica Wreck (concaténation par échantillon)
# Author: Pierre-Louis Stenger
# Date: November 2025
################################################################################

set -eo pipefail

################################################################################
# CONFIGURATION
################################################################################

BASE_DIR="/home/plstenge/seda_DNA_Corsican_wreck_comparison"
BBDUK="/home/plstenge/bbmap/bbduk.sh"
CLUMPIFY="/home/plstenge/bbmap/clumpify.sh"
PHIX="/home/plstenge/bbmap/resources/phix174_ill.ref.fa.gz"
KRAKEN2_DB="/home/plstenge/k2_core_nt_20250609"
KRAKENTOOLS_DIR="/home/plstenge/coprolites/07_kraken2/KrakenTools"
THREADS=36

RAW_HOME="/home/plstenge/seda_DNA_Corsican_wreck/01_raw_data"
RUN1="/storage/groups/gdec/shared_paleo/E1531_final/run1_20250320_AV241601_E1531_Ps5Lane1_Ps6Lane2"
RUN2="/storage/groups/gdec/shared_paleo/E1531_final/run2_20250414_AV241601_E1531_Ps5_Ps6_14042025"
RUN3="/storage/groups/gdec/shared_paleo/E1531_final/run3_20251008_AV241601_E1531_Ps7_Ps8"
RUN4="/storage/groups/gdec/shared_paleo/E1531_final/run4_20251104_AV241601_E1531_Ps7_Ps8_04112025"

SAMPLES=(sed8)

################################################################################
# ACTIVATION ENVIRONNEMENT CONDA
################################################################################

echo "Activation environnement conda..."
module load conda/4.12.0
source ~/.bashrc
conda activate metagenomics

################################################################################
# CRÉATION ARBORESCENCE
################################################################################

for sample in "${SAMPLES[@]}"; do
    mkdir -p "${BASE_DIR}/01_raw_data/${sample}"
    mkdir -p "${BASE_DIR}/02_quality_check_raw/${sample}"
    mkdir -p "${BASE_DIR}/03_bbduk/${sample}"
    mkdir -p "${BASE_DIR}/04_fastuniq/${sample}"
    mkdir -p "${BASE_DIR}/05_clumpify/${sample}"
    mkdir -p "${BASE_DIR}/06_fastp/${sample}"
    mkdir -p "${BASE_DIR}/07_quality_check_clean/${sample}"
    mkdir -p "${BASE_DIR}/08_kraken2/${sample}"
    mkdir -p "${BASE_DIR}/09_krona/${sample}"
    mkdir -p "${BASE_DIR}/10_mpa_tables/${sample}"
    mkdir -p "${BASE_DIR}/12_mapdamage/${sample}"
    mkdir -p "${BASE_DIR}/12_mapdamage/${sample}"
    mkdir -p "${BASE_DIR}/12_mapdamage/${sample}"
    mkdir -p "${BASE_DIR}/00_scripts"
    mkdir -p "${BASE_DIR}/11_summary_tables"
    mkdir -p "${BASE_DIR}/12_mapdamage"
    mkdir -p "${BASE_DIR}/13_final"
    echo "Arborescence créée pour $sample"
done

################################################################################
# CONCATÉNATION DES FICHIERS PAR ÉCHANTILLON
################################################################################


#echo "Concaténation des fichiers sed6..."

# sed6 - R1
#cat \
#  ${RAW_HOME}/1120_sed6_rep3_R1.fastq.gz \
#  ${RUN1}/1120_sed6_rep3/1120_sed6_rep3_R1.fastq.gz \
#  ${RUN2}/1120_sed6_rep3/1120_sed6_rep3_R1.fastq.gz \
#  ${RAW_HOME}/1129_sed6_rep1_R1.fastq.gz \
#  ${RUN1}/1129_sed6_rep1/1129_sed6_rep1_R1.fastq.gz \
#  ${RUN2}/1129_sed6_rep1/1129_sed6_rep1_R1.fastq.gz \
#  ${RUN3}/1129_sed6_rep1/1129_sed6_rep1_R1.fastq.gz \
#  ${RUN4}/1129_sed6_rep1/1129_sed6_rep1_R1.fastq.gz \
#  ${RAW_HOME}/1130_sed6_rep2_R1.fastq.gz \
#  ${RUN1}/1130_sed6_rep2/1130_sed6_rep2_R1.fastq.gz \
#  ${RUN2}/1130_sed6_rep2/1130_sed6_rep2_R1.fastq.gz \
#  ${RUN3}/1130_sed6_rep2/1130_sed6_rep2_R1.fastq.gz \
#  ${RUN4}/1130_sed6_rep2/1130_sed6_rep2_R1.fastq.gz \
#  > "${BASE_DIR}/01_raw_data/sed6/sed6_concat_R1.fastq.gz"
#
## sed6 - R2
#cat \
#  ${RAW_HOME}/1120_sed6_rep3_R2.fastq.gz \
#  ${RUN1}/1120_sed6_rep3/1120_sed6_rep3_R2.fastq.gz \
#  ${RUN2}/1120_sed6_rep3/1120_sed6_rep3_R2.fastq.gz \
#  ${RAW_HOME}/1129_sed6_rep1_R2.fastq.gz \
#  ${RUN1}/1129_sed6_rep1/1129_sed6_rep1_R2.fastq.gz \
#  ${RUN2}/1129_sed6_rep1/1129_sed6_rep1_R2.fastq.gz \
#  ${RUN3}/1129_sed6_rep1/1129_sed6_rep1_R2.fastq.gz \
#  ${RUN4}/1129_sed6_rep1/1129_sed6_rep1_R2.fastq.gz \
#  ${RAW_HOME}/1130_sed6_rep2_R2.fastq.gz \
#  ${RUN1}/1130_sed6_rep2/1130_sed6_rep2_R2.fastq.gz \
#  ${RUN2}/1130_sed6_rep2/1130_sed6_rep2_R2.fastq.gz \
#  ${RUN3}/1130_sed6_rep2/1130_sed6_rep2_R2.fastq.gz \
#  ${RUN4}/1130_sed6_rep2/1130_sed6_rep2_R2.fastq.gz \
#  > "${BASE_DIR}/01_raw_data/sed6/sed6_concat_R2.fastq.gz"
#
#echo "Concaténation sed6 terminée."
#
#echo "Concaténation des fichiers sed8..."
#
## sed8 - R1
#cat \
#  ${RAW_HOME}/1121_sed8_rep1_R1.fastq.gz \
#  ${RUN1}/1121_sed8_rep1/1121_sed8_rep1_R1.fastq.gz \
#  ${RUN2}/1121_sed8_rep1/1121_sed8_rep1_R1.fastq.gz \
#  ${RUN3}/1121_sed8_rep1/1121_sed8_rep1_R1.fastq.gz \
#  ${RUN4}/1121_sed8_rep1/1121_sed8_rep1_R1.fastq.gz \
#  ${RAW_HOME}/1122_sed8_rep2_R1.fastq.gz \
#  ${RUN1}/1122_sed8_rep2/1122_sed8_rep2_R1.fastq.gz \
#  ${RUN2}/1122_sed8_rep2/1122_sed8_rep2_R1.fastq.gz \
#  ${RUN3}/1122_sed8_rep2/1122_sed8_rep2_R1.fastq.gz \
#  ${RUN4}/1122_sed8_rep2/1122_sed8_rep2_R1.fastq.gz \
#  ${RAW_HOME}/1131_sed8_rep3_R1.fastq.gz \
#  ${RUN1}/1131_sed8_rep3/1131_sed8_rep3_R1.fastq.gz \
#  ${RUN2}/1131_sed8_rep3/1131_sed8_rep3_R1.fastq.gz \
#  > "${BASE_DIR}/01_raw_data/sed8/sed8_concat_R1.fastq.gz"
#
## sed8 - R2
#cat \
#  ${RAW_HOME}/1121_sed8_rep1_R2.fastq.gz \
#  ${RUN1}/1121_sed8_rep1/1121_sed8_rep1_R2.fastq.gz \
#  ${RUN2}/1121_sed8_rep1/1121_sed8_rep1_R2.fastq.gz \
#  ${RUN3}/1121_sed8_rep1/1121_sed8_rep1_R2.fastq.gz \
#  ${RUN4}/1121_sed8_rep1/1121_sed8_rep1_R2.fastq.gz \
#  ${RAW_HOME}/1122_sed8_rep2_R2.fastq.gz \
#  ${RUN1}/1122_sed8_rep2/1122_sed8_rep2_R2.fastq.gz \
#  ${RUN2}/1122_sed8_rep2/1122_sed8_rep2_R2.fastq.gz \
#  ${RUN3}/1122_sed8_rep2/1122_sed8_rep2_R2.fastq.gz \
#  ${RUN4}/1122_sed8_rep2/1122_sed8_rep2_R2.fastq.gz \
#  ${RAW_HOME}/1131_sed8_rep3_R2.fastq.gz \
#  ${RUN1}/1131_sed8_rep3/1131_sed8_rep3_R2.fastq.gz \
#  ${RUN2}/1131_sed8_rep3/1131_sed8_rep3_R2.fastq.gz \
#  > "${BASE_DIR}/01_raw_data/sed8/sed8_concat_R2.fastq.gz"
#
#echo "Concaténation sed8 terminée."
#echo "Concaténation terminée pour tous les échantillons."
#
#################################################################################
## QUALITÉ + NETTOYAGE (FastQC, MultiQC, BBDuk, FastUniq, Clumpify, Fastp)
#################################################################################
#for sample in "${SAMPLES[@]}"; do
#  # 1. Contrôle qualité RAW
#  fastqc ${BASE_DIR}/01_raw_data/${sample}/${sample}_concat_R1.fastq.gz ${BASE_DIR}/01_raw_data/${sample}/${sample}_concat_R2.fastq.gz -o ${BASE_DIR}/02_quality_check_raw/${sample} -t 4
#  multiqc ${BASE_DIR}/02_quality_check_raw/${sample} -o ${BASE_DIR}/02_quality_check_raw/${sample} --force
#
#  # 2. Filtrage/adaptateurs BBDuk
#  $BBDUK in1=${BASE_DIR}/01_raw_data/${sample}/${sample}_concat_R1.fastq.gz in2=${BASE_DIR}/01_raw_data/${sample}/${sample}_concat_R2.fastq.gz out1=${BASE_DIR}/03_bbduk/${sample}/${sample}_bbduk_R1.fastq.gz out2=${BASE_DIR}/03_bbduk/${sample}/${sample}_bbduk_R2.fastq.gz ref=$PHIX ktrim=r k=23 mink=11 hdist=1 tpe tbo minlen=25 qtrim=r trimq=20 stats=${BASE_DIR}/03_bbduk/${sample}/${sample}_bbduk_stats.txt
#
#  # 3. Déduplication FastUniq
#  echo "FastUniq pour ${sample}..."
#  mkdir -p ${BASE_DIR}/04_fastuniq/${sample}/tmp
#  # Décompresser les fichiers BBDuk
#  zcat ${BASE_DIR}/03_bbduk/${sample}/${sample}_bbduk_R1.fastq.gz > ${BASE_DIR}/04_fastuniq/${sample}/tmp/${sample}_bbduk_R1.fastq
#  zcat ${BASE_DIR}/03_bbduk/${sample}/${sample}_bbduk_R2.fastq.gz > ${BASE_DIR}/04_fastuniq/${sample}/tmp/${sample}_bbduk_R2.fastq
#  # Créer le fichier liste pour FastUniq (UN FICHIER PAR LIGNE, pas de tabulation)
#  echo "${BASE_DIR}/04_fastuniq/${sample}/tmp/${sample}_bbduk_R1.fastq" > ${BASE_DIR}/04_fastuniq/${sample}/tmp/infile.list
#  echo "${BASE_DIR}/04_fastuniq/${sample}/tmp/${sample}_bbduk_R2.fastq" >> ${BASE_DIR}/04_fastuniq/${sample}/tmp/infile.list
#  # Exécuter FastUniq
#  fastuniq -i ${BASE_DIR}/04_fastuniq/${sample}/tmp/infile.list -t q -o ${BASE_DIR}/04_fastuniq/${sample}/${sample}_fastuniq_R1.fastq -p ${BASE_DIR}/04_fastuniq/${sample}/${sample}_fastuniq_R2.fastq
#  # Compresser les résultats
#  gzip ${BASE_DIR}/04_fastuniq/${sample}/${sample}_fastuniq_R1.fastq
#  gzip ${BASE_DIR}/04_fastuniq/${sample}/${sample}_fastuniq_R2.fastq
#  # Nettoyer les fichiers temporaires
#  rm -rf ${BASE_DIR}/04_fastuniq/${sample}/tmp
#  echo "FastUniq terminé pour ${sample}."
#  
#  # 4. Déduplication optique Clumpify
#  $CLUMPIFY in=${BASE_DIR}/04_fastuniq/${sample}/${sample}_fastuniq_R1.fastq.gz in2=${BASE_DIR}/04_fastuniq/${sample}/${sample}_fastuniq_R2.fastq.gz out=${BASE_DIR}/05_clumpify/${sample}/${sample}_clumpify_R1.fastq.gz out2=${BASE_DIR}/05_clumpify/${sample}/${sample}_clumpify_R2.fastq.gz dedupe=t
#
#  # 5. Nettoyage et merge Fastp
#  fastp -i ${BASE_DIR}/05_clumpify/${sample}/${sample}_clumpify_R1.fastq.gz -I ${BASE_DIR}/05_clumpify/${sample}/${sample}_clumpify_R2.fastq.gz --merged_out ${BASE_DIR}/06_fastp/${sample}/${sample}_fastp_merged.fastq.gz --out1 ${BASE_DIR}/06_fastp/${sample}/${sample}_fastp_R1.fastq.gz --out2 ${BASE_DIR}/06_fastp/${sample}/${sample}_fastp_R2.fastq.gz --json ${BASE_DIR}/06_fastp/${sample}/${sample}_fastp.json --html ${BASE_DIR}/06_fastp/${sample}/${sample}_fastp.html --thread 4 --length_required 30 --qualified_quality_phred 20
#
#  # 6. Contrôle qualité Clean
#  fastqc ${BASE_DIR}/06_fastp/${sample}/${sample}_fastp_*.fastq.gz -o ${BASE_DIR}/07_quality_check_clean/${sample} -t 4
#  multiqc ${BASE_DIR}/07_quality_check_clean/${sample} -o ${BASE_DIR}/07_quality_check_clean/${sample} --force
#
#done
#
#
#################################################################################
## KRAKEN2: CLASSIFICATION TAXONOMIQUE
#################################################################################
#echo "Classification taxonomique Kraken2..."
#
#for sample in "${SAMPLES[@]}"; do
#  echo "Kraken2 pour ${sample}..."
#  FASTPDIR="${BASE_DIR}/06_fastp/${sample}"
#  OUTDIR="${BASE_DIR}/08_kraken2/${sample}"
#  
#  # Analyse des reads merged
#  MERGED="${FASTPDIR}/${sample}_fastp_merged.fastq.gz"
#  if [ -f "$MERGED" ]; then
#    OUTKRAKEN="${OUTDIR}/${sample}_merged.kraken"
#    OUTREPORT="${OUTDIR}/${sample}_merged.report"
#    echo "  ${sample} merged..."
#    kraken2 --confidence 0.2 \
#      --db ${KRAKEN2_DB} \
#      --threads ${THREADS} \
#      --output ${OUTKRAKEN} \
#      --report ${OUTREPORT} \
#      ${MERGED}
#  fi
#  
#  # Analyse des reads unmerged (paired)
#  R1="${FASTPDIR}/${sample}_fastp_R1.fastq.gz"
#  R2="${FASTPDIR}/${sample}_fastp_R2.fastq.gz"
#  if [ -f "$R1" ] && [ -f "$R2" ]; then
#    OUTKRAKEN="${OUTDIR}/${sample}_unmerged.kraken"
#    OUTREPORT="${OUTDIR}/${sample}_unmerged.report"
#    echo "  ${sample} unmerged..."
#    kraken2 --confidence 0.2 \
#      --paired \
#      --db ${KRAKEN2_DB} \
#      --threads ${THREADS} \
#      --output ${OUTKRAKEN} \
#      --report ${OUTREPORT} \
#      ${R1} ${R2}
#  fi
#done
#
#echo "Classification Kraken2 terminée."
#
#################################################################################
## KRONA: VISUALISATION
#################################################################################
#echo "Visualisation Krona..."
#
#CONDA_PREFIX=$(conda info --base)/envs/metagenomics
#KRONA_TAX_DIR="${CONDA_PREFIX}/opt/krona/taxonomy"
#
#if [ ! -d "${KRONA_TAX_DIR}" ] || [ ! -f "${KRONA_TAX_DIR}/taxonomy.tab" ]; then
#  echo "Taxonomie Krona absente. Installation en cours..."
#  ktUpdateTaxonomy.sh "${KRONA_TAX_DIR}"
#  echo "Taxonomie Krona installée avec succès."
#else
#  echo "Taxonomie Krona déjà installée."
#fi
#
#INDIR="${BASE_DIR}/08_kraken2/${SAMPLE}${SUFFIX}"
#OUTDIR="${BASE_DIR}/09_krona/${SAMPLE}${SUFFIX}"
#
#if [ -d "$INDIR" ] && ls ${INDIR}/*.report > /dev/null 2>&1; then
#  cd "${INDIR}"
#  ktImportTaxonomy *.report -o "${OUTDIR}/${SAMPLE}${SUFFIX}_krona.html"
#  echo "Krona HTML généré pour ${SAMPLE}${SUFFIX}"
#else
#  echo "Aucun fichier report trouvé pour ${SAMPLE}${SUFFIX}"
#fi
#
#
#################################################################################
## ÉTAPE 9: Créer des tables MPA à partir de Kraken2 reports
#################################################################################
#
#for sample in "${SAMPLES[@]}"; do
#  KRAKENDIR="${BASE_DIR}/08_kraken2/${sample}"
#  MPA_DIR="${BASE_DIR}/10_mpa_tables/${sample}"
#  mkdir -p "$MPA_DIR"
#  declare -a mpafiles
#
#  for report in ${KRAKENDIR}/*.report; do
#    if [ -f "$report" ]; then
#      base=$(basename "$report" .report)
#      mpafile="${MPA_DIR}/${base}.mpa"
#      python3 ${KRAKENTOOLS_DIR}/kreport2mpa.py -r "$report" -o "$mpafile"
#      mpafiles+=("$mpafile")
#    fi
#  done
#
#  if [ ${#mpafiles[@]} -gt 0 ]; then
#    python3 ${KRAKENTOOLS_DIR}/combine_mpa.py -i ${mpafiles[*]} -o ${MPA_DIR}/combined_${sample}.tsv
#  fi
#  echo "Tables MPA générées pour $sample."
#done
#
#################################################################################
## ÉTAPE 10: Créer un tableau récapitulatif des séquences après chaque étape clé
#################################################################################
#
#SUMMARY_TABLE="${BASE_DIR}/11_summary_tables/sequences_summary.tsv"
#echo -e "Sample\tStage\tNbSequences\tAvgLength\tGCpercent" > $SUMMARY_TABLE
## Fonction pour extraire stats d'un fastq.gz
#extract_fastq_stats() {
#  local fq=$1
#  local sample=$2
#  local stage=$3
#  if [ ! -f "$fq" ]; then return; fi
#  # Nombre de reads
#  nb=$(zcat "$fq" | echo $((`wc -l`/4)))
#  # Moyenne de longueur
#  avg=$(zcat "$fq" | awk 'NR%4==2 {sum+=length($1); n++} END {if(n>0) print sum/n; else print 0}')
#  # %GC
#  gc=$(zcat "$fq" | awk 'NR%4==2 {gc+=gsub(/[GCgc]/,"&") ; t+=length($1)} END {if(t>0) printf "%.2f", 100*gc/t; else print 0}')
#  echo -e "$sample\t$stage\t$nb\t$avg\t$gc" >> $SUMMARY_TABLE
#}
#
## Extraction pour chaque sample après étapes majeures
#for sample in "${SAMPLES[@]}"; do
#  RAW_DIR="${BASE_DIR}/01_raw_data/${sample}"
#  FASTP_DIR="${BASE_DIR}/06_fastp/${sample}"
#  extract_fastq_stats "${RAW_DIR}/${sample}_concat_R1.fastq.gz" "$sample" "RAW_R1"
#  extract_fastq_stats "${RAW_DIR}/${sample}_concat_R2.fastq.gz" "$sample" "RAW_R2"
#  extract_fastq_stats "${FASTP_DIR}/${sample}_fastp_R1.fastq.gz" "$sample" "Clean_R1"
#  extract_fastq_stats "${FASTP_DIR}/${sample}_fastp_R2.fastq.gz" "$sample" "Clean_R2"
#  extract_fastq_stats "${FASTP_DIR}/${sample}_fastp_merged.fastq.gz" "$sample" "Clean_Merged"
#done
#
#echo "Tableau récapitulatif généré: $SUMMARY_TABLE"

################################################################################
# ÉTAPE 11: MAPDAMAGE - ANALYSE DES DOMMAGES DE L'ADN ANCIEN
################################################################################

echo ""
echo "======================================================================"
echo "ÉTAPE 11: MapDamage - Analyse des dommages de l'ADN"
echo "Date de début: $(date)"
echo "======================================================================"
echo ""

KRAKENBASE="${BASE_DIR}/08_kraken2"
FASTQBASE="${BASE_DIR}/06_fastp"
DAMAGEBASE="${BASE_DIR}/12_mapdamage"
LOGFILE="${BASE_DIR}/00_scripts/mapdamage_$(date +%Y%m%d_%H%M%S).txt"
MAPPINGINFO="${BASE_DIR}/11_summary_tables/mapping_bwa_info.tsv"

mkdir -p "${DAMAGEBASE}"

# Activer l'environnement MapDamage
#module load conda/4.12.0
#source ~/.bashrc
conda activate mapdamage_py39

echo "Script MapDamage started at $(date)" | tee -a "${LOGFILE}"

# Initialiser le fichier de mapping info
echo -e "Sample\tSpecies\tType\tTotalReads\tMappedReads\tMappingRate" > "${MAPPINGINFO}"

################################################################################
# DÉFINITION DES GÉNOMES DE RÉFÉRENCE (TAXONS MODIFIÉS)
################################################################################
declare -A TAXONS=(
  ["Vitis_vinifera"]="29760:/storage/groups/gdec/shared_paleo/genomes_REF/12Xv2_grapevine_genome_assembly.fa"
 # ["Homo_sapiens"]="9606:/storage/biodatabanks/ucsc/genomes/hg19/Homo_sapiens-hg19_2012-9-19/fasta/all.fasta"
  ["Mus_musculus"]="10090:/home/plstenge/genomes/Mus_musculus.GRCm39.dna.toplevel.fa"
  ["Melanogrammus_aeglefinus"]="8056:/home/plstenge/genomes/Melanogrammus_aeglefinus_OLKM01.fasta"
  ["Gobiusculus_flavescens"]="257540:/home/plstenge/genomes/Gobiusculus_flavescens_fGobFla1.fasta"
)

################################################################################
# INDEXATION BWA DES GÉNOMES DE RÉFÉRENCE
################################################################################
echo "Vérification des index BWA..." # A ne faire qu'une fois !!
#bwa index /storage/groups/gdec/shared_paleo/genomes_REF/12Xv2_grapevine_genome_assembly.fa
##bwa index /storage/biodatabanks/ucsc/genomes/hg19/Homo_sapiens-hg19_2012-9-19/fasta/all.fasta
#bwa index /home/plstenge/genomes/Mus_musculus.GRCm39.dna.toplevel.fa
#bwa index /home/plstenge/genomes/Melanogrammus_aeglefinus_OLKM01.fasta
#bwa index /home/plstenge/genomes/Gobiusculus_flavescens_fGobFla1.fasta

################################################################################
# FONCTION POUR CALCULER LE TAUX DE MAPPING
################################################################################
calculate_mapping_rate() {
  local bamfile=$1
  local samplename=$2
  local species=$3
  local type=$4

  if [ -f "$bamfile" ]; then
    local totalreads=$(samtools view -c "$bamfile")
    local mappedreads=$(samtools view -c -F 4 "$bamfile")
    local mappingrate=0
    if [ "$totalreads" -gt 0 ]; then
      mappingrate=$(echo "scale=2; $mappedreads * 100 / $totalreads" | bc)
    fi
    echo -e "$samplename\t$species\t$type\t$totalreads\t$mappedreads\t$mappingrate" >> "${MAPPINGINFO}"
    echo "Mapping stats for $samplename/$species/$type: $mappedreads/$totalreads ($mappingrate%)" | tee -a "${LOGFILE}"
  fi
}

################################################################################
# BOUCLE SUR LES ÉCHANTILLONS (sed6 et sed8)
################################################################################
shopt -s nullglob

for sample in "${SAMPLES[@]}"; do
  echo ""
  echo "======================================================================"
  echo "Traitement de l'échantillon: $sample"
  echo "======================================================================"

  KRAKENDIR="${KRAKENBASE}/${sample}"
  FASTQDIR="${FASTQBASE}/${sample}"

  if [ ! -d "$KRAKENDIR" ]; then
    echo "ATTENTION: Répertoire Kraken2 absent pour $sample" | tee -a "${LOGFILE}"
    continue
  fi

  # Boucle sur les fichiers Kraken2 (merged et unmerged)
  for KRAKENFILE in ${KRAKENDIR}/*.kraken; do
    if [ ! -f "$KRAKENFILE" ]; then
      continue
    fi

    KRAKENBASENAME=$(basename "$KRAKENFILE" .kraken)
    echo ""
    echo ">>> Processing $KRAKENBASENAME ($sample)" | tee -a "${LOGFILE}"

    # Extraire le préfixe de base
    PREFIX=$(echo "$KRAKENBASENAME" | sed -E 's/(un)?merged$//')
    echo "Prefix: $PREFIX" | tee -a "${LOGFILE}"

    # Chercher les fichiers FASTQ correspondants
    R1FILE="${FASTQDIR}/${sample}_fastp_R1.fastq.gz"
    R2FILE="${FASTQDIR}/${sample}_fastp_R2.fastq.gz"
    MERGEDFILE="${FASTQDIR}/${sample}_fastp_merged.fastq.gz"

    ############################################################################
    # BOUCLE SUR LES ESPÈCES (5 TAXONS)
    ############################################################################
    for GROUP in "${!TAXONS[@]}"; do
      IFS=':' read -r TAXID REFFASTA <<< "${TAXONS[$GROUP]}"
      DAMAGEDIR="${DAMAGEBASE}/${sample}/${GROUP}"
      mkdir -p "${DAMAGEDIR}"

      echo ""
      echo "--- Espèce: $GROUP (TaxID: $TAXID) ---"

      OUTR1="${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}_R1.fastq"
      OUTR2="${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}_R2.fastq"
      OUTMERGED="${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}_merged.fastq"

      ########################################################################
      # TRAITEMENT DES READS UNMERGED (PAIRED-END)
      ########################################################################
      if [[ "$KRAKENBASENAME" == *"unmerged"* ]] && [ -f "$R1FILE" ] && [ -f "$R2FILE" ]; then
        echo "Extraction des reads unmerged pour $GROUP..." | tee -a "${LOGFILE}"
        python3 ${KRAKENTOOLS_DIR}/extract_kraken_reads.py \
          -k "$KRAKENFILE" \
          -s "$R1FILE" -s2 "$R2FILE" \
          -t "$TAXID" \
          -o "$OUTR1" -o2 "$OUTR2" \
          --fastq-output 2>>"${LOGFILE}"
  
        if [ -f "$OUTR1" ] && [ -f "$OUTR2" ]; then
          echo "Mapping BWA paired-end pour $GROUP..." | tee -a "${LOGFILE}"
    
          # BWA aln
          bwa aln -n 0.08 -l 24 -k 2 -q 20 -t 4 "$REFFASTA" "$OUTR1" > "${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}_R1.sai" 2>>"${LOGFILE}"
          bwa aln -n 0.08 -l 24 -k 2 -q 20 -t 4 "$REFFASTA" "$OUTR2" > "${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}_R2.sai" 2>>"${LOGFILE}"
    
          # BWA sampe
          bwa sampe "$REFFASTA" \
            "${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}_R1.sai" \
            "${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}_R2.sai" \
            "$OUTR1" "$OUTR2" \
            > "${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}.sam" 2>>"${LOGFILE}"
    
          # Conversion SAM → BAM
          samtools view -bS "${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}.sam" > "${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}.bam" 2>>"${LOGFILE}"
    
          # Tri et indexation
          samtools sort -o "${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}.sorted.bam" "${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}.bam" 2>>"${LOGFILE}"
          samtools index "${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}.sorted.bam" 2>>"${LOGFILE}"
    
          # MapDamage
          echo "MapDamage unmerged pour $GROUP..." | tee -a "${LOGFILE}"
          mapDamage -i "${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}.sorted.bam" \
            -r "$REFFASTA" \
            --folder "${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}_mapDamage_unmerged" 2>>"${LOGFILE}"
    
          # Calcul du taux de mapping
          calculate_mapping_rate "${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}.sorted.bam" "$sample" "$GROUP" "unmerged"
    
          # Nettoyage des fichiers intermédiaires
          rm -f "${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}_R1.sai" \
                "${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}_R2.sai" \
                "${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}.sam" \
                "${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}.bam" 2>>"${LOGFILE}"
        fi
      fi

      ########################################################################
      # TRAITEMENT DES READS MERGED (SINGLE-END)
      ########################################################################
      if [[ "$KRAKENBASENAME" == *"merged"* ]] && [ -f "$MERGEDFILE" ]; then
        echo "Extraction des reads merged pour $GROUP..." | tee -a "${LOGFILE}"
        python3 ${KRAKENTOOLS_DIR}/extract_kraken_reads.py \
          -k "$KRAKENFILE" \
          -s "$MERGEDFILE" \
          -t "$TAXID" \
          -o "$OUTMERGED" \
          --fastq-output 2>>"${LOGFILE}"
  
        if [ -f "$OUTMERGED" ]; then
          echo "Mapping BWA single-end pour $GROUP..." | tee -a "${LOGFILE}"
    
          # BWA aln
          bwa aln -n 0.08 -l 24 -k 2 -q 20 -t 4 "$REFFASTA" "$OUTMERGED" > "${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}_merged.sai" 2>>"${LOGFILE}"
    
          # BWA samse
          bwa samse "$REFFASTA" \
            "${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}_merged.sai" \
            "$OUTMERGED" \
            > "${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}_merged.sam" 2>>"${LOGFILE}"
    
          # Conversion SAM → BAM
          samtools view -bS "${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}_merged.sam" > "${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}_merged.bam" 2>>"${LOGFILE}"
    
          # Tri et indexation
          samtools sort -o "${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}_merged.sorted.bam" "${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}_merged.bam" 2>>"${LOGFILE}"
          samtools index "${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}_merged.sorted.bam" 2>>"${LOGFILE}"
    
          # MapDamage
          echo "MapDamage merged pour $GROUP..." | tee -a "${LOGFILE}"
          mapDamage -i "${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}_merged.sorted.bam" \
            -r "$REFFASTA" \
            --folder "${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}_mapDamage_merged" 2>>"${LOGFILE}"
    
          # Calcul du taux de mapping
          calculate_mapping_rate "${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}_merged.sorted.bam" "$sample" "$GROUP" "merged"
    
          # Nettoyage des fichiers intermédiaires
          rm -f "${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}_merged.sai" \
                "${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}_merged.sam" \
                "${DAMAGEDIR}/${KRAKENBASENAME}_${GROUP}_merged.bam" 2>>"${LOGFILE}"
        fi
      fi

    done  # Fin boucle sur les espèces
  done  # Fin boucle sur les fichiers Kraken
done  # Fin boucle sur les échantillons

################################################################################
# FIN MAPDAMAGE
################################################################################
echo ""
echo "======================================================================"
echo "MapDamage terminé!"
echo "Date de fin: $(date)"
echo ""
echo "Résultats MapDamage: ${DAMAGEBASE}"
echo "Statistiques de mapping: ${MAPPINGINFO}"
echo "======================================================================"
echo ""

# Notification par email (optionnel)
echo "MapDamage terminé le $(date). Résultats: ${DAMAGEBASE}" | \
  mail -s "Pipeline Corsica Wreck - MapDamage terminé" pierrelouis.stenger@gmail.com 2>/dev/null || true
