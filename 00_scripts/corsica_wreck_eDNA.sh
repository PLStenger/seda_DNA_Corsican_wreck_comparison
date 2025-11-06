#!/bin/bash

#SBATCH --job-name=corsica_wreck_eDNA
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/seda_DNA_Corsican_wreck_comparison/00_scripts/corsica_wreck_eDNA.err"
#SBATCH --output="/home/plstenge/seda_DNA_Corsican_wreck_comparison/00_scripts/corsica_wreck_eDNA.out"

################################################################################
# Pipeline d'analyse sedaDNA - Comparaison de recettes de séquençage
# Project: Corsica wreck
# Author: Pierre-Louis Stenger
# Date: November 2025
#
# Description:
# Ce pipeline analyse 4 échantillons sedaDNA séquencés avec 2 recettes:
# - Recipe1 (R1): recipes-ShortInsert (2x150bp) - séquençage standard
# - Recipe2 (R2): recipes-v3.3.x-ShortInsert-SmallFragmentRecovery (2x150bp) - récupération fragments courts
# - Combined (R1R2): Fusion des deux recettes pour maximiser la couverture
#
# Échantillons:
# - 1121_sed8_rep1 (sed8, réplicat 1)
# - 1122_sed8_rep2 (sed8, réplicat 2)
# - 1129_sed6_rep1 (sed6, réplicat 1)
# - 1130_sed6_rep2 (sed6, réplicat 2)
################################################################################

set -eo pipefail

# Définition des chemins de base
BASE_DIR="/home/plstenge/seda_DNA_Corsican_wreck_comparison"
RECIPE1_DIR="/home/plstenge/seda_DNA_Corsican_wreck/01_raw_data"
RECIPE2_DIR="/storage/groups/gdec/shared_paleo/E1531_final/run3_20251008_AV241601_E1531_Ps7_Ps8"

# Outils
BBDUK="/home/plstenge/bbmap/bbduk.sh"
CLUMPIFY="/home/plstenge/bbmap/clumpify.sh"
PHIX="/home/plstenge/bbmap/resources/phix174_ill.ref.fa.gz"
KRAKEN2_DB="/home/plstenge/k2_core_nt_20250609"
KRAKENTOOLS_DIR="/home/plstenge/coprolites/07_kraken2/KrakenTools"

# Paramètres
THREADS=36

# Définition des échantillons
declare -a SAMPLES=("1121_sed8_rep1" "1122_sed8_rep2" "1129_sed6_rep1" "1130_sed6_rep2")

echo "=========================================="
echo "Pipeline sedaDNA - Corsica wreck"
echo "Date de début: $(date)"
echo "=========================================="

################################################################################
# ÉTAPE 0: Création de l'arborescence et copie des données
################################################################################

echo ""
echo "=== ÉTAPE 0: Préparation des données ==="
echo ""

# Création de l'arborescence
mkdir -p "${BASE_DIR}/00_scripts"
mkdir -p "${BASE_DIR}/01_raw_data/recipe1_standard"
mkdir -p "${BASE_DIR}/01_raw_data/recipe2_smallfrag"
mkdir -p "${BASE_DIR}/01_raw_data/combined_recipe1_recipe2"
mkdir -p "${BASE_DIR}/02_quality_check"
mkdir -p "${BASE_DIR}/03_bbduk"
mkdir -p "${BASE_DIR}/04_fastuniq"
mkdir -p "${BASE_DIR}/05_clumpify"
mkdir -p "${BASE_DIR}/06_fastp"
mkdir -p "${BASE_DIR}/07_kraken2"
mkdir -p "${BASE_DIR}/08_krona"
mkdir -p "${BASE_DIR}/09_mpa_tables"

## Copie des fichiers Recipe 1 (Standard ShortInsert)
#echo "Copie des fichiers Recipe1 (Standard ShortInsert)..."
#for sample in "${SAMPLES[@]}"; do
#    cp "${RECIPE1_DIR}/${sample}_R1.fastq.gz" "${BASE_DIR}/01_raw_data/recipe1_standard/" 2>/dev/null || echo "Warning: ${sample} R1 not found in Recipe1"
#    cp "${RECIPE1_DIR}/${sample}_R2.fastq.gz" "${BASE_DIR}/01_raw_data/recipe1_standard/" 2>/dev/null || echo "Warning: ${sample} R2 not found in Recipe1"
#done
#
## Copie des fichiers Recipe 2 (SmallFragmentRecovery)
#echo "Copie des fichiers Recipe2 (SmallFragmentRecovery)..."
#for sample in "${SAMPLES[@]}"; do
#    cp "${RECIPE2_DIR}/${sample}/${sample}_R1.fastq.gz" "${BASE_DIR}/01_raw_data/recipe2_smallfrag/" 2>/dev/null || echo "Warning: ${sample} R1 not found in Recipe2"
#    cp "${RECIPE2_DIR}/${sample}/${sample}_R2.fastq.gz" "${BASE_DIR}/01_raw_data/recipe2_smallfrag/" 2>/dev/null || echo "Warning: ${sample} R2 not found in Recipe2"
#done
#
## Création des fichiers combinés (Recipe1 + Recipe2)
#echo "Fusion des fichiers Recipe1 + Recipe2..."
#for sample in "${SAMPLES[@]}"; do
#    R1_RECIPE1="${BASE_DIR}/01_raw_data/recipe1_standard/${sample}_R1.fastq.gz"
#    R2_RECIPE1="${BASE_DIR}/01_raw_data/recipe1_standard/${sample}_R2.fastq.gz"
#    R1_RECIPE2="${BASE_DIR}/01_raw_data/recipe2_smallfrag/${sample}_R1.fastq.gz"
#    R2_RECIPE2="${BASE_DIR}/01_raw_data/recipe2_smallfrag/${sample}_R2.fastq.gz"
#    
#    # Vérification que les fichiers existent
#    if [[ -f "$R1_RECIPE1" ]] && [[ -f "$R1_RECIPE2" ]]; then
#        cat "$R1_RECIPE1" "$R1_RECIPE2" > "${BASE_DIR}/01_raw_data/combined_recipe1_recipe2/${sample}_R1.fastq.gz"
#        echo "  → ${sample}_R1 combiné créé"
#    fi
#    
#    if [[ -f "$R2_RECIPE1" ]] && [[ -f "$R2_RECIPE2" ]]; then
#        cat "$R2_RECIPE1" "$R2_RECIPE2" > "${BASE_DIR}/01_raw_data/combined_recipe1_recipe2/${sample}_R2.fastq.gz"
#        echo "  → ${sample}_R2 combiné créé"
#    fi
#done
#
#echo "Préparation des données terminée."
#
#################################################################################
## ÉTAPE 1: Contrôle qualité avec FastQC et MultiQC
#################################################################################
#
#echo ""
#echo "=== ÉTAPE 1: Contrôle qualité (FastQC/MultiQC) ==="
#echo ""
#
#module load conda/4.12.0
#source ~/.bashrc
#conda activate fastqc
#
## FastQC pour chaque type de recette
#for recipe_type in recipe1_standard recipe2_smallfrag combined_recipe1_recipe2; do
#    echo "FastQC pour ${recipe_type}..."
#    RECIPE_DIR="${BASE_DIR}/01_raw_data/${recipe_type}"
#    OUTPUT_DIR="${BASE_DIR}/02_quality_check/${recipe_type}"
#    mkdir -p "$OUTPUT_DIR"
#    
#    for FILE in "${RECIPE_DIR}"/*.fastq.gz; do
#        if [[ -f "$FILE" ]]; then
#            fastqc "$FILE" -o "$OUTPUT_DIR" -t 4
#        fi
#    done
#done
#
#conda deactivate
#
## MultiQC pour résumer tous les rapports
#conda activate multiqc
#for recipe_type in recipe1_standard recipe2_smallfrag combined_recipe1_recipe2; do
#    echo "MultiQC pour ${recipe_type}..."
#    OUTPUT_DIR="${BASE_DIR}/02_quality_check/${recipe_type}"
#    cd "$OUTPUT_DIR"
#    multiqc . -n "multiqc_${recipe_type}.html"
#done
#conda deactivate
#
#echo "Contrôle qualité terminé."

#################################################################################
## ÉTAPE 2: Filtrage et trimming avec BBDuk
#################################################################################
#
#echo ""
#echo "=== ÉTAPE 2: Filtrage et trimming (BBDuk) ==="
#echo ""
#
#module load conda/4.12.0
#source ~/.bashrc
#conda activate bbduk
#
#for recipe_type in recipe1_standard recipe2_smallfrag combined_recipe1_recipe2; do
#    echo "BBDuk pour ${recipe_type}..."
#    INPUT_DIR="${BASE_DIR}/01_raw_data/${recipe_type}"
#    OUTPUT_DIR="${BASE_DIR}/03_bbduk/${recipe_type}"
#    mkdir -p "$OUTPUT_DIR"
#    
#    cd "$INPUT_DIR"
#    
#    for r1_file in *_R1.fastq.gz; do
#        r2_file="${r1_file/_R1/_R2}"
#        
#        if [[ ! -f "$r2_file" ]]; then
#            echo "ERREUR: Fichier R2 manquant pour $r1_file" >&2
#            continue
#        fi
#        
#        base_name="${r1_file%%_R1.fastq.gz}"
#        
#        echo "  → Traitement de ${base_name}..."
#        
#        $BBDUK -Xmx4g \
#            in1="$r1_file" \
#            in2="$r2_file" \
#            out1="${OUTPUT_DIR}/clean_${r1_file}" \
#            out2="${OUTPUT_DIR}/clean_${r2_file}" \
#            ref=$PHIX \
#            ktrim=rl \
#            k=23 \
#            mink=11 \
#            hdist=1 \
#            tpe \
#            tbo \
#            minlen=25 \
#            qtrim=r \
#            trimq=20 \
#            stats="${OUTPUT_DIR}/${base_name}_bbduk_stats.txt"
#    done
#done
#
#conda deactivate
#
#echo "Filtrage BBDuk terminé."

#################################################################################
## ÉTAPE 3: Déduplication avec FastUniq
#################################################################################
#
#echo ""
#echo "=== ÉTAPE 3: Déduplication (FastUniq) ==="
#echo ""
#
#module load conda/4.12.0
source ~/.bashrc
#conda activate fastuniq
#
#TMP="/tmp/fastuniq_corsica_tmp"
#mkdir -p "$TMP"
#
#for recipe_type in recipe1_standard recipe2_smallfrag combined_recipe1_recipe2; do
#    echo "FastUniq pour ${recipe_type}..."
#    INPUT_DIR="${BASE_DIR}/03_bbduk/${recipe_type}"
#    OUTPUT_DIR="${BASE_DIR}/04_fastuniq/${recipe_type}"
#    mkdir -p "$OUTPUT_DIR"
#    
#    cd "$INPUT_DIR" || continue
#    
#    for R1_gz in clean_*_R1.fastq.gz; do
#        base=$(echo "$R1_gz" | sed 's/_R1\.fastq\.gz//')
#        R2_gz="${base}_R2.fastq.gz"
#        
#        if [[ -f "$R2_gz" ]]; then
#            echo "  → Traitement de ${base}..."
#            
#            R1_tmp="${TMP}/${base}_R1.fastq"
#            R2_tmp="${TMP}/${base}_R2.fastq"
#            listfile="${TMP}/${base}.list"
#            
#            zcat "$INPUT_DIR/$R1_gz" > "$R1_tmp"
#            zcat "$INPUT_DIR/$R2_gz" > "$R2_tmp"
#            
#            echo -e "${R1_tmp}\n${R2_tmp}" > "$listfile"
#            
#            fastuniq -i "$listfile" -t q \
#                -o "${OUTPUT_DIR}/${base}_dedup_R1.fastq" \
#                -p "${OUTPUT_DIR}/${base}_dedup_R2.fastq"
#            
#            rm -f "$R1_tmp" "$R2_tmp" "$listfile"
#        else
#            echo "ATTENTION: fichier R2 manquant pour $base"
#        fi
#    done
#done
#
#rm -rf "$TMP"
#conda deactivate
#
#echo "Déduplication FastUniq terminée."

################################################################################
# ÉTAPE 4: Clumpify - Déduplication optique supplémentaire
################################################################################

echo ""
echo "=== ÉTAPE 4: Clumpify (déduplication optique) ==="
echo ""

module load conda/4.12.0
#source ~/.bashrc
conda activate bbduk

for recipe_type in recipe1_standard recipe2_smallfrag combined_recipe1_recipe2; do
    echo "Clumpify pour ${recipe_type}..."
    INPUT_DIR="${BASE_DIR}/04_fastuniq/${recipe_type}"
    OUTPUT_DIR="${BASE_DIR}/05_clumpify/${recipe_type}"
    mkdir -p "$OUTPUT_DIR"
    
    for R1 in "${INPUT_DIR}"/*_R1.fastq; do
        R2="${R1/_R1.fastq/_R2.fastq}"
        
        if [[ -f "$R2" ]]; then
            base=$(basename "$R1" _R1.fastq)
            echo "  → Traitement de ${base}..."
            
            $CLUMPIFY \
                in="$R1" in2="$R2" \
                out="${OUTPUT_DIR}/${base}_clumpify_R1.fastq.gz" \
                out2="${OUTPUT_DIR}/${base}_clumpify_R2.fastq.gz" \
                dedupe=t
        else
            echo "Fichier R2 manquant pour $R1, ignoré."
        fi
    done
done

conda deactivate

echo "Clumpify terminé."

################################################################################
# ÉTAPE 5: Fastp - Merging et contrôle qualité final
################################################################################

echo ""
echo "=== ÉTAPE 5: Fastp (merging et QC final) ==="
echo ""

module load conda/4.12.0
source ~/.bashrc
conda activate fastp

for recipe_type in recipe1_standard recipe2_smallfrag combined_recipe1_recipe2; do
    echo "Fastp pour ${recipe_type}..."
    INPUT_DIR="${BASE_DIR}/05_clumpify/${recipe_type}"
    OUTPUT_DIR="${BASE_DIR}/06_fastp/${recipe_type}"
    LOG_DIR="${BASE_DIR}/00_scripts/fastp_logs/${recipe_type}"
    mkdir -p "$OUTPUT_DIR"
    mkdir -p "$LOG_DIR"
    
    for R1 in "${INPUT_DIR}"/*_R1.fastq.gz; do
        R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
        BASENAME=$(basename "$R1" _R1.fastq.gz)
        
        echo "  → Traitement de ${BASENAME}..."
        
        OUT_R1="${OUTPUT_DIR}/${BASENAME}_fastp_R1.fastq.gz"
        OUT_R2="${OUTPUT_DIR}/${BASENAME}_fastp_R2.fastq.gz"
        MERGED="${OUTPUT_DIR}/${BASENAME}_fastp_merged.fastq.gz"
        HTML="${LOG_DIR}/${BASENAME}_fastp.html"
        JSON="${LOG_DIR}/${BASENAME}_fastp.json"
        
        fastp \
            --in1 "$R1" --in2 "$R2" \
            --out1 "$OUT_R1" --out2 "$OUT_R2" \
            --merged_out "$MERGED" \
            --length_required 30 \
            --cut_front --cut_tail \
            --cut_window_size 4 \
            --cut_mean_quality 10 \
            --n_base_limit 5 \
            --unqualified_percent_limit 40 \
            --complexity_threshold 30 \
            --qualified_quality_phred 15 \
            --low_complexity_filter \
            --trim_poly_x \
            --poly_x_min_len 10 \
            --merge --correction \
            --overlap_len_require 10 \
            --overlap_diff_limit 5 \
            --overlap_diff_percent_limit 20 \
            --html "$HTML" \
            --json "$JSON" \
            --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            --detect_adapter_for_pe \
            --thread 4
    done
done

conda deactivate

echo "Fastp terminé."

################################################################################
# ÉTAPE 6: Classification taxonomique avec Kraken2
################################################################################

echo ""
echo "=== ÉTAPE 6: Classification taxonomique (Kraken2) ==="
echo ""

module load conda/4.12.0
source ~/.bashrc
conda activate kraken2

for recipe_type in recipe1_standard recipe2_smallfrag combined_recipe1_recipe2; do
    echo "Kraken2 pour ${recipe_type}..."
    FASTP_DIR="${BASE_DIR}/06_fastp/${recipe_type}"
    OUT_DIR="${BASE_DIR}/07_kraken2/${recipe_type}"
    mkdir -p "$OUT_DIR"
    
    # Analyse des reads merged (single-end)
    echo "  → Analyse des reads merged..."
    for MERGED in "${FASTP_DIR}"/*_fastp_merged.fastq.gz; do
        if [[ -f "$MERGED" ]]; then
            SAMPLE=$(basename "$MERGED" _fastp_merged.fastq.gz)
            OUT_KRAKEN="${OUT_DIR}/${SAMPLE}_merged.kraken"
            OUT_REPORT="${OUT_DIR}/${SAMPLE}_merged.report"
            
            echo "    • ${SAMPLE} (merged)"
            
            kraken2 --confidence 0.2 --db "$KRAKEN2_DB" --threads $THREADS \
                --output "$OUT_KRAKEN" --report "$OUT_REPORT" "$MERGED"
        fi
    done
    
    # Analyse des reads unmerged (paired-end)
    echo "  → Analyse des reads unmerged..."
    for R1 in "${FASTP_DIR}"/*_fastp_R1.fastq.gz; do
        if [[ -f "$R1" ]]; then
            SAMPLE=$(basename "$R1" _fastp_R1.fastq.gz)
            R2="${FASTP_DIR}/${SAMPLE}_fastp_R2.fastq.gz"
            
            if [[ -f "$R2" ]]; then
                OUT_KRAKEN="${OUT_DIR}/${SAMPLE}_unmerged.kraken"
                OUT_REPORT="${OUT_DIR}/${SAMPLE}_unmerged.report"
                
                echo "    • ${SAMPLE} (unmerged)"
                
                kraken2 --confidence 0.2 --paired --db "$KRAKEN2_DB" --threads $THREADS \
                    --output "$OUT_KRAKEN" --report "$OUT_REPORT" "$R1" "$R2"
            fi
        fi
    done
done

conda deactivate

echo "Classification Kraken2 terminée."

################################################################################
# ÉTAPE 7: Visualisation avec Krona
################################################################################

echo ""
echo "=== ÉTAPE 7: Visualisation (Krona) ==="
echo ""

module load conda/4.12.0
source ~/.bashrc
conda activate krona

for recipe_type in recipe1_standard recipe2_smallfrag combined_recipe1_recipe2; do
    echo "Krona pour ${recipe_type}..."
    IN_DIR="${BASE_DIR}/07_kraken2/${recipe_type}"
    OUT_DIR="${BASE_DIR}/08_krona/${recipe_type}"
    mkdir -p "$OUT_DIR"
    
    cd "$IN_DIR"
    
    # Krona pour tous les échantillons combinés
    ktImportTaxonomy -t 5 -m 3 -o "${OUT_DIR}/all_samples_krona.html" "${IN_DIR}"/*.report
    
    # Krona individuel pour chaque échantillon
    for report in "${IN_DIR}"/*.report; do
        if [[ -f "$report" ]]; then
            base=$(basename "$report" .report)
            ktImportTaxonomy -t 5 -m 3 -o "${OUT_DIR}/${base}_krona.html" "$report"
        fi
    done
done

conda deactivate

echo "Visualisation Krona terminée."

################################################################################
# ÉTAPE 8: Création des tables d'assignation taxonomique
################################################################################

echo ""
echo "=== ÉTAPE 8: Création des tables MPA ==="
echo ""

for recipe_type in recipe1_standard recipe2_smallfrag combined_recipe1_recipe2; do
    echo "Conversion MPA pour ${recipe_type}..."
    IN_DIR="${BASE_DIR}/07_kraken2/${recipe_type}"
    OUT_DIR="${BASE_DIR}/09_mpa_tables/${recipe_type}"
    mkdir -p "$OUT_DIR"
    
    cd "$IN_DIR"
    
    # Conversion kreport vers format MPA
    declare -a mpa_files=()
    for report in *.report; do
        if [[ -f "$report" ]]; then
            base=$(basename "$report" .report)
            mpa_file="${OUT_DIR}/${base}.mpa"
            
            echo "  → Conversion de ${base}..."
            python3 "${KRAKENTOOLS_DIR}/kreport2mpa.py" -r "$report" -o "$mpa_file"
            
            mpa_files+=("$mpa_file")
        fi
    done
    
    # Combinaison de tous les fichiers MPA en une seule table
    if [[ ${#mpa_files[@]} -gt 0 ]]; then
        echo "  → Combinaison de tous les fichiers MPA..."
        python3 "${KRAKENTOOLS_DIR}/combine_mpa.py" -i "${mpa_files[@]}" -o "${OUT_DIR}/combined_all.tsv"
    fi
done

echo "Création des tables MPA terminée."

################################################################################
# ÉTAPE 9: Génération du rapport de comparaison
################################################################################

echo ""
echo "=== ÉTAPE 9: Génération du rapport de comparaison ==="
echo ""

REPORT_FILE="${BASE_DIR}/00_scripts/pipeline_comparison_report.txt"

cat > "$REPORT_FILE" << EOF
================================================================================
RAPPORT D'ANALYSE - CORSICA WRECK SEDADNA
================================================================================

Date d'analyse: $(date)
Pipeline version: 1.0

================================================================================
DESCRIPTION DU PROJET
================================================================================

Ce pipeline compare trois stratégies de séquençage pour 4 échantillons sedaDNA:

Échantillons:
  - 1121_sed8_rep1 (sédiment 8, réplicat 1)
  - 1122_sed8_rep2 (sédiment 8, réplicat 2)
  - 1129_sed6_rep1 (sédiment 6, réplicat 1)
  - 1130_sed6_rep2 (sédiment 6, réplicat 2)

Stratégies de séquençage:
  1. Recipe1 (recipe1_standard): recipes-ShortInsert (2x150bp)
     → Séquençage standard pour fragments courts
     
  2. Recipe2 (recipe2_smallfrag): recipes-v3.3.x-ShortInsert-SmallFragmentRecovery (2x150bp)
     → Optimisé pour la récupération de très petits fragments (aDNA)
     
  3. Combined (combined_recipe1_recipe2): Recipe1 + Recipe2
     → Fusion des deux datasets pour maximiser la profondeur de séquençage

================================================================================
ÉTAPES DU PIPELINE
================================================================================

1. Contrôle qualité (FastQC/MultiQC)
2. Filtrage et trimming (BBDuk) - Élimination des adaptateurs et phiX
3. Déduplication (FastUniq) - Suppression des duplicats PCR
4. Clumpify - Déduplication optique supplémentaire
5. Fastp - Merging des reads et contrôle qualité final
6. Classification taxonomique (Kraken2) - Base de données nt
7. Visualisation (Krona) - Charts interactifs de composition taxonomique
8. Tables d'assignation (MPA format) - Pour analyses downstream

================================================================================
LOCALISATION DES RÉSULTATS
================================================================================

Répertoire principal: ${BASE_DIR}

Structure des résultats:
  ├── 01_raw_data/                    # Données brutes
  │   ├── recipe1_standard/           # Recipe 1
  │   ├── recipe2_smallfrag/          # Recipe 2
  │   └── combined_recipe1_recipe2/   # Recipe 1+2
  │
  ├── 02_quality_check/               # Rapports FastQC/MultiQC
  ├── 03_bbduk/                       # Reads filtrés
  ├── 04_fastuniq/                    # Reads dédupliqués
  ├── 05_clumpify/                    # Déduplication optique
  ├── 06_fastp/                       # Reads merged et unmerged finaux
  │
  ├── 07_kraken2/                     # Classification taxonomique
  │   ├── recipe1_standard/
  │   ├── recipe2_smallfrag/
  │   └── combined_recipe1_recipe2/
  │
  ├── 08_krona/                       # Visualisations interactives
  │   ├── recipe1_standard/
  │   ├── recipe2_smallfrag/
  │   └── combined_recipe1_recipe2/
  │
  └── 09_mpa_tables/                  # Tables d'assignation
      ├── recipe1_standard/
      ├── recipe2_smallfrag/
      └── combined_recipe1_recipe2/

================================================================================
ANALYSES DE COMPARAISON RECOMMANDÉES
================================================================================

Pour comparer les performances des 3 stratégies, examinez:

1. Nombre de reads à chaque étape:
   - Comparez le taux de rétention entre Recipe1, Recipe2 et Combined
   - Vérifiez le taux de merging (indicateur de qualité aDNA)

2. Composition taxonomique:
   - Ouvrez les fichiers Krona pour visualisation interactive
   - Comparez la diversité détectée entre les recettes
   - Identifiez les taxons uniquement détectés par Recipe2 (fragments courts)

3. Analyse quantitative:
   - Utilisez les fichiers combined_all.tsv dans 09_mpa_tables/
   - Importez dans R pour analyses statistiques et visualisations

4. Évaluation de la valeur ajoutée:
   - Recipe2 vs Recipe1: Récupération de fragments courts supplémentaires?
   - Combined vs Recipe1/Recipe2: Augmentation significative de profondeur?

================================================================================
COMMANDES POUR ANALYSES COMPLÉMENTAIRES EN R
================================================================================

# Chargement et comparaison des tables MPA
library(tidyverse)

# Recipe 1
mpa_r1 <- read_tsv("${BASE_DIR}/09_mpa_tables/recipe1_standard/combined_all.tsv")

# Recipe 2
mpa_r2 <- read_tsv("${BASE_DIR}/09_mpa_tables/recipe2_smallfrag/combined_all.tsv")

# Combined
mpa_combined <- read_tsv("${BASE_DIR}/09_mpa_tables/combined_recipe1_recipe2/combined_all.tsv")

# Comparaison de la richesse taxonomique
richness_comparison <- tibble(
  Recipe = c("Recipe1_Standard", "Recipe2_SmallFrag", "Combined"),
  N_taxa = c(nrow(mpa_r1), nrow(mpa_r2), nrow(mpa_combined))
)

# Visualisation
ggplot(richness_comparison, aes(x = Recipe, y = N_taxa, fill = Recipe)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Richesse taxonomique par stratégie de séquençage",
       y = "Nombre de taxons détectés")

================================================================================
INFORMATIONS SUPPLÉMENTAIRES
================================================================================

Base de données Kraken2: ${KRAKEN2_DB}
KrakenTools: ${KRAKENTOOLS_DIR}

Pour toute question, contacter: pierrelouis.stenger@gmail.com

================================================================================
EOF

cat "$REPORT_FILE"

echo ""
echo "Rapport généré: $REPORT_FILE"

################################################################################
# ÉTAPE 10: Statistiques finales
################################################################################

echo ""
echo "=== ÉTAPE 10: Génération des statistiques ==="
echo ""

STATS_FILE="${BASE_DIR}/00_scripts/pipeline_statistics.txt"

cat > "$STATS_FILE" << EOF
================================================================================
STATISTIQUES D'ANALYSE - CORSICA WRECK
================================================================================

Date: $(date)

EOF

for recipe_type in recipe1_standard recipe2_smallfrag combined_recipe1_recipe2; do
    echo "Statistiques pour ${recipe_type}..." >> "$STATS_FILE"
    echo "----------------------------------------" >> "$STATS_FILE"
    
    # Compter les fichiers à différentes étapes
    RAW_COUNT=$(ls "${BASE_DIR}/01_raw_data/${recipe_type}"/*.fastq.gz 2>/dev/null | wc -l)
    BBDUK_COUNT=$(ls "${BASE_DIR}/03_bbduk/${recipe_type}"/clean_*.fastq.gz 2>/dev/null | wc -l)
    FASTP_COUNT=$(ls "${BASE_DIR}/06_fastp/${recipe_type}"/*_merged.fastq.gz 2>/dev/null | wc -l)
    KRAKEN_COUNT=$(ls "${BASE_DIR}/07_kraken2/${recipe_type}"/*.report 2>/dev/null | wc -l)
    
    echo "Fichiers bruts: ${RAW_COUNT}" >> "$STATS_FILE"
    echo "Fichiers après BBDuk: ${BBDUK_COUNT}" >> "$STATS_FILE"
    echo "Fichiers merged (Fastp): ${FASTP_COUNT}" >> "$STATS_FILE"
    echo "Rapports Kraken2: ${KRAKEN_COUNT}" >> "$STATS_FILE"
    echo "" >> "$STATS_FILE"
done

cat "$STATS_FILE"

################################################################################
# FIN DU PIPELINE
################################################################################

echo ""
echo "=========================================="
echo "PIPELINE TERMINÉ AVEC SUCCÈS"
echo "Date de fin: $(date)"
echo "=========================================="
echo ""
echo "Résultats disponibles dans: ${BASE_DIR}"
echo "Rapport: ${REPORT_FILE}"
echo "Statistiques: ${STATS_FILE}"
echo ""
echo "Prochaines étapes:"
echo "  1. Consultez les rapports MultiQC dans 02_quality_check/"
echo "  2. Explorez les visualisations Krona dans 08_krona/"
echo "  3. Analysez les tables MPA dans 09_mpa_tables/ avec R"
echo "  4. Comparez les performances des 3 stratégies de séquençage"
echo ""

# Envoi d'un email de notification (si configuré)
echo "Pipeline Corsica wreck terminé le $(date). Résultats: ${BASE_DIR}" | \
    mail -s "Pipeline sedaDNA Corsica - Terminé" pierrelouis.stenger@gmail.com 2>/dev/null || true

exit 0
