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
# CORRECTION: Gestion des variables d'environnement bash pour HPC
################################################################################

# CORRECTION: Configuration bash pour environnement non-interactif
set -euo pipefail
export PS1='$ '  # Définir PS1 pour éviter l'erreur "variable sans liaison"

# CORRECTION: Source bashrc de façon sécurisée
if [[ -f ~/.bashrc ]]; then
    set +u  # Temporairement désactiver 'unset variable' check
    source ~/.bashrc
    set -u  # Réactiver le check
fi

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
# FONCTION: Activation conda sécurisée
################################################################################

safe_conda_activate() {
    local env_name=$1
    echo "Activation de l'environnement conda: $env_name"
    
    # Chargement du module conda
    module load conda/4.12.0 2>/dev/null || echo "Module conda déjà chargé"
    
    # Source bashrc seulement si nécessaire
    if [[ ! $(type -t conda) == "function" ]]; then
        set +u
        source ~/.bashrc 2>/dev/null || true
        set -u
    fi
    
    # Activation de l'environnement
    conda activate "$env_name" 2>/dev/null || {
        echo "Erreur lors de l'activation de $env_name"
        conda info --envs
        exit 1
    }
    
    echo "Environnement $env_name activé avec succès"
}

safe_conda_deactivate() {
    conda deactivate 2>/dev/null || true
    echo "Environnement conda désactivé"
}

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

# Copie des fichiers Recipe 1 (Standard ShortInsert)
echo "Copie des fichiers Recipe1 (Standard ShortInsert)..."
for sample in "${SAMPLES[@]}"; do
    if [[ -f "${RECIPE1_DIR}/${sample}_R1.fastq.gz" ]]; then
        cp "${RECIPE1_DIR}/${sample}_R1.fastq.gz" "${BASE_DIR}/01_raw_data/recipe1_standard/"
        echo "  ✓ ${sample}_R1.fastq.gz copié"
    else
        echo "  ✗ ${sample}_R1.fastq.gz non trouvé"
    fi
    
    if [[ -f "${RECIPE1_DIR}/${sample}_R2.fastq.gz" ]]; then
        cp "${RECIPE1_DIR}/${sample}_R2.fastq.gz" "${BASE_DIR}/01_raw_data/recipe1_standard/"
        echo "  ✓ ${sample}_R2.fastq.gz copié"
    else
        echo "  ✗ ${sample}_R2.fastq.gz non trouvé"
    fi
done

# Copie des fichiers Recipe 2 (SmallFragmentRecovery)
echo "Copie des fichiers Recipe2 (SmallFragmentRecovery)..."
for sample in "${SAMPLES[@]}"; do
    if [[ -f "${RECIPE2_DIR}/${sample}/${sample}_R1.fastq.gz" ]]; then
        cp "${RECIPE2_DIR}/${sample}/${sample}_R1.fastq.gz" "${BASE_DIR}/01_raw_data/recipe2_smallfrag/"
        echo "  ✓ ${sample}_R1.fastq.gz copié"
    else
        echo "  ✗ ${sample}_R1.fastq.gz non trouvé dans ${RECIPE2_DIR}/${sample}/"
    fi
    
    if [[ -f "${RECIPE2_DIR}/${sample}/${sample}_R2.fastq.gz" ]]; then
        cp "${RECIPE2_DIR}/${sample}/${sample}_R2.fastq.gz" "${BASE_DIR}/01_raw_data/recipe2_smallfrag/"
        echo "  ✓ ${sample}_R2.fastq.gz copié"
    else
        echo "  ✗ ${sample}_R2.fastq.gz non trouvé dans ${RECIPE2_DIR}/${sample}/"
    fi
done

# Création des fichiers combinés (Recipe1 + Recipe2)
echo "Fusion des fichiers Recipe1 + Recipe2..."
for sample in "${SAMPLES[@]}"; do
    R1_RECIPE1="${BASE_DIR}/01_raw_data/recipe1_standard/${sample}_R1.fastq.gz"
    R2_RECIPE1="${BASE_DIR}/01_raw_data/recipe1_standard/${sample}_R2.fastq.gz"
    R1_RECIPE2="${BASE_DIR}/01_raw_data/recipe2_smallfrag/${sample}_R1.fastq.gz"
    R2_RECIPE2="${BASE_DIR}/01_raw_data/recipe2_smallfrag/${sample}_R2.fastq.gz"
    
    # Fusion R1
    if [[ -f "$R1_RECIPE1" ]] && [[ -f "$R1_RECIPE2" ]]; then
        cat "$R1_RECIPE1" "$R1_RECIPE2" > "${BASE_DIR}/01_raw_data/combined_recipe1_recipe2/${sample}_R1.fastq.gz"
        echo "  ✓ ${sample}_R1 combiné créé"
    else
        echo "  ✗ Impossible de créer ${sample}_R1 combiné"
    fi
    
    # Fusion R2
    if [[ -f "$R2_RECIPE1" ]] && [[ -f "$R2_RECIPE2" ]]; then
        cat "$R2_RECIPE1" "$R2_RECIPE2" > "${BASE_DIR}/01_raw_data/combined_recipe1_recipe2/${sample}_R2.fastq.gz"
        echo "  ✓ ${sample}_R2 combiné créé"
    else
        echo "  ✗ Impossible de créer ${sample}_R2 combiné"
    fi
done

echo "Préparation des données terminée."

################################################################################
# ÉTAPE 1: Contrôle qualité avec FastQC et MultiQC
################################################################################

echo ""
echo "=== ÉTAPE 1: Contrôle qualité (FastQC/MultiQC) ==="
echo ""

safe_conda_activate fastqc

# FastQC pour chaque type de recette
for recipe_type in recipe1_standard recipe2_smallfrag combined_recipe1_recipe2; do
    echo "FastQC pour ${recipe_type}..."
    RECIPE_DIR="${BASE_DIR}/01_raw_data/${recipe_type}"
    OUTPUT_DIR="${BASE_DIR}/02_quality_check/${recipe_type}"
    mkdir -p "$OUTPUT_DIR"
    
    for FILE in "${RECIPE_DIR}"/*.fastq.gz; do
        if [[ -f "$FILE" ]]; then
            fastqc "$FILE" -o "$OUTPUT_DIR" -t 4 --quiet
            echo "  ✓ $(basename "$FILE") analysé"
        fi
    done
done

safe_conda_deactivate

# MultiQC pour résumer tous les rapports
safe_conda_activate multiqc
for recipe_type in recipe1_standard recipe2_smallfrag combined_recipe1_recipe2; do
    echo "MultiQC pour ${recipe_type}..."
    OUTPUT_DIR="${BASE_DIR}/02_quality_check/${recipe_type}"
    cd "$OUTPUT_DIR" || continue
    multiqc . -n "multiqc_${recipe_type}.html" --quiet
    echo "  ✓ Rapport MultiQC généré: multiqc_${recipe_type}.html"
done
safe_conda_deactivate

echo "Contrôle qualité terminé."

################################################################################
# ÉTAPE 2: Filtrage et trimming avec BBDuk
################################################################################

echo ""
echo "=== ÉTAPE 2: Filtrage et trimming (BBDuk) ==="
echo ""

# BBDuk ne nécessite pas conda, utilisation directe
for recipe_type in recipe1_standard recipe2_smallfrag combined_recipe1_recipe2; do
    echo "BBDuk pour ${recipe_type}..."
    INPUT_DIR="${BASE_DIR}/01_raw_data/${recipe_type}"
    OUTPUT_DIR="${BASE_DIR}/03_bbduk/${recipe_type}"
    mkdir -p "$OUTPUT_DIR"
    
    cd "$INPUT_DIR" || continue
    
    for r1_file in *_R1.fastq.gz; do
        if [[ ! -f "$r1_file" ]]; then continue; fi
        
        r2_file="${r1_file/_R1/_R2}"
        
        if [[ ! -f "$r2_file" ]]; then
            echo "  ✗ Fichier R2 manquant pour $r1_file"
            continue
        fi
        
        base_name="${r1_file%%_R1.fastq.gz}"
        
        echo "  → Traitement de ${base_name}..."
        
        $BBDUK -Xmx4g \
            in1="$r1_file" \
            in2="$r2_file" \
            out1="${OUTPUT_DIR}/clean_${r1_file}" \
            out2="${OUTPUT_DIR}/clean_${r2_file}" \
            ref=$PHIX \
            ktrim=rl \
            k=23 \
            mink=11 \
            hdist=1 \
            tpe \
            tbo \
            minlen=25 \
            qtrim=r \
            trimq=20 \
            stats="${OUTPUT_DIR}/${base_name}_bbduk_stats.txt" 2>/dev/null
        
        echo "  ✓ ${base_name} traité avec succès"
    done
done

echo "Filtrage BBDuk terminé."

################################################################################
# ÉTAPE 3: Déduplication avec FastUniq
################################################################################

echo ""
echo "=== ÉTAPE 3: Déduplication (FastUniq) ==="
echo ""

safe_conda_activate fastuniq

TMP="${BASE_DIR}/tmp_fastuniq"
mkdir -p "$TMP"

for recipe_type in recipe1_standard recipe2_smallfrag combined_recipe1_recipe2; do
    echo "FastUniq pour ${recipe_type}..."
    INPUT_DIR="${BASE_DIR}/03_bbduk/${recipe_type}"
    OUTPUT_DIR="${BASE_DIR}/04_fastuniq/${recipe_type}"
    mkdir -p "$OUTPUT_DIR"
    
    cd "$INPUT_DIR" || continue
    
    for R1_gz in clean_*_R1.fastq.gz; do
        if [[ ! -f "$R1_gz" ]]; then continue; fi
        
        base=$(echo "$R1_gz" | sed 's/clean_//' | sed 's/_R1\.fastq\.gz//')
        R2_gz="clean_${base}_R2.fastq.gz"
        
        if [[ -f "$R2_gz" ]]; then
            echo "  → Traitement de ${base}..."
            
            R1_tmp="${TMP}/${base}_R1.fastq"
            R2_tmp="${TMP}/${base}_R2.fastq"
            listfile="${TMP}/${base}.list"
            
            zcat "$R1_gz" > "$R1_tmp"
            zcat "$R2_gz" > "$R2_tmp"
            
            echo -e "${R1_tmp}\n${R2_tmp}" > "$listfile"
            
            fastuniq -i "$listfile" -t q \
                -o "${OUTPUT_DIR}/${base}_dedup_R1.fastq" \
                -p "${OUTPUT_DIR}/${base}_dedup_R2.fastq" 2>/dev/null
            
            rm -f "$R1_tmp" "$R2_tmp" "$listfile"
            echo "  ✓ ${base} dédupliqué"
        fi
    done
done

rm -rf "$TMP"
safe_conda_deactivate

echo "Déduplication FastUniq terminée."

################################################################################
# ÉTAPE 4: Clumpify - Déduplication optique supplémentaire
################################################################################

echo ""
echo "=== ÉTAPE 4: Clumpify (déduplication optique) ==="
echo ""

for recipe_type in recipe1_standard recipe2_smallfrag combined_recipe1_recipe2; do
    echo "Clumpify pour ${recipe_type}..."
    INPUT_DIR="${BASE_DIR}/04_fastuniq/${recipe_type}"
    OUTPUT_DIR="${BASE_DIR}/05_clumpify/${recipe_type}"
    mkdir -p "$OUTPUT_DIR"
    
    for R1 in "${INPUT_DIR}"/*_dedup_R1.fastq; do
        if [[ ! -f "$R1" ]]; then continue; fi
        
        R2="${R1/_R1.fastq/_R2.fastq}"
        
        if [[ -f "$R2" ]]; then
            base=$(basename "$R1" _dedup_R1.fastq)
            echo "  → Traitement de ${base}..."
            
            $CLUMPIFY \
                in="$R1" in2="$R2" \
                out="${OUTPUT_DIR}/${base}_clumpify_R1.fastq.gz" \
                out2="${OUTPUT_DIR}/${base}_clumpify_R2.fastq.gz" \
                dedupe=t 2>/dev/null
            
            echo "  ✓ ${base} clumpifié"
        fi
    done
done

echo "Clumpify terminé."

################################################################################
# ÉTAPE 5: Fastp - Merging et contrôle qualité final
################################################################################

echo ""
echo "=== ÉTAPE 5: Fastp (merging et QC final) ==="
echo ""

safe_conda_activate fastp

for recipe_type in recipe1_standard recipe2_smallfrag combined_recipe1_recipe2; do
    echo "Fastp pour ${recipe_type}..."
    INPUT_DIR="${BASE_DIR}/05_clumpify/${recipe_type}"
    OUTPUT_DIR="${BASE_DIR}/06_fastp/${recipe_type}"
    LOG_DIR="${BASE_DIR}/00_scripts/fastp_logs/${recipe_type}"
    mkdir -p "$OUTPUT_DIR" "$LOG_DIR"
    
    for R1 in "${INPUT_DIR}"/*_clumpify_R1.fastq.gz; do
        if [[ ! -f "$R1" ]]; then continue; fi
        
        R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
        BASENAME=$(basename "$R1" _clumpify_R1.fastq.gz)
        
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
            --thread 4 2>/dev/null
        
        echo "  ✓ ${BASENAME} traité par Fastp"
    done
done

safe_conda_deactivate

echo "Fastp terminé."

################################################################################
# ÉTAPE 6: Classification taxonomique avec Kraken2
################################################################################

echo ""
echo "=== ÉTAPE 6: Classification taxonomique (Kraken2) ==="
echo ""

safe_conda_activate kraken2

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
                --output "$OUT_KRAKEN" --report "$OUT_REPORT" "$MERGED" 2>/dev/null
            
            echo "    ✓ ${SAMPLE} merged analysé"
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
                    --output "$OUT_KRAKEN" --report "$OUT_REPORT" "$R1" "$R2" 2>/dev/null
                
                echo "    ✓ ${SAMPLE} unmerged analysé"
            fi
        fi
    done
done

safe_conda_deactivate

echo "Classification Kraken2 terminée."

################################################################################
# ÉTAPE 7: Visualisation avec Krona
################################################################################

echo ""
echo "=== ÉTAPE 7: Visualisation (Krona) ==="
echo ""

safe_conda_activate krona

for recipe_type in recipe1_standard recipe2_smallfrag combined_recipe1_recipe2; do
    echo "Krona pour ${recipe_type}..."
    IN_DIR="${BASE_DIR}/07_kraken2/${recipe_type}"
    OUT_DIR="${BASE_DIR}/08_krona/${recipe_type}"
    mkdir -p "$OUT_DIR"
    
    cd "$IN_DIR" || continue
    
    # Vérifier qu'il y a des reports
    if ls *.report 1> /dev/null 2>&1; then
        # Krona pour tous les échantillons combinés
        ktImportTaxonomy -t 5 -m 3 -o "${OUT_DIR}/all_samples_krona.html" *.report 2>/dev/null
        echo "  ✓ Krona combiné généré: all_samples_krona.html"
        
        # Krona individuel pour chaque échantillon
        for report in *.report; do
            if [[ -f "$report" ]]; then
                base=$(basename "$report" .report)
                ktImportTaxonomy -t 5 -m 3 -o "${OUT_DIR}/${base}_krona.html" "$report" 2>/dev/null
                echo "  ✓ Krona individuel: ${base}_krona.html"
            fi
        done
    else
        echo "  ✗ Aucun fichier report trouvé dans $IN_DIR"
    fi
done

safe_conda_deactivate

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
    
    cd "$IN_DIR" || continue
    
    # Conversion kreport vers format MPA
    declare -a mpa_files=()
    for report in *.report; do
        if [[ -f "$report" ]]; then
            base=$(basename "$report" .report)
            mpa_file="${OUT_DIR}/${base}.mpa"
            
            echo "  → Conversion de ${base}..."
            python3 "${KRAKENTOOLS_DIR}/kreport2mpa.py" -r "$report" -o "$mpa_file" 2>/dev/null
            
            if [[ -f "$mpa_file" ]]; then
                mpa_files+=("$mpa_file")
                echo "  ✓ ${base}.mpa créé"
            fi
        fi
    done
    
    # Combinaison de tous les fichiers MPA en une seule table
    if [[ ${#mpa_files[@]} -gt 0 ]]; then
        echo "  → Combinaison de ${#mpa_files[@]} fichiers MPA..."
        python3 "${KRAKENTOOLS_DIR}/combine_mpa.py" -i "${mpa_files[@]}" -o "${OUT_DIR}/combined_all.tsv" 2>/dev/null
        echo "  ✓ Table combinée: combined_all.tsv"
    fi
done

echo "Création des tables MPA terminée."

################################################################################
# ÉTAPE 9: Génération du rapport final et statistiques
################################################################################

echo ""
echo "=== ÉTAPE 9: Génération du rapport final ==="
echo ""

REPORT_FILE="${BASE_DIR}/00_scripts/pipeline_report_$(date +%Y%m%d_%H%M%S).txt"

cat > "$REPORT_FILE" << EOF
================================================================================
RAPPORT D'ANALYSE PIPELINE CORSICA WRECK - $(date)
================================================================================

PIPELINE TERMINÉ AVEC SUCCÈS

Échantillons analysés: ${SAMPLES[@]}

Stratégies de séquençage:
  1. recipe1_standard      - recipes-ShortInsert standard
  2. recipe2_smallfrag     - recipes-v3.3.x-ShortInsert-SmallFragmentRecovery  
  3. combined_recipe1_recipe2 - Fusion des deux stratégies

================================================================================
LOCALISATION DES RÉSULTATS
================================================================================

Répertoire principal: ${BASE_DIR}

Résultats par étape:
  • Contrôle qualité:     02_quality_check/
  • Reads filtrés:        03_bbduk/ → 04_fastuniq/ → 05_clumpify/
  • Reads finaux:         06_fastp/
  • Classification:       07_kraken2/
  • Visualisations:       08_krona/
  • Tables d'analyse:     09_mpa_tables/

================================================================================
ANALYSES RECOMMANDÉES
================================================================================

1. Consulter les rapports MultiQC:
   firefox ${BASE_DIR}/02_quality_check/*/multiqc_*.html

2. Explorer les visualisations Krona:
   firefox ${BASE_DIR}/08_krona/*/all_samples_krona.html

3. Analyser avec R:
   library(tidyverse)
   r1 <- read_tsv("${BASE_DIR}/09_mpa_tables/recipe1_standard/combined_all.tsv")
   r2 <- read_tsv("${BASE_DIR}/09_mpa_tables/recipe2_smallfrag/combined_all.tsv")
   combined <- read_tsv("${BASE_DIR}/09_mpa_tables/combined_recipe1_recipe2/combined_all.tsv")

================================================================================
STATISTIQUES FINALES
================================================================================

EOF

# Ajout des statistiques
for recipe_type in recipe1_standard recipe2_smallfrag combined_recipe1_recipe2; do
    echo "Statistiques ${recipe_type}:" >> "$REPORT_FILE"
    
    RAW=$(find "${BASE_DIR}/01_raw_data/${recipe_type}" -name "*.fastq.gz" 2>/dev/null | wc -l)
    MERGED=$(find "${BASE_DIR}/06_fastp/${recipe_type}" -name "*_merged.fastq.gz" 2>/dev/null | wc -l)
    REPORTS=$(find "${BASE_DIR}/07_kraken2/${recipe_type}" -name "*.report" 2>/dev/null | wc -l)
    
    echo "  - Fichiers bruts: ${RAW}" >> "$REPORT_FILE"
    echo "  - Fichiers merged: ${MERGED}" >> "$REPORT_FILE"
    echo "  - Rapports Kraken2: ${REPORTS}" >> "$REPORT_FILE"
    echo "" >> "$REPORT_FILE"
done

cat >> "$REPORT_FILE" << EOF
================================================================================
PIPELINE TERMINE - $(date)
================================================================================
EOF

cat "$REPORT_FILE"

################################################################################
# FIN DU PIPELINE
################################################################################

echo ""
echo "=========================================="
echo "PIPELINE TERMINÉ AVEC SUCCÈS"
echo "Date de fin: $(date)"
echo "=========================================="
echo ""
echo "Rapport final: $REPORT_FILE"
echo ""
echo "Pour analyser les résultats:"
echo "  cd ${BASE_DIR}"
echo "  ls -lah */recipe*/"
echo ""

# Notification email (si configuré)
echo "Pipeline Corsica wreck terminé - $(date)" | \
    mail -s "Pipeline sedaDNA - Terminé" pierrelouis.stenger@gmail.com 2>/dev/null || true

exit 0
