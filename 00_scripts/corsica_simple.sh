#!/bin/bash

#SBATCH --job-name=corsica_wreck
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=/home/plstenge/seda_DNA_Corsican_wreck_comparison/corsica_wreck_%j.out
#SBATCH --error=/home/plstenge/seda_DNA_Corsican_wreck_comparison/corsica_wreck_%j.err

################################################################################
# Pipeline sedaDNA - Corsica wreck
# Comparaison: Recipe1 vs Recipe2 vs Combined
################################################################################

set -euo pipefail

# ===== CONFIGURATION =====
BASE_DIR="/home/plstenge/seda_DNA_Corsican_wreck_comparison"
RECIPE1_DIR="/home/plstenge/seda_DNA_Corsican_wreck/01_raw_data"
RECIPE2_DIR="/storage/groups/gdec/shared_paleo/E1531_final/run3_20251008_AV241601_E1531_Ps7_Ps8"

BBDUK="/home/plstenge/bbmap/bbduk.sh"
CLUMPIFY="/home/plstenge/bbmap/clumpify.sh"
PHIX="/home/plstenge/bbmap/resources/phix174_ill.ref.fa.gz"
KRAKEN2_DB="/home/plstenge/k2_core_nt_20250609"
KRAKENTOOLS="/home/plstenge/coprolites/07_kraken2/KrakenTools"
THREADS=36

# Créer BASE_DIR en premier (sinon SLURM ne peut pas écrire les logs)
mkdir -p "$BASE_DIR"

echo "=========================================="
echo "Pipeline sedaDNA - Corsica wreck"
echo "Démarrage: $(date)"
echo "=========================================="

# ===== ÉTAPE 0: Préparation =====
echo ""
echo "=== ÉTAPE 0: Préparation des données ==="

mkdir -p "$BASE_DIR"/{01_raw_data/{recipe1_standard,recipe2_smallfrag,combined_recipe1_recipe2},02_qc,03_bbduk,04_fastuniq,05_clumpify,06_fastp,07_kraken2,08_krona,09_mpa}

SAMPLES=("1121_sed8_rep1" "1122_sed8_rep2" "1129_sed6_rep1" "1130_sed6_rep2")

echo "Copie Recipe1..."
for s in "${SAMPLES[@]}"; do
    [[ -f "$RECIPE1_DIR/${s}_R1.fastq.gz" ]] && cp "$RECIPE1_DIR/${s}_R1.fastq.gz" "$BASE_DIR/01_raw_data/recipe1_standard/"
    [[ -f "$RECIPE1_DIR/${s}_R2.fastq.gz" ]] && cp "$RECIPE1_DIR/${s}_R2.fastq.gz" "$BASE_DIR/01_raw_data/recipe1_standard/"
done

echo "Copie Recipe2..."
for s in "${SAMPLES[@]}"; do
    [[ -f "$RECIPE2_DIR/$s/${s}_R1.fastq.gz" ]] && cp "$RECIPE2_DIR/$s/${s}_R1.fastq.gz" "$BASE_DIR/01_raw_data/recipe2_smallfrag/"
    [[ -f "$RECIPE2_DIR/$s/${s}_R2.fastq.gz" ]] && cp "$RECIPE2_DIR/$s/${s}_R2.fastq.gz" "$BASE_DIR/01_raw_data/recipe2_smallfrag/"
done

echo "Création Combined (R1+R2)..."
for s in "${SAMPLES[@]}"; do
    [[ -f "$BASE_DIR/01_raw_data/recipe1_standard/${s}_R1.fastq.gz" ]] && \
    [[ -f "$BASE_DIR/01_raw_data/recipe2_smallfrag/${s}_R1.fastq.gz" ]] && \
        cat "$BASE_DIR/01_raw_data/recipe1_standard/${s}_R1.fastq.gz" "$BASE_DIR/01_raw_data/recipe2_smallfrag/${s}_R1.fastq.gz" \
        > "$BASE_DIR/01_raw_data/combined_recipe1_recipe2/${s}_R1.fastq.gz"
    
    [[ -f "$BASE_DIR/01_raw_data/recipe1_standard/${s}_R2.fastq.gz" ]] && \
    [[ -f "$BASE_DIR/01_raw_data/recipe2_smallfrag/${s}_R2.fastq.gz" ]] && \
        cat "$BASE_DIR/01_raw_data/recipe1_standard/${s}_R2.fastq.gz" "$BASE_DIR/01_raw_data/recipe2_smallfrag/${s}_R2.fastq.gz" \
        > "$BASE_DIR/01_raw_data/combined_recipe1_recipe2/${s}_R2.fastq.gz"
done

# ===== ÉTAPE 1: FastQC =====
echo ""
echo "=== ÉTAPE 1: FastQC ==="

module load conda/4.12.0
source ~/.bashrc 2>/dev/null || true

conda activate fastqc

for recipe in recipe1_standard recipe2_smallfrag combined_recipe1_recipe2; do
    echo "FastQC - $recipe..."
    mkdir -p "$BASE_DIR/02_qc/$recipe"
    for f in "$BASE_DIR/01_raw_data/$recipe"/*.fastq.gz; do
        [[ -f "$f" ]] && fastqc "$f" -o "$BASE_DIR/02_qc/$recipe" -q
    done
done

conda deactivate

# ===== ÉTAPE 2: BBDuk =====
echo ""
echo "=== ÉTAPE 2: BBDuk (filtrage) ==="

for recipe in recipe1_standard recipe2_smallfrag combined_recipe1_recipe2; do
    echo "BBDuk - $recipe..."
    cd "$BASE_DIR/01_raw_data/$recipe"
    
    for r1 in *_R1.fastq.gz; do
        [[ ! -f "$r1" ]] && continue
        r2="${r1/_R1/_R2}"
        [[ ! -f "$r2" ]] && continue
        
        out_dir="$BASE_DIR/03_bbduk/$recipe"
        mkdir -p "$out_dir"
        
        $BBDUK -Xmx4g in1="$r1" in2="$r2" \
            out1="$out_dir/c_$r1" out2="$out_dir/c_$r2" \
            ref=$PHIX ktrim=rl k=23 mink=11 hdist=1 tpe tbo minlen=25 qtrim=r trimq=20 \
            stats="$out_dir/${r1%.fastq.gz}.stats" 2>/dev/null
    done
done

# ===== ÉTAPE 3: FastUniq =====
echo ""
echo "=== ÉTAPE 3: FastUniq (déduplication) ==="

conda activate fastuniq

for recipe in recipe1_standard recipe2_smallfrag combined_recipe1_recipe2; do
    echo "FastUniq - $recipe..."
    in_dir="$BASE_DIR/03_bbduk/$recipe"
    out_dir="$BASE_DIR/04_fastuniq/$recipe"
    mkdir -p "$out_dir"
    tmp="/tmp/fu_$recipe"
    mkdir -p "$tmp"
    
    cd "$in_dir"
    for r1 in c_*_R1.fastq.gz; do
        [[ ! -f "$r1" ]] && continue
        r2="c_${r1#c_}_R2.fastq.gz"
        [[ ! -f "$r2" ]] && continue
        
        base="${r1#c_}"
        base="${base%_R1.fastq.gz}"
        
        zcat "$r1" > "$tmp/${base}_R1.fastq"
        zcat "$r2" > "$tmp/${base}_R2.fastq"
        echo -e "$tmp/${base}_R1.fastq\n$tmp/${base}_R2.fastq" > "$tmp/${base}.list"
        
        fastuniq -i "$tmp/${base}.list" -t q -o "$out_dir/${base}_dedup_R1.fastq" -p "$out_dir/${base}_dedup_R2.fastq" 2>/dev/null
        rm -f "$tmp/${base}_"*.fastq "$tmp/${base}.list"
    done
    rm -rf "$tmp"
done

conda deactivate

# ===== ÉTAPE 4: Clumpify =====
echo ""
echo "=== ÉTAPE 4: Clumpify (dédup optique) ==="

module load conda/4.12.0
source ~/.bashrc 2>/dev/null || true
conda activate bbduk

for recipe in recipe1_standard recipe2_smallfrag combined_recipe1_recipe2; do
    echo "Clumpify - $recipe..."
    in_dir="$BASE_DIR/04_fastuniq/$recipe"
    out_dir="$BASE_DIR/05_clumpify/$recipe"
    mkdir -p "$out_dir"
    
    for r1 in "$in_dir"/*_R1.fastq; do
        [[ ! -f "$r1" ]] && continue
        r2="${r1/_R1/_R2}"
        [[ ! -f "$r2" ]] && continue
        
        base=$(basename "$r1" _dedup_R1.fastq)
        $CLUMPIFY in="$r1" in2="$r2" out="$out_dir/${base}_clump_R1.fastq.gz" out2="$out_dir/${base}_clump_R2.fastq.gz" dedupe=t 2>/dev/null
    done
done

conda deactivate

# ===== ÉTAPE 5: Fastp =====
echo ""
echo "=== ÉTAPE 5: Fastp (merging + QC final) ==="

module load conda/4.12.0
source ~/.bashrc 2>/dev/null || true
conda activate fastp

for recipe in recipe1_standard recipe2_smallfrag combined_recipe1_recipe2; do
    echo "Fastp - $recipe..."
    in_dir="$BASE_DIR/05_clumpify/$recipe"
    out_dir="$BASE_DIR/06_fastp/$recipe"
    mkdir -p "$out_dir"
    
    for r1 in "$in_dir"/*_R1.fastq.gz; do
        [[ ! -f "$r1" ]] && continue
        r2="${r1/_R1.fastq.gz/_R2.fastq.gz}"
        [[ ! -f "$r2" ]] && continue
        
        base=$(basename "$r1" _clump_R1.fastq.gz)
        fastp --in1 "$r1" --in2 "$r2" --out1 "$out_dir/${base}_1.fq.gz" --out2 "$out_dir/${base}_2.fq.gz" \
            --merged_out "$out_dir/${base}_merged.fq.gz" \
            --length_required 30 --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 10 \
            --n_base_limit 5 --unqualified_percent_limit 40 --complexity_threshold 30 \
            --qualified_quality_phred 15 --low_complexity_filter --trim_poly_x --poly_x_min_len 10 \
            --merge --correction --overlap_len_require 10 --overlap_diff_limit 5 --overlap_diff_percent_limit 20 \
            --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            --detect_adapter_for_pe --thread 4 --html "$out_dir/${base}_fastp.html" --json "$out_dir/${base}_fastp.json" 2>/dev/null
    done
done

conda deactivate

# ===== ÉTAPE 6: Kraken2 =====
echo ""
echo "=== ÉTAPE 6: Kraken2 (classification) ==="

module load conda/4.12.0
source ~/.bashrc 2>/dev/null || true
conda activate kraken2

for recipe in recipe1_standard recipe2_smallfrag combined_recipe1_recipe2; do
    echo "Kraken2 - $recipe..."
    in_dir="$BASE_DIR/06_fastp/$recipe"
    out_dir="$BASE_DIR/07_kraken2/$recipe"
    mkdir -p "$out_dir"
    
    # Merged
    for merged in "$in_dir"/*_merged.fq.gz; do
        [[ ! -f "$merged" ]] && continue
        base=$(basename "$merged" _merged.fq.gz)
        kraken2 --confidence 0.2 --db "$KRAKEN2_DB" --threads $THREADS \
            --output "$out_dir/${base}_m.kraken" --report "$out_dir/${base}_m.report" "$merged" 2>/dev/null
    done
    
    # Unmerged
    for r1 in "$in_dir"/*_1.fq.gz; do
        [[ ! -f "$r1" ]] && continue
        r2="${r1/_1.fq.gz/_2.fq.gz}"
        [[ ! -f "$r2" ]] && continue
        
        base=$(basename "$r1" _1.fq.gz)
        kraken2 --confidence 0.2 --paired --db "$KRAKEN2_DB" --threads $THREADS \
            --output "$out_dir/${base}_u.kraken" --report "$out_dir/${base}_u.report" "$r1" "$r2" 2>/dev/null
    done
done

conda deactivate

# ===== ÉTAPE 7: Krona =====
echo ""
echo "=== ÉTAPE 7: Krona (visualisation) ==="

module load conda/4.12.0
source ~/.bashrc 2>/dev/null || true
conda activate krona

for recipe in recipe1_standard recipe2_smallfrag combined_recipe1_recipe2; do
    echo "Krona - $recipe..."
    in_dir="$BASE_DIR/07_kraken2/$recipe"
    out_dir="$BASE_DIR/08_krona/$recipe"
    mkdir -p "$out_dir"
    
    cd "$in_dir"
    [[ $(ls *.report 2>/dev/null | wc -l) -gt 0 ]] && \
        ktImportTaxonomy -t 5 -m 3 -o "$out_dir/all.html" *.report 2>/dev/null
done

conda deactivate

# ===== ÉTAPE 8: MPA Tables =====
echo ""
echo "=== ÉTAPE 8: MPA Tables ==="

for recipe in recipe1_standard recipe2_smallfrag combined_recipe1_recipe2; do
    echo "MPA - $recipe..."
    in_dir="$BASE_DIR/07_kraken2/$recipe"
    out_dir="$BASE_DIR/09_mpa/$recipe"
    mkdir -p "$out_dir"
    
    cd "$in_dir"
    for report in *.report; do
        [[ -f "$report" ]] && python3 "$KRAKENTOOLS/kreport2mpa.py" -r "$report" -o "$out_dir/${report%.report}.mpa" 2>/dev/null
    done
    
    mpa_files=("$out_dir"/*.mpa)
    if [[ -f "${mpa_files[0]}" ]]; then
        python3 "$KRAKENTOOLS/combine_mpa.py" -i "${mpa_files[@]}" -o "$out_dir/combined.tsv" 2>/dev/null
    fi
done

echo ""
echo "=========================================="
echo "PIPELINE TERMINÉ"
echo "Fin: $(date)"
echo "Résultats: $BASE_DIR"
echo "=========================================="

exit 0
