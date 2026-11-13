#!/usr/bin/env bash
#SBATCH --job-name=triticum_k35_PCA
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --partition=smp
#SBATCH --mem=120G
#SBATCH --output=/home/plstenge/seda_DNA_Corsican_wreck_comparison/00_scripts/13_triticum_PCA_%j.out
#SBATCH --error=/home/plstenge/seda_DNA_Corsican_wreck_comparison/00_scripts/13_triticum_PCA_%j.err
# Only k35, merged shotgun FASTQs; no ITS selection and no k29/k25 contribution.
# Modes: preflight (default; safe inspection), run (explicitly enabled after review).
set -Eeo pipefail
module load conda/4.12.0
source ~/.bashrc
conda activate "${CONDA_ENV:-mapdamage_py39}"
set -Eeuo pipefail
export LC_ALL=C
trap 'printf "ERROR line %s: %s\n" "$LINENO" "$BASH_COMMAND" >&2' ERR
MODE=${MODE:-preflight}
PROJECT=${PROJECT:-/home/plstenge/seda_DNA_Corsican_wreck_comparison}
PIPELINE_DIR=${PIPELINE_DIR:-/home/plstenge/adna_analysis_new-main}
INPUT_RUN=${INPUT_RUN:-$PROJECT/11_triticum_ITS_cascade/runs/run_20261112_115451_9952692}
REF=${REF:-/storage/groups/gdec/shared_paleo/genomes_REF/160802_Svevo_v2_pseudomolecules.fasta}
VCF=${VCF:-/storage/groups/gdec/shared_paleo/VCFWhealbi/WhealbiV1_4x_all_chr.vcf.gz}
POP=${POP:-/storage/groups/gdec/shared_paleo/VCFWhealbi/Whealbiv1_popOrigin.txt}
BED=${BED:-/storage/groups/gdec/shared_paleo/VCFWhealbi/WhealbiV1_4x_target_positions.bed}
ADNA_GIT=${ADNA_GIT:-$PIPELINE_DIR/adna_variousscripts-master}
THREADS=${SLURM_CPUS_PER_TASK:-36}
MIN_MAPQ=${MIN_MAPQ:-20}
MIN_BASEQ=${MIN_BASEQ:-20}
OUT="$PROJECT/13_triticum_PCA_k35/runs/run_$(date +%Y%m%d_%H%M%S)_${SLURM_JOB_ID:-manual_$$}"
mkdir -p "$OUT"/{qc,reference,mapping,pileup,PCA}
exec > >(tee -a "$OUT/pipeline.log") 2>&1
echo "Output: $OUT; mode=$MODE"
cp -- "$0" "$OUT/qc/submitted_script.sh"
conda list > "$OUT/qc/conda_packages.txt"
# Locate exact helper basenames; do not silently pick between differing copies.
resolve_helper() {
 local name=$1 override=$2 fallback=$3
 if [[ -n "$override" ]]; then printf '%s\n' "$override"; return; fi
 local paths=()
 mapfile -t paths < <(find "$PIPELINE_DIR" -type f -name "$name" | sort)
 if [[ ${#paths[@]} == 1 ]]; then printf '%s\n' "${paths[0]}"
 elif [[ ${#paths[@]} == 0 && -r "$fallback" ]]; then printf '%s\n' "$fallback"
 else printf 'Cannot uniquely locate %s; set PCA_SCRIPT or SUMMARY_SCRIPT explicitly. Candidates:\n%s\n' "$name" "${paths[*]}" >&2; return 1
 fi
}
PCA_SCRIPT=$(resolve_helper annotateVCFFromMPileup_andRunPCA.sh "${PCA_SCRIPT:-}" /home/ayouncha/aDNA_analysis/Pipeline/Snakemake/PCA_script/annotateVCFFromMPileup_andRunPCA.sh)
SUMMARY_SCRIPT=$(resolve_helper SummerizePileup.py "${SUMMARY_SCRIPT:-}" /home/ayouncha/aDNA_analysis/Pipeline/Snakemake/PCA_script/SummerizePileup.py)
for path in "$REF" "$VCF" "$POP" "$BED" "$PCA_SCRIPT" "$SUMMARY_SCRIPT"; do
 [[ -s "$path" && -r "$path" ]] || { echo "Missing or unreadable: $path"; exit 1; }
done
[[ -d "$ADNA_GIT" ]] || { echo "Missing ADNA_GIT: $ADNA_GIT"; exit 1; }
for sample in sed6 sed8; do
 [[ -s "$INPUT_RUN/extracted/$sample/Triticum_k35.fastq.gz" && -s "$INPUT_RUN/qc/${sample}_extraction.tsv" ]] || { echo "Missing k35 extraction: $sample"; exit 1; }
done
printf 'parameter\tvalue\nmode\t%s\ninput_run\t%s\nreference\t%s\nvcf\t%s\nbed\t%s\npop\t%s\npca_script\t%s\nsummary_script\t%s\nadna_git\t%s\nmin_mapq\t%s\nmin_baseq\t%s\n' "$MODE" "$INPUT_RUN" "$REF" "$VCF" "$BED" "$POP" "$PCA_SCRIPT" "$SUMMARY_SCRIPT" "$ADNA_GIT" "$MIN_MAPQ" "$MIN_BASEQ" > "$OUT/qc/parameters.tsv"
cp "$PCA_SCRIPT" "$OUT/qc/original_PCA_helper.sh"
cp "$SUMMARY_SCRIPT" "$OUT/qc/original_summary_helper.py"
# These copies + report allow review of hidden dependencies, -T semantics, SNP calling.
{
 echo '=== Explicit dependencies ==='
 for cmd in bwa samtools python3 bcftools smartpca; do command -v "$cmd" || true; done
 echo '=== Existing Svevo index ==='
 for ext in amb ann bwt pac sa fai; do ls -lh "$REF.$ext" 2>/dev/null || true; done
 echo '=== PCA helper (first 320 lines) ==='; sed -n '1,320p' "$PCA_SCRIPT"
 echo '=== Summary helper (first 240 lines) ==='; sed -n '1,240p' "$SUMMARY_SCRIPT"
 echo '=== Repository shell entry points ==='; find "$ADNA_GIT" -maxdepth 3 -type f -name '*.sh' | head -n 60 || true
 echo '=== BED first rows ==='; head -n 5 "$BED"
 echo '=== Population first rows ==='; head -n 5 "$POP"
} > "$OUT/qc/preflight_report.txt"
if [[ "$MODE" == preflight ]]; then
 echo "Preflight completed. Review qc/preflight_report.txt and both copied helpers before MODE=run."
 date -Is > "$OUT/PREFLIGHT_COMPLETE"; exit 0
fi
[[ "$MODE" == run ]] || { echo 'MODE must be preflight or run'; exit 1; }
for cmd in bwa samtools python3 bcftools smartpca; do command -v "$cmd" >/dev/null || { echo "Missing $cmd in active environment"; exit 1; }; done
# No automatic huge index build: avoid repeating the storage failure of script 12.
BWA_PREFIX=${BWA_PREFIX:-$REF}
for ext in amb ann bwt pac sa; do
 [[ -s "$BWA_PREFIX.$ext" ]] || { echo "Missing Svevo index: $BWA_PREFIX.$ext; supply BWA_PREFIX for an existing index of EXACTLY this FASTA"; exit 1; }
done
LOCALREF="$OUT/reference/reference.fa"; ln -s "$REF" "$LOCALREF"
for ext in amb ann bwt pac sa; do ln -s "$BWA_PREFIX.$ext" "$LOCALREF.$ext"; done
if [[ -s "$REF.fai" ]]; then ln -s "$REF.fai" "$LOCALREF.fai"; else samtools faidx "$LOCALREF"; fi
# Check BED contigs/ranges against Svevo; no guessed renaming/liftover.
python3 - "$LOCALREF.fai" "$BED" <<'PY'
import sys
lengths={f[0]:int(f[1]) for f in (l.split('\t') for l in open(sys.argv[1]))}
with open(sys.argv[2]) as h:
    for n,line in enumerate(h,1):
        if not line.strip() or line.startswith(('#','track','browser')): continue
        f=line.split(); c=f[0]; a,b=map(int,f[1:3])
        if c not in lengths or not (0<=a<b<=lengths[c]):
            raise ValueError(f'BED/Svevo incompatibility line {n}: {line.strip()}')
PY
bcftools view -h "$VCF" > "$OUT/qc/modern_vcf_header.txt"
for sample in sed6 sed8; do
 d="$OUT/mapping/$sample"; mkdir -p "$d"
 fq="$INPUT_RUN/extracted/$sample/Triticum_k35.fastq.gz"
 # Same alignment settings as original PCA pipeline (not script 12).
 bwa aln -l 1024 -n 0.01 -o 2 -t "$THREADS" "$LOCALREF" "$fq" > "$d/reads.sai" 2> "$d/bwa_aln.log"
 bwa samse "$LOCALREF" "$d/reads.sai" "$fq" 2> "$d/bwa_samse.log" \
   | samtools view -u -F 2308 -q "$MIN_MAPQ" - \
   | samtools sort -@ 4 -m 2G -T "$d/sorttmp" -o "$d/mapped.bam" -
 samtools index -c "$d/mapped.bam"
 samtools flagstat "$d/mapped.bam" > "$d/flagstat.txt"
 # No extra deduplication: original input FASTQs already deduplicated.
 # No mapDamage rescaling here; explicit deviation from original pipeline.
 # Full pileup retained to match original helper interface (BED supplied to helper).
 samtools mpileup -B -s -q "$MIN_MAPQ" -Q "$MIN_BASEQ" -f "$LOCALREF" "$d/mapped.bam" > "$OUT/pileup/${sample}_merged.pileup" 2> "$d/mpileup.log"
 [[ -s "$OUT/pileup/${sample}_merged.pileup" ]] || { echo "No usable pileup for $sample"; exit 1; }
 python3 "$SUMMARY_SCRIPT" -i "$OUT/pileup/${sample}_merged.pileup" -P 2 \
   -o "$OUT/pileup/${sample}_merged_PotentialSNPs.txt" -s "$OUT/pileup/${sample}_merged_chr_stats.txt" > "$d/summary.log" 2>&1
 # Separate workdir isolates helper temporary filenames.
 work="$OUT/PCA/$sample"; mkdir -p "$work"
 ( cd "$work"
   bash "$PCA_SCRIPT" -v "$VCF" -p "$POP" -t "$BED" \
     -m "$OUT/pileup/${sample}_merged.pileup" -s "$OUT/pileup/${sample}_merged_PotentialSNPs.txt" \
     -n "$sample" -c aDNA -o "$work" -P "$ADNA_GIT" -T
 ) > "$work/helper.log" 2>&1
 for ext in geno ind snp; do
   [[ -s "$work/$sample.PCAmatrix_withSample.$ext" ]] || { echo "Missing helper output: $sample/$ext"; exit 1; }
 done
 rm -f "$d/reads.sai"
done
# Shared modern-only axes: merge original EIGENSTRAT outputs only if compatible.
cat > "$OUT/qc/merge_eigenstrat.py" <<'PY'
from pathlib import Path
import sys
root=Path(sys.argv[1]); dest=root/'PCA/joint'; dest.mkdir()
samples=['sed6','sed8']; data=[]
for sample in samples:
    prefix=root/'PCA'/sample/f'{sample}.PCAmatrix_withSample'
    inds=[l.split() for l in Path(str(prefix)+'.ind').read_text().splitlines() if l.strip()]
    anc=[i for i,r in enumerate(inds) if r[0]==sample or r[2]=='aDNA']
    if len(anc)!=1: raise ValueError(f'Expected exactly one ancient individual in {sample}.ind')
    ai=anc[0]; modern=[r for i,r in enumerate(inds) if i!=ai]
    records={}
    with open(str(prefix)+'.snp') as sn, open(str(prefix)+'.geno') as ge:
        from itertools import zip_longest
        for sl,gl in zip_longest(sn,ge):
            if sl is None or gl is None: raise ValueError('SNP/GENO line count mismatch')
            sf=sl.split(); g=gl.strip()
            if len(g)!=len(inds) or set(g)-set('0129'): raise ValueError('Not text EIGENSTRAT genotypes')
            if sf[0] in records: raise ValueError('Duplicate SNP ID')
            records[sf[0]]=(sf,g[:ai]+g[ai+1:],g[ai])
    data.append((modern,records))
if data[0][0]!=data[1][0]: raise ValueError('Modern individuals/order/populations differ; refusing joint PCA')
modern=data[0][0]
if any(r[2]=='aDNA' for r in modern): raise ValueError('Ambiguous modern population labels')
with (dest/'joint.ind').open('w') as h:
    for row in modern: h.write(' '.join(row)+'\n')
    for sample in samples: h.write(f'{sample} U aDNA\n')
pops=sorted({r[2] for r in modern})
(dest/'modern.poplist').write_text('\n'.join(pops)+'\n')
kept=0; observed=[0,0]
with (dest/'joint.snp').open('w') as sn, (dest/'joint.geno').open('w') as ge:
    for key,(sf,mg,ag) in data[0][1].items():
        if key not in data[1][1]: continue
        sf2,mg2,ag2=data[1][1][key]
        if sf!=sf2 or mg!=mg2: raise ValueError('Modern SNP coding/genotypes differ: '+key)
        sn.write(' '.join(sf)+'\n'); ge.write(mg+ag+ag2+'\n'); kept+=1
        observed[0]+=ag!='9'; observed[1]+=ag2!='9'
with (dest/'overlap_QC.tsv').open('w') as h:
    h.write('sample\tshared_panel_SNPs\tancient_nonmissing_SNPs\n')
    for sample,n in zip(samples,observed): h.write(f'{sample}\t{kept}\t{n}\n')
if not kept or min(observed)<100:
    raise ValueError(f'Too few observed SNPs for a defensible projection: {observed}; QC saved; no joint PCA run')
params={'genotypename':dest/'joint.geno','snpname':dest/'joint.snp','indivname':dest/'joint.ind',
'evecoutname':dest/'joint.evec','evaloutname':dest/'joint.eval','poplistname':dest/'modern.poplist',
'numoutevec':4,'lsqproject':'YES','numoutlieriter':0}
(dest/'smartpca.par').write_text(''.join(f'{k}: {v}\n' for k,v in params.items()))
PY
python3 "$OUT/qc/merge_eigenstrat.py" "$OUT"
smartpca -p "$OUT/PCA/joint/smartpca.par" > "$OUT/PCA/joint/smartpca.log" 2>&1
[[ -s "$OUT/PCA/joint/joint.evec" ]] || exit 1
date -Is > "$OUT/SUCCESS"
echo "Completed: $OUT; joint PCA coordinates in PCA/joint/joint.evec"
