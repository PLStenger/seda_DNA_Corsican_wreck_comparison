#!/usr/bin/env bash
#SBATCH --job-name=triticum_cascade_PCA
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --partition=smp
#SBATCH --mem=120G
#SBATCH --output=/home/plstenge/seda_DNA_Corsican_wreck_comparison/00_scripts/13_triticum_PCA_k35_k29_k25_%j.out
#SBATCH --error=/home/plstenge/seda_DNA_Corsican_wreck_comparison/00_scripts/13_triticum_PCA_k35_k29_k25_%j.err
# Cascade k35 -> k29 -> k25, merged shotgun FASTQs; no ITS selection.
# Re-extract Triticum from original Kraken outputs, as in script 10; no Kraken rerun.
# Modern-only smartpca axes, joint projection of sed6 and sed8, R/ggplot2 output.
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
KRAKEN_RUN=${KRAKEN_RUN:-${1:-$PROJECT/09_kraken_cascade_k35_k29_k25/runs/run_20261030_162223_9307771}}
INPUT_DIR=${INPUT_DIR:-/storage/groups/gdec/shared_paleo/CapoSagro/fastp}
MIN_OBS_SNPS=${MIN_OBS_SNPS:-100}
REF=${REF:-/storage/groups/gdec/shared_paleo/genomes_REF/160802_Svevo_v2_pseudomolecules.fasta}
VCF=${VCF:-/storage/groups/gdec/shared_paleo/VCFWhealbi/WhealbiV1_4x_all_chr.vcf.gz}
POP=${POP:-/storage/groups/gdec/shared_paleo/VCFWhealbi/Whealbiv1_popOrigin.txt}
BED=${BED:-/storage/groups/gdec/shared_paleo/VCFWhealbi/WhealbiV1_4x_target_positions.bed}
ADNA_GIT=${ADNA_GIT:-$PIPELINE_DIR/adna_variousscripts-master}
THREADS=${SLURM_CPUS_PER_TASK:-36}
MIN_MAPQ=${MIN_MAPQ:-20}
MIN_BASEQ=${MIN_BASEQ:-20}
OUT="$PROJECT/13_triticum_PCA_k35_k29_k25/runs/run_$(date +%Y%m%d_%H%M%S)_${SLURM_JOB_ID:-manual_$$}"
mkdir -p "$OUT"/{qc,reference,extracted,mapping,pileup,PCA}
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
[[ "$MIN_OBS_SNPS" =~ ^[1-9][0-9]*$ ]] || { echo 'Invalid MIN_OBS_SNPS'; exit 1; }
for sample in sed6 sed8; do
 [[ -s "$INPUT_DIR/clean_${sample}_concat_dedup_fastp_merged.fastq.gz" ]] || { echo "Missing original FASTQ: $sample"; exit 1; }
 for stage in k35 k29 k25; do
  for ext in kraken report; do
   [[ -r "$KRAKEN_RUN/kraken/$stage/${sample}_${stage}.$ext" ]] || { echo "Missing $sample/$stage/$ext"; exit 1; }
  done
 done
done
printf 'parameter\tvalue\nmode\t%s\nkraken_run\t%s\nreference\t%s\nvcf\t%s\nbed\t%s\npop\t%s\npca_script\t%s\nsummary_script\t%s\nadna_git\t%s\nmin_mapq\t%s\nmin_baseq\t%s\n' "$MODE" "$KRAKEN_RUN" "$REF" "$VCF" "$BED" "$POP" "$PCA_SCRIPT" "$SUMMARY_SCRIPT" "$ADNA_GIT" "$MIN_MAPQ" "$MIN_BASEQ" > "$OUT/qc/parameters.tsv"
printf 'input_dir\t%s\nprofile\tk35_k29_k25\nmin_obs_snps\t%s\n' "$INPUT_DIR" "$MIN_OBS_SNPS" >> "$OUT/qc/parameters.tsv"
cp "$PCA_SCRIPT" "$OUT/qc/original_PCA_helper.sh"
cp "$SUMMARY_SCRIPT" "$OUT/qc/original_summary_helper.py"
# These copies + report allow review of hidden dependencies, -T semantics, SNP calling.
{
 echo '=== Explicit dependencies ==='
 for cmd in bwa samtools python3 bcftools smartpca Rscript gzip; do command -v "$cmd" || true; done
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
for cmd in bwa samtools python3 bcftools smartpca Rscript gzip; do command -v "$cmd" >/dev/null || { echo "Missing $cmd in active environment"; exit 1; }; done
Rscript -e 'if (!requireNamespace("ggplot2", quietly=TRUE)) stop("Install ggplot2 in the R used by this environment")'
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
# Reconstruct exact cumulative Triticum input; script 10 deletes extracted FASTQs by default.
printf 'reference_id\ttaxon\ttaxid\treference_path\nTriticum\tTriticum\t4564\tNA\n' > "$OUT/qc/targets.tsv"
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
for sample in sed6 sed8; do
 python3 "$OUT/qc/extract_cascade.py" "$KRAKEN_RUN" "$INPUT_DIR" "$OUT" "$sample"
 cat "$OUT/extracted/$sample/Triticum_k35.fastq" \
     "$OUT/extracted/$sample/Triticum_k29.fastq" \
     "$OUT/extracted/$sample/Triticum_k25.fastq" \
     | gzip -c > "$OUT/extracted/$sample/Triticum_k35_k29_k25.fastq.gz"
 [[ -s "$OUT/extracted/$sample/Triticum_k35_k29_k25.fastq.gz" ]]
 n=$(awk -F '\t' '$2=="Triticum" && $3=="k35_k29_k25" {print $4}' "$OUT/qc/${sample}_extraction.tsv")
 [[ "$n" =~ ^[1-9][0-9]*$ ]] || { echo "No Triticum cascade reads: $sample"; exit 1; }
 rm -- "$OUT/extracted/$sample/"Triticum_k{35,29,25}.fastq
done
for sample in sed6 sed8; do
 d="$OUT/mapping/$sample"; mkdir -p "$d"
 fq="$OUT/extracted/$sample/Triticum_k35_k29_k25.fastq.gz"
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
root=Path(sys.argv[1]); min_obs=int(sys.argv[2]); dest=root/'PCA/joint'; dest.mkdir()
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
if not kept or min(observed)<min_obs:
    raise ValueError(f'Observed SNPs below configured technical threshold (not a reliability guarantee): {observed}; QC saved; no joint PCA run')
params={'genotypename':dest/'joint.geno','snpname':dest/'joint.snp','indivname':dest/'joint.ind',
'evecoutname':dest/'joint.evec','evaloutname':dest/'joint.eval','poplistname':dest/'modern.poplist',
'numoutevec':4,'lsqproject':'YES','numoutlieriter':0}
(dest/'smartpca.par').write_text(''.join(f'{k}: {v}\n' for k,v in params.items()))
PY
python3 "$OUT/qc/merge_eigenstrat.py" "$OUT" "$MIN_OBS_SNPS"
smartpca -p "$OUT/PCA/joint/smartpca.par" > "$OUT/PCA/joint/smartpca.log" 2>&1
[[ -s "$OUT/PCA/joint/joint.evec" ]] || exit 1
cat > "$OUT/qc/plot_PCA_whealbi.R" <<'RSCRIPT'
#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript plot_PCA_whealbi.R <PCA/joint>")
work <- normalizePath(args[1], mustWork = TRUE)
ind <- read.table(file.path(work, "joint.ind"), header = FALSE,
                  col.names = c("sample", "sex", "population"),
                  colClasses = "character", comment.char = "", quote = "")
evec <- read.table(file.path(work, "joint.evec"), header = FALSE,
                   comment.char = "#", quote = "", stringsAsFactors = FALSE)
npc <- ncol(evec) - 2L
if (npc < 2L) stop("Fewer than two PCs in joint.evec")
names(evec) <- c("sample", paste0("PC", seq_len(npc)), "evec_population")
if (anyDuplicated(ind$sample) || anyDuplicated(evec$sample)) stop("Duplicate individual IDs")
if (!setequal(ind$sample, evec$sample)) stop("IND/EVEC individual mismatch")
evec$population <- ind$population[match(evec$sample, ind$sample)]
evec$type <- ifelse(evec$sample %in% c("sed6", "sed8"), "Ancient", "Modern")
if (sum(evec$type == "Ancient") != 2L) stop("Expected sed6 and sed8")
for (pc in paste0("PC", seq_len(npc))) {
 evec[[pc]] <- as.numeric(evec[[pc]])
 if (any(!is.finite(evec[[pc]]))) stop("Non-finite coordinates: ", pc)
}
modern <- evec[evec$type == "Modern", , drop = FALSE]
ancient <- evec[evec$type == "Ancient", , drop = FALSE]
qc <- read.delim(file.path(work, "overlap_QC.tsv"))
ancient$n_snps <- qc$ancient_nonmissing_SNPs[match(ancient$sample, qc$sample)]
if (anyNA(ancient$n_snps)) stop("Missing ancient SNP counts")
ancient$label <- paste0(ancient$sample, " (", format(ancient$n_snps, big.mark = ",", trim = TRUE), " SNPs)")
write.table(evec, file.path(work, "PCA_coordinates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ancient, file.path(work, "PCA_ancient_coordinates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
p <- ggplot() +
 geom_hline(yintercept = 0, colour = "grey90", linewidth = 0.3) +
 geom_vline(xintercept = 0, colour = "grey90", linewidth = 0.3) +
 geom_point(data = modern, aes(PC1, PC2, colour = population), size = 2, alpha = 0.7) +
 geom_point(data = ancient, aes(PC1, PC2, shape = sample), size = 4.5,
            colour = "black", fill = "#FFD34E", stroke = 1.1) +
 scale_shape_manual(values = c(sed6 = 23, sed8 = 24), name = "Ancient DNA") +
 labs(title = "Triticum: ancient DNA projected onto Whealbi",
      subtitle = paste0("Axes fitted on ", nrow(modern), " modern individuals; cascade k35 + k29 + k25"),
      x = "PC1", y = "PC2", colour = "Whealbi population",
      caption = "Ancient samples excluded from axis estimation. Labels: non-missing SNPs before smartpca filtering.") +
 coord_equal(clip = "off") +
 scale_x_continuous(expand = expansion(mult = 0.15)) +
 scale_y_continuous(expand = expansion(mult = 0.15)) +
 theme_bw(base_size = 12) +
 theme(panel.grid.minor = element_blank(), legend.position = "right",
       plot.title = element_text(face = "bold"),
       plot.margin = margin(12, 35, 12, 12))
if (requireNamespace("ggrepel", quietly = TRUE)) {
 p <- p + ggrepel::geom_text_repel(data = ancient, aes(PC1, PC2, label = label),
       seed = 42, colour = "black", size = 4, min.segment.length = 0, max.overlaps = Inf)
} else {
 p <- p + geom_text(data = ancient, aes(PC1, PC2, label = label),
                    nudge_y = 0.035 * max(diff(range(evec$PC2)), 0.01), size = 3.5)
}
ggsave(file.path(work, "PCA_Whealbi_sed6_sed8.pdf"), p, width = 12, height = 8, device = "pdf")
ggsave(file.path(work, "PCA_Whealbi_sed6_sed8.png"), p, width = 12, height = 8, dpi = 300)
saveRDS(p, file.path(work, "PCA_Whealbi_sed6_sed8.rds"))
writeLines(capture.output(sessionInfo()), file.path(work, "R_sessionInfo.txt"))
RSCRIPT
Rscript "$OUT/qc/plot_PCA_whealbi.R" "$OUT/PCA/joint" > "$OUT/PCA/joint/ggplot.log" 2>&1
date -Is > "$OUT/SUCCESS"
echo "Completed: $OUT; joint PCA coordinates in PCA/joint/joint.evec"
