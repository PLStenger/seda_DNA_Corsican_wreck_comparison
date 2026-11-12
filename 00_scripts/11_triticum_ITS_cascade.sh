#!/usr/bin/env bash
#SBATCH --job-name=triticum_ITS_cascade
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --partition=smp
#SBATCH --mem=32G
#SBATCH --output=/home/plstenge/seda_DNA_Corsican_wreck_comparison/00_scripts/11_triticum_ITS_%j.out
#SBATCH --error=/home/plstenge/seda_DNA_Corsican_wreck_comparison/00_scripts/11_triticum_ITS_%j.err

# sbatch 11_triticum_ITS_cascade.sh [KRAKEN_RUN]
# Optional environment variables: INDEX, CONDA_ENV, MIN_MAPQ, THREADS,
# PROJECT, INPUT_DIR, PIPELINE_DIR, CONTIG_MAP, SAMPLES_STR.
# No Kraken rerun, no trimming, no additional deduplication, no index modification.
set -Eeo pipefail
module load conda/4.12.0
source ~/.bashrc
conda activate "${CONDA_ENV:-mapdamage_py39}"
set -Eeuo pipefail
export LC_ALL=C
trap 'printf "ERROR line %s: %s\n" "$LINENO" "$BASH_COMMAND" >&2' ERR
PROJECT=${PROJECT:-/home/plstenge/seda_DNA_Corsican_wreck_comparison}
INPUT_DIR=${INPUT_DIR:-/storage/groups/gdec/shared_paleo/CapoSagro/fastp}
PIPELINE_DIR=${PIPELINE_DIR:-/home/plstenge/adna_analysis_new-main}
INDEX=${INDEX:-$PIPELINE_DIR/Bowtie2_index/ITS_ABDGH}
KRAKEN_RUN=${1:-/home/plstenge/seda_DNA_Corsican_wreck_comparison/09_kraken_cascade_k35_k29_k25/runs/run_20261030_162223_9307771}
THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-36}}
MIN_MAPQ=${MIN_MAPQ:-20}
read -r -a SAMPLES <<< "${SAMPLES_STR:-sed6 sed8}"
PROFILES=(k35 k29 k25 k35_k29_k25)
[[ "$MIN_MAPQ" =~ ^[0-9]+$ && "$MIN_MAPQ" -le 255 ]] || { echo 'Invalid MIN_MAPQ'; exit 1; }
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || exit 1
for cmd in python3 bowtie2 bowtie2-inspect samtools; do
    command -v "$cmd" >/dev/null || { echo "Missing executable in active environment: $cmd" >&2; exit 1; }
done
for sample in "${SAMPLES[@]}"; do
    [[ -s "$INPUT_DIR/clean_${sample}_concat_dedup_fastp_merged.fastq.gz" ]] || { echo "Missing FASTQ: $sample"; exit 1; }
    for stage in k35 k29 k25; do
        for ext in kraken report; do
            [[ -r "$KRAKEN_RUN/kraken/$stage/${sample}_${stage}.$ext" ]] || { echo "Missing Kraken input: $sample/$stage/$ext"; exit 1; }
        done
    done
done
RUN_ID="run_$(date +%Y%m%d_%H%M%S)_${SLURM_JOB_ID:-manual_$$}"
OUT="$PROJECT/11_triticum_ITS_cascade/runs/$RUN_ID"
mkdir -p "$OUT"/{qc,extracted,results}
exec > >(tee -a "$OUT/pipeline.log") 2>&1
printf 'Output: %s\nKraken run: %s\nIndex: %s\n' "$OUT" "$KRAKEN_RUN" "$INDEX"
# sbatch may expose its spool copy through $0; this still saves the submitted script.
cp -- "$0" "$OUT/qc/submitted_script.sh"
conda list > "$OUT/qc/conda_packages.txt"
bowtie2 --version > "$OUT/qc/bowtie2_version.txt"
samtools --version > "$OUT/qc/samtools_version.txt"
printf 'parameter\tvalue\nkraken_run\t%s\ninput_dir\t%s\nindex\t%s\nmin_mapq\t%s\nthreads\t%s\nmode\t--very-sensitive-local; single-end merged; seed=0\n' \
    "$KRAKEN_RUN" "$INPUT_DIR" "$INDEX" "$MIN_MAPQ" "$THREADS" > "$OUT/qc/parameters.tsv"
# Inspect actual contig names, rather than guessing subgenomes from substrings.
bowtie2-inspect -n "$INDEX" > "$OUT/qc/index_contigs.txt" 2> "$OUT/qc/index_inspect.log"
[[ -s "$OUT/qc/index_contigs.txt" ]] || { echo 'Empty index'; exit 1; }
if [[ -n "${CONTIG_MAP:-}" ]]; then
    cp -- "$CONTIG_MAP" "$OUT/qc/contig_subgenome.tsv"
else
    printf 'contig\tsubgenome\nrDNA_6BS_Langdon_extraction_3\tB\nrDNA_5DS_tauschi_extraction_3\tD\nMZ547657_extraction_2\tG\nrDNA_1BS_Langdon_extraction_3\tB\nMZ547658_extraction_2\tA\nHordeum\tHord\n' > "$OUT/qc/contig_subgenome.tsv"
fi
cat > "$OUT/qc/worker.py" <<'PY'
import csv, gzip, sys
from pathlib import Path
from contextlib import ExitStack
from collections import Counter

def load_map(root):
    with (root/'qc/contig_subgenome.tsv').open() as h:
        rows=list(csv.DictReader(h, delimiter='\t'))
    mapping={r['contig']:r['subgenome'] for r in rows}
    if len(mapping)!=len(rows) or not mapping or any(not k or not v for k,v in mapping.items()):
        raise ValueError('Invalid/duplicate entries in contig_subgenome.tsv')
    names=(root/'qc/index_contigs.txt').read_text().splitlines()
    missing=set(names)-set(mapping)
    if missing:
        raise ValueError('Unmapped index contigs: '+', '.join(sorted(missing))+
                         '. Supply CONTIG_MAP: tab-separated header contig/subgenome; never guess labels.')
    return {k:mapping[k] for k in names}

def extract(root, run, fastp, sample):
    stages=('k35','k29','k25')
    allowed={}
    for stage in stages:
        stack=[]; taxa={'4564'}; found=False
        with (run/'kraken'/stage/f'{sample}_{stage}.report').open() as h:
            for line in h:
                f=line.rstrip('\r\n').split('\t')
                if len(f)!=6:
                    raise ValueError(f'Expected standard 6-column Kraken report: {sample}/{stage}')
                taxid=f[4].strip(); name=f[5]; depth=len(name)-len(name.lstrip(' '))
                while stack and stack[-1][0]>=depth:
                    stack.pop()
                if taxid=='4564': found=True
                if taxid=='4564' or any(t=='4564' for _,t in stack): taxa.add(taxid)
                stack.append((depth,taxid))
        allowed[stage]=taxa
        if not found: print(f'NOTE: Triticum node absent from {sample}/{stage} report', flush=True)
    folder=root/'extracted'/sample; folder.mkdir()
    counts=Counter(); seen=Counter(); classified=Counter(); n=0
    with ExitStack() as ctx:
        kr={s:ctx.enter_context((run/'kraken'/s/f'{sample}_{s}.kraken').open()) for s in stages}
        outputs={s:ctx.enter_context(gzip.open(folder/f'Triticum_{s}.fastq.gz','wt')) for s in stages}
        combined=ctx.enter_context(gzip.open(folder/'Triticum_k35_k29_k25.fastq.gz','wt'))
        provenance=ctx.enter_context(gzip.open(folder/'read_assignment.tsv.gz','wt'))
        provenance.write('read_id\tstage\ttaxid\n')
        fq=ctx.enter_context(gzip.open(fastp/f'clean_{sample}_concat_dedup_fastp_merged.fastq.gz','rt'))
        while True:
            header=fq.readline()
            if not header: break
            seq,plus,qual=fq.readline(),fq.readline(),fq.readline(); n+=1
            if not header.startswith('@') or not header[1:].split() or not plus.startswith('+') or not qual or not seq.strip() or len(seq.rstrip('\r\n'))!=len(qual.rstrip('\r\n')):
                raise ValueError(f'Invalid FASTQ record: {sample}/{n}')
            name=header[1:].split()[0]
            for stage in stages:
                f=kr[stage].readline().rstrip('\r\n').split('\t')
                if len(f)<5 or f[0] not in ('C','U') or f[1]!=name:
                    raise ValueError(f'Kraken/FASTQ desynchronization: {sample}/{stage}/{n}/{name}')
                if f[3]!=str(len(seq.rstrip('\r\n'))):
                    raise ValueError(f'Unexpected Kraken read length: {sample}/{stage}/{name}')
                seen[stage]+=1
                if f[0]=='C':
                    classified[stage]+=1
                    if f[2] in allowed[stage]:
                        record=header+seq+plus+qual
                        outputs[stage].write(record); combined.write(record)
                        provenance.write(f'{name}\t{stage}\t{f[2]}\n'); counts[stage]+=1
                    break
        for s in stages:
            if kr[s].readline(): raise ValueError(f'Extra Kraken records: {sample}/{s}')
        if seen['k29']!=seen['k35']-classified['k35'] or seen['k25']!=seen['k29']-classified['k29']:
            raise ValueError('Invalid cascade counts')
    with (root/'qc'/f'{sample}_extraction.tsv').open('w') as h:
        h.write('sample\tprofile\textracted_reads\n')
        for s in (*stages,'k35_k29_k25'):
            h.write(f'{sample}\t{s}\t{sum(counts.values()) if s=="k35_k29_k25" else counts[s]}\n')
    with (root/'qc'/f'{sample}_cascade.tsv').open('w') as h:
        h.write('stage\tinput_reads\tclassified_reads\n')
        for s in stages: h.write(f'{s}\t{seen[s]}\t{classified[s]}\n')
    print(f'{sample}: {n} input reads; {sum(counts.values())} Triticum reads',flush=True)

def stats(root, sample, profile, n, q, folder):
    mapping=load_map(root); counts=Counter(); groups=Counter(); total=mapped=passed=0
    with gzip.open(folder/'primary_assignments.tsv.gz','wt') as h:
        h.write('read_id\tcontig\tsubgenome\tmapq\tAS\tXS\tpass_mapq\n')
        for line in sys.stdin:
            if line.startswith('@'): continue
            f=line.rstrip('\n').split('\t'); flag=int(f[1]); total+=1
            if flag & (256|2048): raise ValueError('Unexpected nonprimary alignment')
            if flag & 4: continue
            contig=f[2]; group=mapping[contig]; mq=int(f[4]); keep=mq>=q
            tags={t.split(':',2)[0]:t.split(':',2)[2] for t in f[11:] if t.count(':')>=2}
            mapped+=1; passed+=keep
            counts[contig,'all']+=1; groups[group,'all']+=1
            if keep: counts[contig,'filtered']+=1; groups[group,'filtered']+=1
            h.write(f'{f[0]}\t{contig}\t{group}\t{mq}\t{tags.get("AS","NA")}\t{tags.get("XS","NA")}\t{int(keep)}\n')
    if total!=n: raise ValueError(f'BAM/FASTQ count mismatch: {total} != {n}')
    with (folder/'per_contig.tsv').open('w') as h:
        h.write('contig\tsubgenome\tprimary_mapped_all_mapq\tprimary_mapped_mapq_ge_threshold\n')
        for c,g in mapping.items(): h.write(f'{c}\t{g}\t{counts[c,"all"]}\t{counts[c,"filtered"]}\n')
    def pct(a,b): return f'{100*a/b:.6f}' if b else 'NA'
    with (root/'subgenome_summary.tsv').open('a') as h:
        for g in sorted(set(mapping.values())):
            a,b=groups[g,'all'],groups[g,'filtered']
            h.write('\t'.join(map(str,[sample,profile,g,n,a,b,q,pct(a,n),pct(b,n),pct(b,passed)]))+'\n')
    with (root/'mapping_summary.tsv').open('a') as h:
        h.write('\t'.join(map(str,[sample,profile,n,total,mapped,passed,q,pct(mapped,n),pct(passed,n),'NO_READS' if not n else 'OK']))+'\n')

mode=sys.argv[1]; root=Path(sys.argv[2])
if mode=='validate': load_map(root)
elif mode=='extract': extract(root,Path(sys.argv[3]),Path(sys.argv[4]),sys.argv[5])
elif mode=='stats': stats(root,sys.argv[3],sys.argv[4],int(sys.argv[5]),int(sys.argv[6]),Path(sys.argv[7]))
else: raise ValueError(mode)
PY
python3 "$OUT/qc/worker.py" validate "$OUT"
printf 'sample\tprofile\textracted_reads\tprimary_records\tprimary_mapped_all_mapq\tprimary_mapped_mapq_ge_threshold\tmin_mapq\tmapping_percent_all\tmapping_percent_filtered\tstatus\n' > "$OUT/mapping_summary.tsv"
printf 'sample\tprofile\tsubgenome\textracted_reads\tprimary_mapped_all_mapq\tprimary_mapped_mapq_ge_threshold\tmin_mapq\tpercent_input_all\tpercent_input_filtered\tpercent_filtered_mapped\n' > "$OUT/subgenome_summary.tsv"
for sample in "${SAMPLES[@]}"; do
    python3 "$OUT/qc/worker.py" extract "$OUT" "$KRAKEN_RUN" "$INPUT_DIR" "$sample"
    for profile in "${PROFILES[@]}"; do
        d="$OUT/results/$sample/$profile"; mkdir -p "$d"
        n=$(awk -F '\t' -v p="$profile" 'NR>1 && $2==p {print $3}' "$OUT/qc/${sample}_extraction.tsv")
        [[ "$n" =~ ^[0-9]+$ ]] || exit 1
        if [[ "$n" == 0 ]]; then
            python3 "$OUT/qc/worker.py" stats "$OUT" "$sample" "$profile" 0 "$MIN_MAPQ" "$d" < /dev/null
            continue
        fi
        printf 'Mapping %s / %s: %s reads\n' "$sample" "$profile" "$n"
        # Keep unmapped records to verify the total input count. Default one primary/read.
        bowtie2 -x "$INDEX" -U "$OUT/extracted/$sample/Triticum_${profile}.fastq.gz" \
            --very-sensitive-local --seed 0 -p "$THREADS" 2> "$d/bowtie2.log" \
            | samtools view -u -F 2304 - \
            | samtools sort -@ 2 -m 1G -T "$d/sorttmp" -o "$d/primary.sorted.bam" -
        samtools quickcheck "$d/primary.sorted.bam"
        samtools index "$d/primary.sorted.bam"
        samtools flagstat "$d/primary.sorted.bam" > "$d/flagstat.txt"
        samtools idxstats "$d/primary.sorted.bam" > "$d/idxstats_all.tsv"
        samtools view -b -F 2308 -q "$MIN_MAPQ" -o "$d/mapped.q${MIN_MAPQ}.bam" "$d/primary.sorted.bam"
        samtools index "$d/mapped.q${MIN_MAPQ}.bam"
        samtools idxstats "$d/mapped.q${MIN_MAPQ}.bam" > "$d/idxstats_q${MIN_MAPQ}.tsv"
        samtools view "$d/primary.sorted.bam" \
            | python3 "$OUT/qc/worker.py" stats "$OUT" "$sample" "$profile" "$n" "$MIN_MAPQ" "$d"
    done
done
date -Is > "$OUT/SUCCESS"
printf 'Completed: %s\n' "$OUT"
