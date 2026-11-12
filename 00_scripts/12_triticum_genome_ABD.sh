#!/usr/bin/env bash
#SBATCH --job-name=triticum_ABD
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --partition=smp
#SBATCH --mem=120G
#SBATCH --output=/home/plstenge/seda_DNA_Corsican_wreck_comparison/00_scripts/12_triticum_ABD_%j.out
#SBATCH --error=/home/plstenge/seda_DNA_Corsican_wreck_comparison/00_scripts/12_triticum_ABD_%j.err
# Usage: sbatch 12_triticum_genome_ABD.sh /path/to/SUCCESSFUL_11_run
# Inputs: all Kraken Triticum reads, NOT just ITS-mapped reads.
set -Eeo pipefail
module load conda/4.12.0
source ~/.bashrc
conda activate "${CONDA_ENV:-mapdamage_py39}"
set -Eeuo pipefail
export LC_ALL=C
trap 'printf "ERROR line %s: %s\n" "$LINENO" "$BASH_COMMAND" >&2' ERR
PROJECT=${PROJECT:-/home/plstenge/seda_DNA_Corsican_wreck_comparison}
REF=${REF:-/home/plstenge/genomes/Triticum_aestivum/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic.fna}
INPUT_RUN=${1:?Usage: sbatch 12_triticum_genome_ABD.sh /path/to/successful_11_run}
THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-36}}
MIN_MAPQ=${MIN_MAPQ:-25}
WINDOW=${WINDOW:-1000000}
read -r -a SAMPLES <<< "${SAMPLES_STR:-sed6 sed8}"
read -r -a PROFILES <<< "${PROFILES_STR:-k35_k29_k25}"
for v in THREADS WINDOW; do [[ ${!v} =~ ^[1-9][0-9]*$ ]] || exit 1; done
[[ $MIN_MAPQ =~ ^[0-9]+$ && $MIN_MAPQ -le 255 ]] || exit 1
[[ -s "$REF" && "$REF" != *.gz ]] || { echo 'Reference must be readable, uncompressed FASTA'; exit 1; }
for cmd in bwa samtools python3; do command -v "$cmd" >/dev/null || { echo "Missing: $cmd"; exit 1; }; done
for sample in "${SAMPLES[@]}"; do
 for profile in "${PROFILES[@]}"; do
  [[ -s "$INPUT_RUN/extracted/$sample/Triticum_${profile}.fastq.gz" ]] || { echo "Missing FASTQ: $sample/$profile"; exit 1; }
 done
done
OUT="$PROJECT/12_triticum_genome_ABD/runs/run_$(date +%Y%m%d_%H%M%S)_${SLURM_JOB_ID:-manual_$$}"
mkdir -p "$OUT"/{qc,reference,results}
exec > >(tee -a "$OUT/pipeline.log") 2>&1
printf 'Output: %s\nInput run: %s\nReference: %s\n' "$OUT" "$INPUT_RUN" "$REF"
cp -- "$0" "$OUT/qc/submitted_script.sh"
conda list > "$OUT/qc/conda_packages.txt"
samtools --version > "$OUT/qc/samtools_version.txt"
printf 'parameter\tvalue\ninput_run\t%s\nreference\t%s\nmin_mapq\t%s\nwindow_bp\t%s\nbwa_aln\t-n 0.08 -l 24 -k 2 -q 20\n' "$INPUT_RUN" "$REF" "$MIN_MAPQ" "$WINDOW" > "$OUT/qc/parameters.tsv"
LOCALREF="$OUT/reference/reference.fa"
ln -s "$REF" "$LOCALREF"
# Reuse complete existing BWA indexes; otherwise build locally, never beside source.
complete=1
for ext in amb ann bwt pac sa; do [[ -s "$REF.$ext" ]] || complete=0; done
if [[ $complete == 1 ]]; then
 for ext in amb ann bwt pac sa; do ln -s "$REF.$ext" "$LOCALREF.$ext"; done
else
 echo 'Building local BWA index (large wheat reference; may take hours).'
 bwa index -a bwtsw "$LOCALREF" > "$OUT/reference/bwa_index.log" 2>&1
fi
if [[ -s "$REF.fai" ]]; then ln -s "$REF.fai" "$LOCALREF.fai"; else samtools faidx "$LOCALREF"; fi
cat > "$OUT/qc/analyze.py" <<'PY'
import sys, csv, re, gzip
from pathlib import Path
from collections import Counter

def references(root):
    rows=[]
    with (root/'reference/reference.fa').open('rb') as fa, (root/'reference/reference.fa.fai').open() as ix:
        for line in ix:
            f=line.rstrip().split('\t'); name=f[0]; length=int(f[1]); offset=int(f[2])
            start=max(0,offset-1048576); fa.seek(start)
            block=fa.read(offset-start)
            heads=[x.decode() for x in block.splitlines() if x.startswith(b'>')]
            if not heads or heads[-1][1:].split()[0]!=name:
                raise ValueError('Cannot recover FASTA header for '+name)
            header=heads[-1][1:]
            match=re.search(r'\bchromosome\s+([1-7])([ABD])\b',header,re.I)
            chrom=match.group(1)+match.group(2).upper() if match else 'OTHER'
            group=match.group(2).upper() if match else 'OTHER'
            rows.append(dict(contig=name,length_bp=length,chromosome=chrom,subgenome=group,header=header))
    expected={f'{i}{g}' for i in range(1,8) for g in 'ABD'}
    observed=[r['chromosome'] for r in rows if r['chromosome']!='OTHER']
    if set(observed)!=expected or len(observed)!=21:
        raise ValueError('Expected exactly one chromosome for each 1A..7D; inspect reference headers')
    with (root/'qc/contig_metadata.tsv').open('w') as h:
        w=csv.DictWriter(h,fieldnames=rows[0].keys(),delimiter='\t'); w.writeheader(); w.writerows(rows)

root=Path(sys.argv[2])
if sys.argv[1]=='reference':
    references(root); sys.exit()
sample,profile,q,window,folder=sys.argv[3:]; q=int(q); window=int(window); folder=Path(folder)
with (root/'qc/contig_metadata.tsv').open() as h:
    meta={r['contig']:r for r in csv.DictReader(h,delimiter='\t')}
counts=Counter(); bases=Counter(); starts=Counter(); winbases=Counter(); hist=Counter()
n=0; mapped=0; kept=0; cigar_re=re.compile(r'(\d+)([MIDNSHP=X])')
for line in sys.stdin:
    if line.startswith('@'): continue
    f=line.rstrip().split('\t'); flag=int(f[1]); n+=1
    if flag & 2304: raise ValueError('Nonprimary record')
    if flag & 4: continue
    c=f[2]; mq=int(f[4]); mapped+=1
    if c not in meta: raise ValueError('Unknown contig '+c)
    hist[mq]+=1; filters=['all']
    if mq>=q: filters.append('filtered'); kept+=1
    position=int(f[3])-1
    for mode in filters:
        counts[c,mode]+=1; starts[c,position//window,mode]+=1
    cursor=position
    for length,op in cigar_re.findall(f[5]):
        length=int(length)
        if op in 'M=X':
            end=cursor+length
            for mode in filters: bases[c,mode]+=length
            while cursor<end:
                b=cursor//window; stop=min(end,(b+1)*window)
                for mode in filters: winbases[c,b,mode]+=stop-cursor
                cursor=stop
        elif op in 'DN': cursor+=length

def write(name,header,rows):
    with (folder/name).open('w') as h:
        w=csv.writer(h,delimiter='\t'); w.writerow(header); w.writerows(rows)
def ratio(a,b): return a/b if b else 'NA'
write('mapping_summary.tsv',['sample','profile','input_reads','mapped_all','mapped_filtered','min_mapq'],[[sample,profile,n,mapped,kept,q]])
write('mapq_histogram.tsv',['mapq','mapped_reads'],sorted(hist.items()))
rows=[]
for c,r in meta.items():
    L=int(r['length_bp'])
    for mode in ('all','filtered'):
        rows.append([c,r['chromosome'],r['subgenome'],L,mode,counts[c,mode],1e6*counts[c,mode]/L,bases[c,mode]/L])
write('per_chromosome.tsv',['contig','chromosome','subgenome','length_bp','filter','reads','reads_per_Mb','mean_depth'],rows)
rows=[]
for mode in ('all','filtered'):
    nuclear_reads=sum(counts[c,mode] for c,r in meta.items() if r['subgenome'] in 'ABD')
    nuclear_L=sum(int(r['length_bp']) for r in meta.values() if r['subgenome'] in 'ABD')
    for g in ('A','B','D','OTHER'):
        cc=[c for c,r in meta.items() if r['subgenome']==g]
        L=sum(int(meta[c]['length_bp']) for c in cc)
        N=sum(counts[c,mode] for c in cc); B=sum(bases[c,mode] for c in cc)
        expected=nuclear_reads*L/nuclear_L if g!='OTHER' else None
        rows.append([sample,profile,mode,g,L,N,ratio(N*1e6,L),ratio(B,L),ratio(100*N,nuclear_reads) if g!='OTHER' else 'NA',expected if expected is not None else 'NA',ratio(N,expected) if expected is not None else 'NA'])
write('subgenome_balance.tsv',['sample','profile','filter','subgenome','length_bp','reads','reads_per_Mb','mean_depth','percent_ABD_reads','expected_reads_length_model','observed_expected'],rows)
with gzip.open(folder/'windows.tsv.gz','wt') as h:
    w=csv.writer(h,delimiter='\t')
    w.writerow(['contig','chromosome','subgenome','start_0based','end_exclusive','filter','read_starts','aligned_bases','mean_depth'])
    for c,r in meta.items():
        L=int(r['length_bp'])
        for start in range(0,L,window):
            end=min(start+window,L); b=start//window
            for mode in ('all','filtered'):
                w.writerow([c,r['chromosome'],r['subgenome'],start,end,mode,starts[c,b,mode],winbases[c,b,mode],winbases[c,b,mode]/(end-start)])
# Optional plots; analysis remains usable if matplotlib is unavailable.
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    print('matplotlib unavailable: TSV outputs complete, figures skipped',file=sys.stderr)
else:
    groups='ABD'; colors=['#4477AA','#EE6677','#228833']
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for ax,mode in zip(axes,('all','filtered')):
        values=[]
        for g in groups:
            cc=[c for c,r in meta.items() if r['subgenome']==g]
            L=sum(int(meta[c]['length_bp']) for c in cc)
            values.append(sum(counts[c,mode] for c in cc)*1e6/L)
        ax.bar(list(groups),values,color=colors); ax.set_ylabel('Primary mapped reads / Mb')
        ax.set_title('All MAPQ' if mode=='all' else f'MAPQ >= {q}')
    fig.suptitle(f'{sample} / {profile}'); fig.tight_layout(); fig.savefig(folder/'subgenome_balance.png',dpi=180); plt.close(fig)
    fig,axes=plt.subplots(7,3,figsize=(16,15),squeeze=False)
    for c,r in meta.items():
        if r['subgenome'] not in groups: continue
        ax=axes[int(r['chromosome'][0])-1,groups.index(r['subgenome'])]; L=int(r['length_bp'])
        x=[]; y=[]
        for start in range(0,L,window):
            x.append(start/1e6); y.append(winbases[c,start//window,'filtered']/(min(start+window,L)-start))
        ax.plot(x,y,color=colors[groups.index(r['subgenome'])],linewidth=.8)
        ax.set_title(r['chromosome']); ax.set_xlabel('Position (Mb)'); ax.set_ylabel('Mean depth')
    ymax=max(ax.get_ylim()[1] for row in axes for ax in row)
    for row in axes:
        for ax in row: ax.set_ylim(0,ymax)
    fig.suptitle(f'{sample}: MAPQ >= {q}; windows = {window:,} bp'); fig.tight_layout(rect=(0,0,1,.97)); fig.savefig(folder/'chromosome_windows.png',dpi=160); plt.close(fig)
PY
python3 "$OUT/qc/analyze.py" reference "$OUT"
for sample in "${SAMPLES[@]}"; do
 for profile in "${PROFILES[@]}"; do
    d="$OUT/results/$sample/$profile"; mkdir -p "$d"
    fq="$INPUT_RUN/extracted/$sample/Triticum_${profile}.fastq.gz"
    echo "Mapping $sample / $profile"
    bwa aln -n 0.08 -l 24 -k 2 -q 20 -t "$THREADS" "$LOCALREF" "$fq" > "$d/reads.sai" 2> "$d/bwa_aln.log"
    bwa samse "$LOCALREF" "$d/reads.sai" "$fq" 2> "$d/bwa_samse.log" \
        | samtools view -u -F 2304 - \
        | samtools sort -@ 4 -m 2G -T "$d/sorttmp" -o "$d/primary.sorted.bam" -
    samtools quickcheck "$d/primary.sorted.bam"
    # CSI avoids the 512 Mb coordinate limit of BAI for wheat chromosomes.
    samtools index -c "$d/primary.sorted.bam"
    samtools flagstat "$d/primary.sorted.bam" > "$d/flagstat.txt"
    samtools idxstats "$d/primary.sorted.bam" > "$d/idxstats.tsv"
    samtools view -b -F 2308 -q "$MIN_MAPQ" -o "$d/mapped.q${MIN_MAPQ}.bam" "$d/primary.sorted.bam"
    samtools index -c "$d/mapped.q${MIN_MAPQ}.bam"
    samtools view "$d/primary.sorted.bam" \
       | python3 "$OUT/qc/analyze.py" stats "$OUT" "$sample" "$profile" "$MIN_MAPQ" "$WINDOW" "$d"
    rm -f "$d/reads.sai"
 done
done
python3 - "$OUT" <<'PY'
import sys,csv
from pathlib import Path
root=Path(sys.argv[1])
for name in ('subgenome_balance.tsv','mapping_summary.tsv'):
    files=sorted((root/'results').glob('*/*/'+name))
    with (root/name).open('w') as h:
        w=csv.writer(h,delimiter='\t'); first=True
        for p in files:
            with p.open() as src:
                r=csv.reader(src,delimiter='\t'); header=next(r)
                if first: w.writerow(header); first=False
                w.writerows(r)
PY
date -Is > "$OUT/SUCCESS"
echo "Completed: $OUT"
