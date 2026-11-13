#!/usr/bin/env bash
#SBATCH --job-name=triticum_ITS2_ASV
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --partition=smp
#SBATCH --output=triticum_ITS2_%j.out
#SBATCH --error=triticum_ITS2_%j.err
# Usage: bash 13_triticum_ITS2_ASV.sh ASVs.fasta [abundances.tsv]
# Activate an environment containing blastn, makeblastdb, bowtie2-inspect, python3.
# Initialisation de Conda, comme dans le script shotgun.
set -Eeo pipefail

module load conda/4.12.0
source ~/.bashrc
conda activate "${CONDA_ENV:-mapdamage_py39}"

set -Eeuo pipefail

echo "CONDA_DEFAULT_ENV=${CONDA_DEFAULT_ENV:-non_defini}"
echo "CONDA_PREFIX=${CONDA_PREFIX:-non_defini}"
export LC_ALL=C
trap 'printf "ERROR line %s: %s\n" "$LINENO" "$BASH_COMMAND" >&2' ERR
ASVS=${1:?Usage: script ASVs.fasta [abundances.tsv]}
TABLE=${2:-}
INDEX=${INDEX:-/home/plstenge/adna_analysis_new-main/Bowtie2_index/ITS_ABDGH}
THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-4}}
MIN_ID=${MIN_ID:-98}
MIN_COV=${MIN_COV:-95}
MIN_SCORE_GAP=${MIN_SCORE_GAP:-5}
OUT_BASE=${OUT_BASE:-$PWD/13_triticum_ITS2_ASV/runs}
OUT="$OUT_BASE/run_$(date +%Y%m%d_%H%M%S)_${SLURM_JOB_ID:-manual_$$}"
[[ -s "$ASVS" ]] || { echo "Missing ASV FASTA" >&2; exit 1; }
[[ -z "$TABLE" || -s "$TABLE" ]] || { echo "Missing abundance TSV" >&2; exit 1; }
for cmd in python3 bowtie2-inspect blastn makeblastdb; do
    command -v "$cmd" >/dev/null || { echo "Missing executable: $cmd" >&2; exit 1; }
done
mkdir -p "$OUT"/{qc,reference,results}
exec > >(tee -a "$OUT/pipeline.log") 2>&1
cp -- "$0" "$OUT/qc/submitted_script.sh"
printf 'Output: %s\nIndex: %s\n' "$OUT" "$INDEX"
bowtie2-inspect "$INDEX" > "$OUT/reference/index_sequences.fasta" 2> "$OUT/qc/inspect.log"
if [[ -n "${CONTIG_MAP:-}" ]]; then
    cp -- "$CONTIG_MAP" "$OUT/qc/contig_subgenome.tsv"
else
    printf 'contig\tsubgenome\nrDNA_6BS_Langdon_extraction_3\tB\nrDNA_5DS_tauschi_extraction_3\tD\nMZ547657_extraction_2\tG\nrDNA_1BS_Langdon_extraction_3\tB\nMZ547658_extraction_2\tA\nHordeum\tHord\n' > "$OUT/qc/contig_subgenome.tsv"
fi
cat > "$OUT/qc/worker.py" <<'PY'
import csv, sys, math, hashlib
from pathlib import Path
from collections import defaultdict

def fasta(path):
    records = {}; name = None; seq = []
    def save():
        if name is None: return
        s = ''.join(seq).upper()
        if not s or set(s)-set('ACGTRYSWKMBDHVN'):
            raise ValueError('Empty/invalid nucleotide sequence: '+name)
        if name in records: raise ValueError('Duplicate FASTA ID: '+name)
        records[name] = s
    with open(path) as h:
        for line in h:
            line = line.strip()
            if not line: continue
            if line.startswith('>'):
                save(); name = line[1:].split()[0]; seq = []
            else:
                if name is None: raise ValueError('FASTA sequence without header')
                seq.append(line)
    save()
    if not records: raise ValueError('Empty FASTA')
    return records

def write_tsv(path, header, rows):
    with open(path, 'w', newline='') as h:
        w = csv.writer(h, delimiter='\t'); w.writerow(header); w.writerows(rows)

mode, root = sys.argv[1], Path(sys.argv[2])
if mode == 'prepare':
    asvs = fasta(sys.argv[3]); refs = fasta(root/'reference/index_sequences.fasta')
    with open(root/'qc/contig_subgenome.tsv') as h:
        rows = list(csv.DictReader(h, delimiter='\t'))
    mapping = {r['contig']: r['subgenome'] for r in rows}
    if len(mapping) != len(rows) or any(not k or not v for k,v in mapping.items()):
        raise ValueError('Invalid/duplicate contig map')
    missing = set(refs)-set(mapping)
    if missing: raise ValueError('Unmapped reference IDs; supply CONTIG_MAP: '+', '.join(sorted(missing)))
    for filename, records in [('asvs.fasta', asvs), ('references.fasta', refs)]:
        with open(root/'reference'/filename, 'w') as h:
            for i,(k,v) in enumerate(records.items(), 1):
                ident = k if filename == 'asvs.fasta' else f'REF{i:06d}'
                h.write(f'>{ident}\n{v}\n')
    write_tsv(root/'reference/reference_ids.tsv', ['blast_id','contig','subgenome'],
              [(f'REF{i:06d}', k, mapping[k]) for i,k in enumerate(refs,1)])
    write_tsv(root/'qc/asv_lengths.tsv', ['ASV_ID','length'], [(k,len(v)) for k,v in asvs.items()])
    with open(root/'qc/reference_count.txt','w') as h: h.write(str(len(refs)))
    write_tsv(root/'qc/input_sha256.tsv', ['file','sha256'],
              [(str(p),hashlib.sha256(Path(p).read_bytes()).hexdigest()) for p in sys.argv[3:] if p])
elif mode == 'summarize':
    min_id, min_cov, gap = map(float, sys.argv[3:6]); table = sys.argv[6]
    if not (0 <= min_id <= 100 and 0 <= min_cov <= 100 and math.isfinite(gap) and gap >= 0):
        raise ValueError('Invalid thresholds')
    asvs = fasta(root/'reference/asvs.fasta')
    with open(root/'reference/reference_ids.tsv') as h:
        refs = {r['blast_id']: (r['contig'],r['subgenome']) for r in csv.DictReader(h,delimiter='\t')}
    hits = defaultdict(list); hitrows = []
    with open(root/'results/blast.raw.tsv') as h:
        for line in h:
            f = line.rstrip('\n').split('\t')
            q, ref = f[:2]; contig, group = refs[ref]
            pid, length, qlen, qs, qe, ss, se, ev, bits, score, mism, gaps = map(float, f[2:])
            cov = 100*(abs(qe-qs)+1)/qlen
            passed = pid >= min_id and cov >= min_cov
            hit = dict(contig=contig, group=group, pid=pid, cov=cov, score=score, bits=bits)
            hits[q].append((passed, hit))
            hitrows.append([q,contig,group,*f[2:],cov,int(passed)])
    write_tsv(root/'results/all_hits.tsv',
        ['ASV_ID','contig','subgenome','pident','alignment_length','query_length','qstart','qend',
         'sstart','send','evalue','bitscore','raw_score','mismatches','gap_openings','query_coverage','pass_thresholds'], hitrows)
    assignments = {}; rows = []; categories = set()
    for q in asvs:
        eligible = [v for ok,v in hits[q] if ok]
        if not hits[q]: status, category, best, candidates, delta = 'NO_HIT','NO_HIT',None,[], 'NA'
        elif not eligible: status, category, best, candidates, delta = 'LOW_QUALITY','LOW_QUALITY',None,[], 'NA'
        else:
            bygroup = {}
            for v in eligible:
                g = v['group']
                if g not in bygroup or v['score'] > bygroup[g]['score']: bygroup[g] = v
            ranked = sorted(bygroup.values(),key=lambda v:(-v['score'],v['group']))
            best = ranked[0]
            # Equal scores are ALWAYS ambiguous, even if MIN_SCORE_GAP=0.
            candidates = sorted(v['group'] for v in ranked
                                if v['score'] == best['score'] or best['score']-v['score'] < gap)
            delta = best['score']-ranked[1]['score'] if len(ranked)>1 else 'NA'
            status = 'UNIQUE_PANEL' if len(candidates)==1 else 'AMBIGUOUS'
            category = candidates[0] if status=='UNIQUE_PANEL' else 'AMBIGUOUS_'+'_'.join(candidates)
        assignments[q] = category; categories.add(category)
        rows.append([q,len(asvs[q]),status,category,','.join(candidates) or 'NA',
                     best['contig'] if best else 'NA',best['pid'] if best else 'NA',
                     best['cov'] if best else 'NA',best['score'] if best else 'NA',delta])
    write_tsv(root/'results/asv_assignments.tsv',
              ['ASV_ID','length','status','category','candidate_groups','representative_best_contig',
               'best_pident','best_query_coverage','best_raw_score','score_gap_next_group'], rows)
    all_categories = sorted(categories | {g for c,g in refs.values()} | {'NO_HIT','LOW_QUALITY'})
    write_tsv(root/'results/asv_summary.tsv', ['category','n_ASVs'],
              [(c,sum(v==c for v in assignments.values())) for c in all_categories])
    if table:
        with open(table,encoding='utf-8-sig',newline='') as h:
            reader = csv.reader(h,delimiter='\t')
            header = next(reader)
            if len(header)<2 or header[0] not in ('#ASV_ID','ASV_ID','Feature ID','feature-id'):
                raise ValueError('Expected tab-separated table with ASV_ID or #ASV_ID as first column')
            sample_idx = [i for i in range(1,len(header)) if header[i].lower() not in ('taxonomy','taxon')]
            if not sample_idx or len(header)!=len(set(header)): raise ValueError('Invalid/duplicate sample columns')
            abundance = {}; ignored = 0
            for row in reader:
                if not row: continue
                if len(row)!=len(header): raise ValueError('Wrong number of TSV columns for '+row[0])
                if row[0] not in asvs: ignored += 1; continue
                if row[0] in abundance: raise ValueError('Duplicate ASV table row: '+row[0])
                values=[]
                for i in sample_idx:
                    try: x=int(row[i])
                    except ValueError: raise ValueError('Non-integer count: '+row[0]+'/'+header[i])
                    if x<0: raise ValueError('Negative abundance')
                    values.append(x)
                abundance[row[0]] = values
        missing = set(asvs)-set(abundance)
        if missing: raise ValueError('ASVs missing from count table: '+', '.join(sorted(missing)))
        weighted=[]
        for j,i in enumerate(sample_idx):
            total = sum(v[j] for v in abundance.values())
            for c in all_categories:
                n = sum(abundance[q][j] for q in asvs if assignments[q]==c)
                weighted.append([header[i],c,n,total,100*n/total if total else 'NA'])
        write_tsv(root/'results/subgenome_summary.tsv',
                  ['sample','category','ASV_read_count','total_selected_ASV_reads','percent_selected_ASV_reads'],weighted)
        write_tsv(root/'qc/table_qc.tsv',['metric','value'],
                  [('selected_ASVs',len(abundance)),('other_table_rows_ignored',ignored)])
else: raise ValueError(mode)
PY
python3 "$OUT/qc/worker.py" prepare "$OUT" "$ASVS" "$TABLE"
NREF=$(cat "$OUT/qc/reference_count.txt")
makeblastdb -in "$OUT/reference/references.fasta" -dbtype nucl -parse_seqids \
    -out "$OUT/reference/ITS" > "$OUT/qc/makeblastdb.log" 2>&1
blastn -query "$OUT/reference/asvs.fasta" -db "$OUT/reference/ITS" \
    -task blastn -word_size 7 -dust no -soft_masking false -strand both \
    -reward 2 -penalty -3 -gapopen 5 -gapextend 2 -evalue 1e-10 \
    -max_target_seqs "$NREF" -max_hsps 1 -num_threads "$THREADS" \
    -outfmt '6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore score mismatch gapopen' \
    -out "$OUT/results/blast.raw.tsv" 2> "$OUT/qc/blast.log"
python3 "$OUT/qc/worker.py" summarize "$OUT" "$MIN_ID" "$MIN_COV" "$MIN_SCORE_GAP" "$TABLE"
printf 'parameter\tvalue\nindex\t%s\nmin_identity\t%s\nmin_query_coverage\t%s\nmin_raw_score_gap\t%s\nthreads\t%s\n' \
    "$INDEX" "$MIN_ID" "$MIN_COV" "$MIN_SCORE_GAP" "$THREADS" > "$OUT/qc/parameters.tsv"
blastn -version > "$OUT/qc/blast_version.txt"
bowtie2-inspect --version > "$OUT/qc/bowtie2_version.txt"
date -Is > "$OUT/SUCCESS"
printf 'Completed: %s\n' "$OUT"
