# mouse notes

# find . -mindepth 1 -maxdepth 1 -type d | xargs -n 1 -P 10 -I {} sh -c 'n={}; g=${n:2}; cd ${g} && STAR --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles ${g}.fa --sjdbGTFfile ${g}.gtf'
# python bin/STAR_from_fastq.py  --reference STAR_references/C57B6J/ --out_dir aligned_reference  --jobTree align_jobTree_ref --batchSystem parasol --parasolCommand ~/bin/remparasol --rsem --include_tissues liver &> log
# find aligned_reference -type f -name "*Coord*bam" | xargs -n 1 -P 20 -I{} sh -c  "htseq-count -m union -r pos -i gene_id -a 10 --stranded=no -f bam '{}' /hive/users/ifiddes/ihategit/pipeline/gencode_vm8/C57B6J.gtf > '{}'.htseq.counts"

# python bin/STAR_from_fastq.py  --ref_dir STAR_references_v4/ --out_dir aligned_new_reference  --jobTree align_jobTree_v4 --batchSystem parasol --parasolCommand ~/bin/remparasol --htseq --include_tissues liver &> log
# find aligned_new_reference -type f -name "*Coord*bam" | xargs -n 1 -P 20 -I{} sh -c  'n="{}"; g=`dirname $n | cut -d "/" -f 2`; ref=STAR_references/${g}/${g}.gtf; htseq-count -m union -r pos -i gene_id -a 10 --stranded=no -f bam "{}" ${ref} > "{}".htseq.counts'


# load a transcript-gene map
from pycbio.bio.bio import *
from pycbio.bio.transcripts import *
from pycbio.sys.fileOps import *
tx_map = {}
recs = [x.split() for x in open("/hive/users/ifiddes/ihategit/pipeline/gencode_vm8/C57B6J.tsv")]
recs = recs[1:]
for x in recs:
    tx_map[x[3]] = x[0]


# also load a common name map, which is what we will operate on
common_name_map = {}
for x in recs:
    common_name_map[x[1]] = x[0]

# for each genome, load the counts seen in both reference and target

def load_dir(data, mode, path, genomes):
    for g in genomes:
        for base_dir, d, files in os.walk(os.path.join(path, g)):
            if files:
                files = [os.path.join(base_dir, x) for x in files if x.endswith('cds')]
                assert len(files) > 0
                for p in files:
                    for l in iterLines(p, skipLines=2):
                        gene_name, chrom, start, end, strand, length, count = l.split()
                        if mode == 'ref':
                            gene_name = common_name_map[gene_name]
                        data[g][mode][gene_name].append(count)



#genomes = ['PWK_PhJ', 'NOD_ShiLtJ', 'A_J', 'WSB_EiJ', 'NZO_HlLtJ', '129S1_SvImJ', 'CAST_EiJ']
genomes = ['C57BL_6NJ']
base_dir = '/hive/users/ifiddes/comparativeRNAseq/star_alignments_v4'
ref_dir = 'mm10'
tgt_dir = 'strain_specific'
experiment = 'brain_sanger'
data = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
ref_path = os.path.join(base_dir, ref_dir, experiment)
tgt_path = os.path.join(base_dir, tgt_dir, experiment)
load_dir(data, 'ref', ref_path, genomes)
load_dir(data, 'tgt', tgt_path, genomes)

# generate table for R
# for this, we need two tsv files, one with counts and one with conditions

def filter_data(d):
    ref_txs = {ens_id for ens_id, val in d['ref'].iteritems()}
    tgt_txs = {ens_id for ens_id, val in d['tgt'].iteritems()}
    keep_names = ref_txs & tgt_txs
    r = {}
    for n in keep_names:
        r[n] = [d['ref'][n], d['tgt'][n]]
    return r


for g, d in data.iteritems():
    cleaned_d = filter_data(d)
    with open('deseq_counts_tables/{}.counts.txt'.format(g), 'w') as counts:
        first_line = cleaned_d.values()[0]
        assert len(first_line[0]) == len(first_line[1])
        header = ['ref{}'.format(i) for i in range(len(first_line[0]))] + ['tgt{}'.format(i) for i in range(len(first_line[0]))]
        counts.write('\t'.join(header) + "\n")
        for ens_id, vals in cleaned_d.iteritems():
            f = [item for sublist in vals for item in sublist]
            l = [ens_id] + f
            counts.write('\t'.join(l) + '\n')
    with open('deseq_counts_tables/{}.conds.txt'.format(g), 'w') as outf:
        outf.write('condition\n')
        for x in header:
            outf.write('\t'.join([x, ''.join([q for q in x if q.isalpha()])]) + "\n")


# R code below to run DESeq, saving a normalized table to plot.

#!/usr/bin/env Rscript
library("DESeq2")
args <- commandArgs(TRUE)
genome <- args[1]
out_path <- args[2]
data_path <- paste0('~/hive/comparativeRNAseq/deseq_counts_tables/', genome, '.counts.txt')
conds_path <- paste0('~/hive/comparativeRNAseq/deseq_counts_tables/', genome, '.conds.txt')
conds <- read.table(conds_path, sep="\t", header=T, row.names=1)
data <- read.table(data_path, sep="\t", header=T, row.names=1)
dds <- DESeqDataSetFromMatrix(countData=data, colData=conds, design=~condition)
dds$condition <- relevel(dds$condition, ref='ref')
dds <- DESeq(dds)
res <- results(dds)
pdf(paste0(out_path, '/', genome, '.maplot.pdf'))
plotMA(res, main=genome)
dev.off()
p <- paste0(out_path, '/', genome, '.results.rData')
save(dds, file=p)
baseMeanPerLvl <- sapply( levels(dds$condition), function(lvl) rowMeans( counts(dds,normalized=TRUE)[,dds$condition == lvl] ) )
p2 <- paste0(out_path, '/', genome, '.data_table.tsv')
q <- cbind(baseMeanPerLvl, res$padj)
write.table(q, p2)


for g in genomes:
    runProc(['Rscript', 'run_DEseq.R', g, 'deseq_results'])

# scatterplots

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from math import log

def strip_version(df):
    df.index = [x.split('.')[0] for x in df.index]
    return df


def load_paralogy(g, df, para_dir='paralogy_filtered_sets'):
    unique_genes = [x.rstrip() for x in open(os.path.join(para_dir, g + '.txt'))]
    gene_r = {gene_id: 'Not Unique' for gene_id in df.index}
    for gene_id in unique_genes:
        gene_r[gene_id] = 'Unique'
    return df.join(pd.Series(gene_r, name='Paralogy'), how='inner')


def calculate_pval_vector(df, col='Unnamed: 2'):
    num_sig = len([x for x in df[col] if x <= 0.05])
    num_non_sig = len(df) - num_sig
    df['ref'] = map(log, df['ref'])
    df['tgt'] = map(log, df['tgt'])
    p_val_vector = ['<= 0.05 (n = {})'.format(num_sig) if x <= 0.05 else '> 0.05 (n = {})'.format(num_non_sig) for x in df[col]]
    df[col] = p_val_vector
    df.columns = [ref_col_name, tgt_col_name, 'Adjusted p-value']
    return df


for g in genomes:
    ref_col_name = 'Normalized counts (aligned to mm10)'
    tgt_col_name = 'Normalized counts (aligned to {})'.format(g)
    df = pd.read_table(os.path.join('deseq_results', g + '.data_table.tsv'), sep=' ', header=0, index_col=0)
    df = df[df['ref'] >= 1]
    df = df[df['tgt'] >= 1]
    df = strip_version(df)
    df = calculate_pval_vector(df)
    p = sns.lmplot(x=ref_col_name, y=tgt_col_name, data=df, scatter_kws={"s": 6, "alpha": 0.7}, hue='Adjusted p-value', fit_reg=False)
    p.set(xlim=(0, 10), ylim=(0, 10))
    p.savefig('scatterplots/{}_nonzero_log.pdf'.format(g), format='pdf')
    df2 = load_paralogy(g, df)
    p2 = sns.lmplot(x=ref_col_name, y=tgt_col_name, data=df2, scatter_kws={"s": 6, "alpha": 0.7}, col='Paralogy', hue='Adjusted p-value', fit_reg=False)
    p2.set(xlim=(0, 10), ylim=(0, 10))
    p2.savefig('scatterplots/{}_nonzero_log_biotype_paralogy.pdf'.format(g), format='pdf')