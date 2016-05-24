"""
Calculates plots
"""
import os
import sys
sys.path.extend(['/Users/ifiddes/hive/ihategit/pipeline/', '/Users/ifiddes/hive/ihategit/pipeline/submodules', '/Users/ifiddes/hive/ihategit/pipeline/submodules/pycbio', '/Users/ifiddes/hive/ihategit/pipeline/submodules/comparativeAnnotator'])
import argparse
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from math import log
from collections import Counter, defaultdict
from pycbio.sys.fileOps import iterRows, ensureDir, iterLines, ensureFileDir
from pycbio.sys.procOps import runProc, callProcLines


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--attrs', default="/Users/ifiddes/hive/ihategit/pipeline/gencode_vm8/C57B6J.tsv")
    parser.add_argument('--counts_dir', default='counts_tables_deseq')
    parser.add_argument('--deseq_dir', default='deseq_results')
    parser.add_argument('--genomes', nargs='+', default=None)
    parser.add_argument('--ref_dir', default='star_alignments_v4/mm10')
    parser.add_argument('--tgt_dir', default='star_alignments_v4/strain_specific')
    parser.add_argument('--out_dir', default='scatterplots')
    parser.add_argument('--experiment', required=True)
    return parser.parse_args()


def get_tx_map(attrs):
    tx_map = {}
    for x in iterRows(attrs, skipLines=1):
        tx_map[x[3]] = x[0]
    return tx_map


def get_common_name_map(attrs):
    common_name_map = {}
    for x in iterRows(attrs, skipLines=1):
        common_name_map[x[1]] = x[0]
    return common_name_map


def load_dir(data, mode, path, genomes, common_name_map):
    for g in genomes:
        for base_dir, d, files in os.walk(os.path.join(path, g)):
            if files:
                files = [os.path.join(base_dir, x) for x in files if x.endswith('.filtered.bam.counts.cds')]
                for p in files:
                    for l in iterLines(p, skipLines=2):
                        gene_name, chrom, start, end, strand, length, count = l.split()
                        if mode == 'ref':
                            gene_name = common_name_map[gene_name]
                        data[g][mode][gene_name].append(count)


def filter_data(d):
    ref_txs = {ens_id for ens_id, val in d['ref'].iteritems()}
    tgt_txs = {ens_id for ens_id, val in d['tgt'].iteritems()}
    keep_names = ref_txs & tgt_txs
    r = {}
    for n in keep_names:
        r[n] = [d['ref'][n], d['tgt'][n]]
    return r


def construct_counts_tables(data, base_path):
    for g, d in data.iteritems():
        counts_path = os.path.join(base_path, g + '.counts.txt')
        conds_path = os.path.join(base_path, g + '.conds.txt')
        ensureFileDir(counts_path)
        cleaned_d = filter_data(d)
        with open(counts_path, 'w') as counts:
            first_line = cleaned_d.values()[0]
            assert len(first_line[0]) == len(first_line[1])
            header = ['ref{}'.format(i) for i in range(len(first_line[0]))] + ['tgt{}'.format(i) for i in range(len(first_line[0]))]
            counts.write('\t'.join(header) + "\n")
            for ens_id, vals in cleaned_d.iteritems():
                f = [item for sublist in vals for item in sublist]
                l = [ens_id] + f
                counts.write('\t'.join(l) + '\n')
        with open(conds_path, 'w') as outf:
            outf.write('condition\n')
            for x in header:
                outf.write('\t'.join([x, ''.join([q for q in x if q.isalpha()])]) + "\n")



def strip_version(df):
    df.index = [x.split('.')[0] for x in df.index]
    return df


def load_paralogy(g, df, para_dir='paralogy_filtered_sets'):
    unique_genes = [x.rstrip() for x in open(os.path.join(para_dir, g + '.txt'))]
    gene_r = {gene_id: 'Not Unique' for gene_id in df.index}
    for gene_id in unique_genes:
        gene_r[gene_id] = 'Unique'
    return df.join(pd.Series(gene_r, name='Paralogy'), how='inner')


def calculate_pval_vector(df, ref_col_name, tgt_col_name, col='Unnamed: 2'):
    num_sig = len([x for x in df[col] if x <= 0.05])
    num_non_sig = len(df) - num_sig
    df['ref'] = map(log, df['ref'])
    df['tgt'] = map(log, df['tgt'])
    p_val_vector = ['<= 0.05 (n = {})'.format(num_sig) if x <= 0.05 else '> 0.05 (n = {})'.format(num_non_sig) for x in df[col]]
    df[col] = p_val_vector
    df.columns = [ref_col_name, tgt_col_name, 'Adjusted p-value']
    return df


def generate_plots(genomes, deseq_path, out_dir):
    for g in genomes:
        base_out = os.path.join(out_dir, g)
        ensureDir(base_out)
        main_plot = os.path.join(base_out, 'nonzero_log.pdf')
        para_plot = os.path.join(base_out, 'nonzero_log_paralogy.pdf')
        ref_col_name = 'Normalized counts (aligned to mm10)'
        tgt_col_name = 'Normalized counts (aligned to {})'.format(g)
        df = pd.read_table(os.path.join(deseq_path, g + '.data_table.tsv'), sep=' ', header=0, index_col=0)
        df = df[df['ref'] >= 1]
        df = df[df['tgt'] >= 1]
        df = strip_version(df)
        df = calculate_pval_vector(df, ref_col_name, tgt_col_name)
        p = sns.lmplot(x=ref_col_name, y=tgt_col_name, data=df, scatter_kws={"s": 6, "alpha": 0.7}, hue='Adjusted p-value', fit_reg=False)
        p.set(xlim=(0, 10), ylim=(0, 10))
        p.savefig(main_plot, format='pdf')
        df2 = load_paralogy(g, df)
        p2 = sns.lmplot(x=ref_col_name, y=tgt_col_name, data=df2, scatter_kws={"s": 6, "alpha": 0.7}, col='Paralogy', hue='Adjusted p-value', fit_reg=False)
        p2.set(xlim=(0, 10), ylim=(0, 10))
        p2.savefig(para_plot, format='pdf')
        plt.close('all')


def main():
    args = parse_args()
    data = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    ref_path = os.path.join(args.ref_dir, args.experiment)
    tgt_path = os.path.join(args.tgt_dir, args.experiment)
    if args.genomes is None:
        args.genomes = os.listdir(ref_path)
    common_name_map = get_common_name_map(args.attrs)
    load_dir(data, 'ref', ref_path, args.genomes, common_name_map)
    load_dir(data, 'tgt', tgt_path, args.genomes, common_name_map)
    base_counts_dir = os.path.join(args.counts_dir, args.experiment)
    construct_counts_tables(data, base_counts_dir)
    deseq_path = os.path.join(args.deseq_dir, args.experiment)
    ensureDir(deseq_path)
    for g in args.genomes:
        runProc(['Rscript', 'run_DEseq.R', g, deseq_path, base_counts_dir])
    out_dir = os.path.join(args.out_dir, args.experiment)
    generate_plots(args.genomes, deseq_path, out_dir)


if __name__ == '__main__':
    main()
