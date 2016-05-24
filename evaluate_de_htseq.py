"""
Generate intersection comparison tables for each genome from htseq counts
"""
import sys
import os
from collections import defaultdict
ref_dir = 'aligned_reference/'
tgt_dir = 'aligned_new_reference/'

# TODO: this will not work if there is more than one institute/tissue. You need to fix this.


# we need a map of transcript names, because the gencode names are common names
# some transcript names have more than one gencode name. In those cases, we will add counts.
attrs = [x.split() for x in open('/hive/users/ifiddes/ihategit/pipeline/gencode_vm8/C57B6J.tsv')]
attrs_map = {x[0]: x[1] for x in attrs[1:]}

# while we loop over the files, also store the changes in metrics, which are the final lines with underscores
metrics_map = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

ref_data_map = defaultdict(lambda: defaultdict(dict))
for base_path, dirs, files in os.walk(ref_dir):
    if files:
        try:
            genome, institute, tissue, experiment = base_path.replace(ref_dir, "").split("/")
        except ValueError:
            continue
        h = open(os.path.join(base_path, 'Aligned.sortedByCoord.out.bam.htseq.counts'))
        for l in h:
            if l.startswith('_'):
                metric, val = l.split()
                metric = metric[2:]
                metrics_map[genome]['ref'][metric].append(val)
            else:
                gene_id, count = l.split()
                ref_data_map[genome][experiment][gene_id] = count


tgt_data_map = defaultdict(lambda: defaultdict(dict))
for base_path, dirs, files in os.walk(tgt_dir):
    if files:
        try:
            genome, institute, tissue, experiment = base_path.replace(tgt_dir, "").split("/")
        except ValueError:
            continue
        h = open(os.path.join(base_path, 'Aligned.sortedByCoord.out.bam.htseq.counts'))
        for l in h:
            if l.startswith('_'):
                metric, val = l.split()
                metric = metric[2:]
                metrics_map[genome]['ref'][metric].append(val)
            gene_id, count = l.split()
            if ',' in gene_id or 'jg' in gene_id:
                continue  # deal with CGP genes
            gene_name = attrs_map[gene_id]
            exp_count = int(exp_count)
            if gene_id in ref_data_map[genome]:
                tgt_data_map[genome][experiment][gene_name] += exp_count
            else:
                tgt_data_map[genome][experiment][gene_name] = exp_count


# find intersections
shared_genes = defaultdict(dict)
for genome in tgt_data_map:
    for experiment in tgt_data_map[genome]:
        shared_genes[genome][experiment] = tgt_data_map[genome][experiment].viewkeys() & ref_data_map[genome][experiment].viewkeys()


# build database of intersections
matrices = {}
for genome in shared_genes:
    matrix = defaultdict(dict)
    for experiment in shared_genes[genome]:
        ref = ref_data_map[genome][experiment]
        tgt = tgt_data_map[genome][experiment]
        for gene_id in shared_genes[genome][experiment]:
            matrix[experiment + '_REF'][gene_id] = ref[gene_id]
            matrix[experiment + '_TGT'][gene_id] = tgt[gene_id]
    matrices[genome] = matrix


dataframes = {}
for genome, matrix in matrices.iteritems():
    matrix = pd.DataFrame.from_dict(matrix)
    dataframes[genome] = matrix


out_dir = 'matrices_for_deseq'
os.makedirs(out_dir)
for genome, matrix in dataframes.iteritems():
    out_path = os.path.join(out_dir, genome + ".expression_matrix.txt")
    matrix.to_csv(out_path, sep='\t', header=True, index=True)


# R code for DESeq2
#!/usr/bin/env Rscript
library("DESeq2")
args <- commandArgs(TRUE)
data_mat <- data.matrix(read.table(args[1]))
n <- dim(data_mat)[2]
conds <- as.factor(rep(c('REF', 'TGT'), n / 2))
sizes <- MedianNorm(data_mat)
eb_out <- EBTest(Data=data_mat, Conditions=conds, sizeFactors=sizes, maxround=5)
PP <- as.data.frame(GetPPMat(eb_out))
fc_res <- PostFC(eb_out)

results <- cbind(PP, fc_res$PostFC, fc_res$RealFC,unlist(eb_out$C1Mean)[rownames(PP)], unlist(eb_out$C2Mean)[rownames(PP)])
colnames(results) <- c("PPEE", "PPDE", "PostFC", "RealFC","C1Mean","C2Mean")
results <- results[order(results[,"PPDE"], decreasing = TRUE),]
write.table(results, file = args[2], sep = "\t")


from pycbio.sys.procOps import runProc
for genome in dataframes:
    in_path = os.path.join(out_dir, genome + ".expression_matrix.txt")
    out_path = os.path.join(out_dir, genome + ".de_results.txt")
    cmd = ['Rscript', 'run_EBSeq.R', in_path, out_path]
    runProc(cmd)