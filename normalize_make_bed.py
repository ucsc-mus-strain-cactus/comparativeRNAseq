#!/usr/bin/env python
"""
Digs through kallisto output to find median expression for each tissue type.
Correlates this expression with the genome positions in the reference gencode database 
"""

import sys
import os
import argparse
import pysam
import errno
import numpy as np
from itertools import izip
from collections import defaultdict
from math import log

from sonLib.bioio import popenCatch

databases = ["wgEncodeGencodeCompVM4", "wgEncodeGencodeBasicVM4", "wgEncodeGencodePseudoGeneVM4"]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target_dir", help="kallisto output location (base directory)", 
                        default="/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/rnaseq/expression")
    parser.add_argument("--out_dir", help="output dir", default="/cluster/home/ifiddes/mus_strain_data/pipeline_data/rnaseq/log_expression_bedfiles")
    return parser.parse_args()


def get_gene_map():
    db_map = {}
    for db in databases:
        db_map[db] = {}
        result = popenCatch("hgsql -e 'select * from {}' mm10".format(db)).rstrip().split("\n")
        for entry in result[1:]:
            entry = entry.split("\t")
            name, chrom, strand, start, stop, thick_start, thick_stop, block_count, exon_starts, exon_ends = entry[1:11]
            exon_starts = [int(x) for x in exon_starts.split(",") if x != ""]
            exon_ends = [int(x) for x in exon_ends.split(",") if x != ""]
            block_sizes = ",".join(map(str, [e - s for e, s in izip(exon_ends, exon_starts)]))
            block_starts = ",".join(map(str, [x - int(start) for x in exon_starts]))
            db_map[db][name] = [chrom, start, stop, name, 0, strand, thick_start, thick_stop, 0, block_count, 
                                block_sizes, block_starts]
    return db_map


def find_abundance_files(target_dir):
    sample_map = {}
    for base_path, dirs, files in os.walk(target_dir):
        if "abundance.txt" in files:
            try:
                genome, institute, tissue, sample = base_path.replace(target_dir, "").split("/")[1:]
            except ValueError:
                raise RuntimeError("Looks like the directory structure is not correct.")
            if genome not in sample_map:
                sample_map[genome] = {}
            if tissue not in sample_map[genome]:
                sample_map[genome][tissue] = []
            path = os.path.join(target_dir, base_path, "abundance.txt")
            assert os.path.exists(path)
            sample_map[genome][tissue].append(path)
    return sample_map


def munge_abundance_files(abundances):
    gene_map = defaultdict(list)
    for abundance_file in abundances:
        f_h = open(abundance_file)
        _ = f_h.next()  # avoid header
        for line in f_h:
            transcript_id, length, eff_length, est_counts, tpm = line.split()
            gene_map[transcript_id].append(float(tpm))
    median_gene_map = {}
    for transcript, vals in gene_map.iteritems():
        med = np.median(np.array(vals))
        median_gene_map[transcript] = med
    return median_gene_map


def munge_wrapper(sample_map):
    expression_map = {}
    for genome in sample_map:
        if genome not in expression_map:
            expression_map[genome] = {}
        for tissue in sample_map[genome]:
            abundances = sample_map[genome][tissue]
            median_genes = munge_abundance_files(abundances)
            expression_map[genome][tissue] = median_genes
    return expression_map


def rescale_dict(d, new_min=0.0, new_max=1000.0):
    old_min = min(d.itervalues())
    old_max = max(d.itervalues())
    def _rescale(old_min, old_max, new_min, new_max, old_val):
        return int(round((((old_val - old_min) * (new_max - new_min)) / (old_max - old_min)) + new_min))
    return {x: _rescale(old_min, old_max, new_min, new_max, y) for x, y in d.iteritems()}


def fix_dict_sum(d, scale_val):
    return {x: y * scale_val for x, y in d.iteritems()}


def take_log_of_dict(d):
    return {x: log(y + 0.01) for x, y in d.iteritems()}


def normalize_median_genes(median_genes, orig_sum=10 ** 6):
    # first, we expect each sample to sum to 10^6. Rescale if necessary.
    scale_val = orig_sum / sum(median_genes.itervalues())
    if not 0.99 < scale_val < 1.01:
        median_genes = fix_dict_sum(median_genes, scale_val)
    # transform to log scale
    log_transformed = take_log_of_dict(median_genes)
    rescaled_log_transformed = rescale_dict(log_transformed)
    return rescaled_log_transformed


def normalize_expression_map(expression_map):
    normalized_expression_map = {}
    for genome in expression_map:
        if genome not in normalized_expression_map:
            normalized_expression_map[genome] = {}
        for tissue in expression_map[genome]:
            median_genes = expression_map[genome][tissue]
            rescaled_log_transformed = normalize_median_genes(median_genes)
            normalized_expression_map[genome][tissue] = rescaled_log_transformed
    return normalized_expression_map


def filter_normalized_expression_map_by_database(db, normalized_expression_map):
    final_map = {}
    for genome in normalized_expression_map:
        final_map[genome] = {}
        for tissue in normalized_expression_map[genome]:
            if tissue not in final_map[genome]:
                final_map[genome][tissue] = {}
            val_map = normalized_expression_map[genome][tissue]
            genes_to_keep = val_map.viewkeys() & db.viewkeys()
            for gene in genes_to_keep:
                expression = normalized_expression_map[genome][tissue][gene]
                bed_record = db[gene]
                bed_record[4] = expression
                final_map[genome][tissue][gene] = bed_record
    return final_map


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise


def main():
    args = parse_args()
    mkdir_p(args.out_dir)
    db_map = get_gene_map()
    sample_map = find_abundance_files(args.target_dir)
    expression_map = munge_wrapper(sample_map)
    normalized_expression_map = normalize_expression_map(expression_map)
    final_results = {}
    for db, record_map in db_map.iteritems():
        final_results[db] = filter_normalized_expression_map_by_database(record_map, normalized_expression_map)
    for db in final_results:
        for genome in final_results[db]:
            for tissue in final_results[db][genome]:
                out_path = os.path.join(args.out_dir, db, genome, tissue)
                mkdir_p(out_path)
                with open(os.path.join(out_path, "expression.bed"), "w") as outf:
                    for gene, record in final_results[db][genome][tissue].iteritems():
                        outf.write("\t".join(map(str, record)) + "\n")


if __name__ == "__main__":
    main()
