#!/usr/bin/env python
"""
Kallisto from fastqs.

Hunts through --dir (which is a FTP download from ftp://ftp-mouse.sanger.ac.uk/) and runs kallisto. 
Designed with 1505 fastq release in mind. The starting directory structure will be copied.

Attempts to determine if the input RNAseq was paired or unpaired and acts accordingly.
"""

import sys
import os
import argparse
import pysam
from itertools import izip

from src.helperFunctions import mkdir_p, fastq_read_size, find_paired_fastqs

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import setLoggingFromOptions, system, logger, fastqRead

def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source_dir", help="Directory containing files.",
                        default="/hive/groups/recon/projs/mus_strain_cactus/data/assembly_rel_1505/fastqs/"
                                "ftp-mouse.sanger.ac.uk/REL-1505-RNA-Seq/fastq")
    parser.add_argument("--reference", help="kallisto reference", 
                        default="gencode_VM4_Reference/kallisto_reference/tgencode.VM4.idx")
    parser.add_argument("--out_dir", help="output location (base directory)", 
                        default="/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/rnaseq/kallisto_expression")
    return parser


def wrapper(target, source_dir, reference, out_dir):
    for base_path, dirs, files in os.walk(source_dir):
        if files:
            try:
                genome, rna_seq, institute, tissue = base_path.replace(source_dir, "").split("/")[1:]
            except ValueError:
                raise RuntimeError("Looks like the directory structure is not what we expected.")
            fastq_files = find_paired_fastqs(source_dir, base_path, files)
            for experiment, fastq_path in fastq_files.iteritems():
                try:
                    fwd_fastq_path, rev_fastq_path = fastq_path
                    target.addChildTargetFn(kallisto_paired, args=(genome, institute, tissue, reference, out_dir, 
                                                                   experiment, fwd_fastq_path, rev_fastq_path))
                except ValueError:
                    read_size = fastq_read_size(fastq_path)
                    target.addChildTargetFn(kallisto_single, args=(genome, institute, tissue, reference, out_dir, 
                                                                   experiment, fastq_path, read_size))


def kallisto_paired(target, genome, institute, tissue, reference, out_dir, experiment, fwd_fastq_path, rev_fastq_path):
    out_path = build_out_dirs(out_dir, genome, institute, tissue, experiment)
    system("kallisto quant -i {} -o {} {} {}".format(reference, out_path, fwd_fastq_path, rev_fastq_path))


def kallisto_single(target, genome, institute, tissue, reference, out_dir, experiment, fastq_path, read_size):
    out_path = build_out_dirs(out_dir, genome, institute, tissue, experiment)
    system("kallisto quant --single -l {} -i {} -o {} {}".format(read_size, reference, out_path, fastq_path))


def build_out_dirs(out_dir, genome, institute, tissue, experiment):
    out_path = os.path.join(out_dir, genome, institute, tissue, experiment)
    mkdir_p(out_path)
    return out_path


def main():
    parser = build_parser()
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    setLoggingFromOptions(args)

    i = Stack(Target.makeTargetFn(wrapper, args=(args.source_dir, args.reference, args.out_dir))).startJobTree(args)

    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == "__main__":
    from bin.kallisto_from_fastq import *
    main()
