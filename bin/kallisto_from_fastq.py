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

from src.helperFunctions import mkdir_p, fastq_read_size, find_paired_fastqs, name_map

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import setLoggingFromOptions, system, logger, fastqRead

def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source_dir", help="Directory containing files.",
                        default="/hive/groups/recon/projs/mus_strain_cactus/data/assembly_rel_1505/fastqs/"
                                "ftp-mouse.sanger.ac.uk/REL-1505-RNA-Seq/fastq")
    parser.add_argument("--out_dir", help="output location (base directory)", required=True)
    parser.add_argument('--include_tissues', nargs='+', default=None,
                        help='Space separated list of tissues to include. If set, will only do these.')
    refs = parser.add_mutually_exclusive_group(required=True)
    refs.add_argument("--reference", help="reference directory for STAR")
    refs.add_argument('--ref_dir', help='reference directory of directories. Mutually exclusive with --reference. '
                                        'Each directory should be in the format genome/ref_files')
    return parser


def wrapper(target, args):
    for base_path, dirs, files in os.walk(args.source_dir):
        if files:
            try:
                genome, rna_seq, institute, tissue = base_path.replace(args.source_dir, "").split("/")[1:]
            except ValueError:
                raise RuntimeError("Looks like the directory structure is not what we expected.")
            if args.include_tissues is not None and tissue not in args.include_tissues:
                continue
            fastq_files = find_paired_fastqs(args.source_dir, base_path, files)
            if args.reference is not None:
                reference = args.reference
            else:
                reference = os.path.join(args.ref_dir, genome)
            assert reference is not None
            for experiment, fastq_path in fastq_files.iteritems():
                try:
                    fwd_fastq_path, rev_fastq_path = fastq_path
                    target.addChildTargetFn(kallisto_paired, args=(genome, institute, tissue, reference, args.out_dir,
                                                                   experiment, fwd_fastq_path, rev_fastq_path))
                except ValueError:
                    read_size = fastq_read_size(fastq_path)
                    target.addChildTargetFn(kallisto_single, args=(genome, institute, tissue, reference, args.out_dir,
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

    i = Stack(Target.makeTargetFn(wrapper, args=(args,))).startJobTree(args)

    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == "__main__":
    from bin.kallisto_from_fastq import *
    main()
