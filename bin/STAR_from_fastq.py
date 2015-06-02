#!/usr/bin/env python
"""
Aligns fastq to transcriptomes using STAR.

Hunts through --dir (which is a FTP download from ftp://ftp-mouse.sanger.ac.uk/) and runs STAR aligning to the
transcriptome. Designed with 1505 fastq release in mind. The starting directory structure will be copied.

Attempts to determine if the input RNAseq was paired or unpaired and acts accordingly.
"""

import sys
import os
import argparse

from src.helperFunctions import mkdir_p, fastq_read_size, find_paired_fastqs, name_map

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import setLoggingFromOptions, system, logger, getRandomAlphaNumericString

def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source_dir", help="Directory containing fastq files.",
                        default="/hive/groups/recon/projs/mus_strain_cactus/data/assembly_rel_1505/fastqs/"
                                "ftp-mouse.sanger.ac.uk/REL-1505-RNA-Seq/fastq")
    parser.add_argument("--out_dir", help="output location (base directory)", 
                        default="/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/rnaseq/STAR_output")
    parser.add_argument("--reference", help="reference directory for STAR",
                        default="gencode_VM4_Reference/STAR_reference/")
    return parser


# star_flags variable should not need changing per run
star_flags = (" --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate  --outBAMcompression -1 "
             "--outFilterType BySJout --outSAMattributes NH HI AS NM MD --outFilterMultimapNmax 20 "
             "--outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 "
             "--alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 "
             "--sjdbScore 1 --readFilesCommand zcat --quantMode TranscriptomeSAM --quantTranscriptomeBan Singleend "
             "--outWigType wiggle --limitBAMsortRAM 80000000000")

# format star_base_cmd for each STAR call
star_paired_cmd = "STAR --genomeDir {} --readFilesIn {} {} --outFileNamePrefix {} --outTmpDir {}"
star_single_cmd = "STAR --genomeDir {} --readFilesIn {} --outFileNamePrefix {} --outTmpDir {}"


def wrapper(target, source_dir, reference, out_dir):
    for base_path, dirs, files in os.walk(source_dir):
        if files:
            try:
                genome, rna_seq, institute, tissue = base_path.replace(source_dir, "").split("/")[1:]
            except ValueError:
                raise RuntimeError("Looks like the directory structure is not what we expected.")
            fastq_files = find_paired_fastqs(source_dir, base_path, files)
            genome = name_map[genome]
            for experiment, fastq_path in fastq_files.iteritems():
                try:
                    fwd_fastq_path, rev_fastq_path = fastq_path
                    target.addChildTargetFn(run_paired_star, args=(genome, institute, tissue, reference, out_dir, 
                                                                   experiment, fwd_fastq_path, rev_fastq_path))
                except ValueError:
                    target.addChildTargetFn(run_single_star, args=(genome, institute, tissue, reference, out_dir, 
                                                                   experiment, fastq_path))


def run_paired_star(target, genome, institute, tissue, reference, out_dir, experiment, fwd_fastq_path, rev_fastq_path):
    tmp_dir = os.path.join(target.getLocalTempDir(), "tmp_" + getRandomAlphaNumericString())
    out_path = build_out_dirs(out_dir, genome, institute, tissue, experiment) + "/"
    this_star_base_cmd = star_paired_cmd.format(reference, fwd_fastq_path, rev_fastq_path, out_path, tmp_dir)
    system(this_star_base_cmd + star_flags)


def run_single_star(target, genome, institute, tissue, reference, out_dir, experiment, fastq_path):
    tmp_dir = os.path.join(target.getLocalTempDir(), "tmp_" + getRandomAlphaNumericString())
    out_path = build_out_dirs(out_dir, genome, institute, tissue, experiment) + "/"
    this_star_base_cmd = star_single_cmd.format(reference, fastq_path, out_path, tmp_dir)
    system(this_star_base_cmd + star_flags)


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
    from bin.STAR_from_fastq import *
    main()