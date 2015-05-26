#!/usr/bin/env python
"""
Kallisto from aligned BAMs.

Hunts through --dir (which is a FTP download from ftp://ftp-mouse.sanger.ac.uk/) and runs kallisto. 
Designed with 1505 BAM release in mind. The starting directory structure will be copied.

Attempts to determine if the input RNAseq was paired or unpaired and acts accordingly.
"""

import sys
import os
import argparse
import pysam
from itertools import izip

from src.helperFunctions import mkdir_p, fastq_read_size

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import setLoggingFromOptions, system, logger, fastqRead

def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source_dir", help="Directory containing files.",
                        default="/hive/groups/recon/projs/mus_strain_cactus/data/assembly_rel_1505/bam/"
                                "ftp-mouse.sanger.ac.uk/REL-1505-RNA-Seq/GRCm38")
    parser.add_argument("--reference", help="kallisto reference", 
                        default="gencode_VM4_Reference/kallisto_reference/tgencode.VM4.idx")
    parser.add_argument("--out_dir", help="output location (base directory)", 
                        default="/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/rnaseq/expression")
    return parser


def wrapper(target, source_dir, reference, out_dir):
    for base_path, dirs, files in os.walk(source_dir):
        if files:
            bams = [os.path.join(source_dir, base_path, x) for x in files if x.endswith(".bam")]
            assert bams
            try:
                genome, institute, tissue = base_path.replace(source_dir, "").split("/")[1:]
            except ValueError:
                raise RuntimeError("Looks like the directory structure is not correct.")
            for bam in bams:
                target.addChildTargetFn(extract_fastq, args=(genome, institute, tissue, bam, reference, out_dir))


def extract_fastq(target, genome, institute, tissue, bam, reference, out_dir):
    out_path = build_out_dirs(out_dir, genome, institute, tissue, bam)
    if is_paired_sequencing(bam) is True:
        logger.debug("Extracting paired sequencing for bam {} with out_path {}".format(bam, out_path))
        name_sorted_sam_path = os.path.join(out_path, "name_sorted.sam")
        system("samtools sort -n -O sam -T {} {} > {}".format(target.getLocalTempDir(), bam, name_sorted_sam_path))
        fwd_fastq_path = os.path.join(out_path, os.path.basename(bam).split(".")[0] + ".fwd.fastq")
        rev_fastq_path = os.path.join(out_path, os.path.basename(bam).split(".")[0] + ".rev.fastq")
        target.setFollowOnTargetFn(get_paired_fastqs, args=(genome, institute, tissue, bam, reference, out_dir, 
                                                            name_sorted_sam_path, fwd_fastq_path, rev_fastq_path))
    else:
        logger.debug("Extracting UNPAIRED sequencing for bam {} with out_path {}".format(bam, out_path))
        fastq_path = os.path.join(out_path, os.path.basename(bam).split(".")[0] + ".fastq")
        system(r"""samtools view {} | grep -v ^@ | awk '{{print "@"$1"\n"$10"\n+\n"$11}}' > {}""".format(bam, 
                                                                                                        fastq_path))
        read_size = fastq_read_size(fastq_path)
        target.setFollowOnTargetFn(kallisto_single, args=(genome, institute, tissue, bam, reference, out_dir, 
                                                          fastq_path, read_size))


def get_paired_fastqs(target, genome, institute, tissue, bam, reference, out_dir, name_sorted_sam_path, fwd_fastq_path,
                      rev_fastq_path):
    logger.info("Extracting paired fastqs")
    target.addChildTargetFn(get_fwd, args=(name_sorted_sam_path, fwd_fastq_path))
    target.addChildTargetFn(get_rev, args=(name_sorted_sam_path, rev_fastq_path))
    target.setFollowOnTargetFn(kallisto_paired, args=(genome, institute, tissue, bam, reference, out_dir, 
                                                      name_sorted_sam_path, fwd_fastq_path, rev_fastq_path))


def get_fwd(target, name_sorted_sam_path, fwd_fastq_path):
    logger.info("Extracting fwd from {} to {}".format(name_sorted_sam_path, fwd_fastq_path))
    system(r"""cat {} | grep -v ^@ | awk 'NR%2==0 {{print "@"$1"\n"$10"\n+\n"$11}}' > {}""".format(name_sorted_sam_path,
                                                                                                   fwd_fastq_path))


def get_rev(target, name_sorted_sam_path, rev_fastq_path):
    logger.info("Extracting fwd from {} to {}".format(name_sorted_sam_path, rev_fastq_path))
    system(r"""cat {} | grep -v ^@ | awk 'NR%2==1 {{print "@"$1"\n"$10"\n+\n"$11}}' > {}""".format(name_sorted_sam_path,
                                                                                                rev_fastq_path))


def kallisto_paired(target, genome, institute, tissue, bam, reference, out_dir, name_sorted_sam_path, fwd_fastq_path,
                      rev_fastq_path):
    out_path = build_out_dirs(out_dir, genome, institute, tissue, bam)
    system("kallisto quant -i {} -o {} {} {}".format(reference, out_path, fwd_fastq_path, rev_fastq_path))


def kallisto_single(target, genome, institute, tissue, bam, reference, out_dir, fastq_path, read_size):
    out_path = build_out_dirs(out_dir, genome, institute, tissue, bam)
    system("kallisto quant --single -l {} -i {} -o {} {}".format(read_size, reference, out_path, fastq_path))


def build_out_dirs(out_dir, genome, institute, tissue, bam):
    out_path = os.path.join(out_dir, genome, institute, tissue, os.path.basename(bam).split(".")[0])
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
    from kallisto_remap import *
    main()
