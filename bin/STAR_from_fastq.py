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
from src.helperFunctions import mkdir_p, fastq_read_size, find_paired_fastqs
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import setLoggingFromOptions, system, logger, getRandomAlphaNumericString


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source_dir", help="Directory containing fastq files.",
                        default="/hive/groups/recon/projs/mus_strain_cactus/data/assembly_rel_1505/fastqs/"
                                "ftp-mouse.sanger.ac.uk/REL-1505-RNA-Seq/fastq")
    parser.add_argument("--out_dir", help="output location (base directory)", required=True)
    parser.add_argument('--num_threads', default=8)
    parser.add_argument('--rsem', action='store_true', help='Run RSEM as well?. Assumes RSEM index is built in same folder as STAR index.')
    parser.add_argument('--include_tissues', nargs='+', default=None,
                        help='Space separated list of tissues to include. If set, will only do these.')
    refs = parser.add_mutually_exclusive_group(required=True)
    refs.add_argument("--reference", help="reference directory for STAR")
    refs.add_argument('--ref_dir', help='reference directory of directories. Mutually exclusive with --reference. '
                                        'Each directory should be in the format genome/ref_files')
    return parser


# star_flags variable should not need changing per run
star_flags = (" --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate  --outBAMcompression -1 "
             "--outFilterType BySJout --outSAMattributes NH HI AS NM MD --outFilterMultimapNmax 20 "
             "--outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 "
             "--alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 "
             "--sjdbScore 1 --readFilesCommand zcat --quantMode TranscriptomeSAM "
             "--outWigType wiggle --limitBAMsortRAM 80000000000")

# format star_base_cmd for each STAR call
star_paired_cmd = "STAR --genomeDir {} --readFilesIn {} {} --outFileNamePrefix {} --outTmpDir {} --runThreadN {}"
star_single_cmd = "STAR --genomeDir {} --readFilesIn {} --outFileNamePrefix {} --outTmpDir {} --runThreadN {}"

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
                    target.addChildTargetFn(run_paired_star, args=(genome, institute, tissue, reference, args.out_dir,
                                                                   experiment, fwd_fastq_path, rev_fastq_path,
                                                                   args.num_threads, args.rsem))
                except ValueError:
                    target.addChildTargetFn(run_single_star, args=(genome, institute, tissue, reference, args.out_dir,
                                                                   experiment, fastq_path,
                                                                   args.num_threads, args.rsem))


def run_paired_star(target, genome, institute, tissue, reference, out_dir, experiment, fwd_fastq_path, rev_fastq_path,
                    num_threads, rsem):
    # STAR wants a temp dir that doesn't exist, so we have to give it a fresh path because jobTree makes localTempDir()
    tmp_dir = os.path.join(target.getLocalTempDir(), "tmp_" + getRandomAlphaNumericString())
    out_path = build_out_dirs(out_dir, genome, institute, tissue, experiment) + "/"
    this_star_base_cmd = star_paired_cmd.format(reference, fwd_fastq_path, rev_fastq_path, out_path, tmp_dir, num_threads)
    system(this_star_base_cmd + star_flags)
    if rsem is True:
        target.setFollowOnTargetFn(run_rsem, args=(genome, out_path, reference, True), cpu=1)


def run_single_star(target, genome, institute, tissue, reference, out_dir, experiment, fastq_path, num_threads, rsem):
    tmp_dir = os.path.join(target.getLocalTempDir(), "tmp_" + getRandomAlphaNumericString())
    out_path = build_out_dirs(out_dir, genome, institute, tissue, experiment) + "/"
    this_star_base_cmd = star_single_cmd.format(reference, fastq_path, out_path, tmp_dir, num_threads)
    system(this_star_base_cmd + star_flags)
    if rsem is True:
        target.setFollowOnTargetFn(run_rsem, args=(genome, out_path, reference, False), cpu=1)


def run_rsem(target, genome, out_path, reference, is_paired):
    ref = os.path.join(reference, genome)
    bam = os.path.join(out_path, 'Aligned.toTranscriptome.out.bam')
    out_rsem = os.path.join(out_path, 'RSEM')
    cmd = ['rsem-calculate-expression', '--bam', bam, '--temporary-folder', target.getLocalTempDir(), ref, out_rsem]
    if is_paired is True:
        cmd.append('--paired-end')
    system(' '.join(map(str, cmd)))


def build_out_dirs(out_dir, genome, institute, tissue, experiment):
    out_path = os.path.join(out_dir, genome, institute, tissue, experiment)
    mkdir_p(out_path)
    return out_path


def main():
    parser = build_parser()
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    setLoggingFromOptions(args)
    args.defaultCpu = args.num_threads
    args.defaultMemory = 38 * 1024 ** 3
    i = Stack(Target.makeTargetFn(wrapper, args=(args,), memory=args.defaultMemory, cpu=args.defaultCpu)).startJobTree(args)

    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == "__main__":
    from bin.STAR_from_fastq import *
    main()