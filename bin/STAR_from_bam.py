#!/usr/bin/env python
"""
Realigns BAM files to new transcriptomes using STAR.

Hunts through --dir (which is a FTP download from ftp://ftp-mouse.sanger.ac.uk/) and runs STAR aligning to the
transcriptome. Designed with 1505 fastq release in mind. The starting directory structure will be copied.

Attempts to determine if the input RNAseq was paired or unpaired and acts accordingly.
"""

import sys
import os
import argparse
import pysam
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import setLoggingFromOptions, getRandomAlphaNumericString
from pycbio.sys.procOps import runProc
from pycbio.sys.fileOps import tmpFileGet, ensureDir


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source_dir", help="Directory containing BAM files.",
                        default="/hive/groups/recon/projs/mus_strain_cactus/data/rel-1509-rna-seq/ftp-mouse.sanger.ac.uk/REL-1509-Assembly-RNA-Seq")
    parser.add_argument("--out_dir", help="output location (base directory)", required=True)
    parser.add_argument('--num_threads', default=8)
    parser.add_argument('--htseq', action='store_true', help='Run HTSeq-count as well? Assumes a GTF with the genome name is in the same folder as the STAR index.')
    parser.add_argument('--include_tissues', nargs='+', default=None,
                        help='Space separated list of tissues to include. If set, will only do these.')
    parser.add_argument('--include_genomes', nargs='+', default=None,
                        help='Space separated list of genomes to include.')
    parser.add_argument('--include_institutes', nargs='+', default=None,
                        help='Space separated list of institutes to include.')
    refs = parser.add_mutually_exclusive_group(required=True)
    refs.add_argument("--reference", help="reference directory for STAR")
    refs.add_argument('--ref_dir', help='reference directory of directories. Mutually exclusive with --reference. '
                                        'Each directory should be in the format genome/ref_files')
    return parser


# star_flags variable should not need changing per run
star_flags = ['--outSAMunmapped', 'Within', '--outSAMtype', 'BAM', 'Unsorted', 'SortedByCoordinate', '--outBAMcompression', '-1',
              '--outFilterType', 'BySJout', '--outSAMattributes', 'NH', 'HI', 'AS', 'NM', 'MD',
              '--outFilterMultimapNmax', '20', '--outFilterMismatchNmax', '999', '--outFilterMismatchNoverReadLmax',
              '0.04', '--alignIntronMin', '20', '--alignIntronMax', '1000000', '--alignMatesGapMax', '1000000',
              '--alignSJoverhangMin', '8', '--alignSJDBoverhangMin', '1', '--sjdbScore', '1', '--limitBAMsortRAM', '80000000000']


def is_paired_sequencing(bamfile):
    # TODO: this is scary. Should check for unpaired being 0, and number paired == total number
    r = pysam.flagstat(bamfile)
    paired = int(r[5].split()[0])
    if paired != 0:
        return True
    else:
        return False


def fastq_read_size(fastq_path, num_reads=10000):
    sizes = []
    fastq_handle = fastqRead(gzip.open(fastq_path))
    for i in xrange(num_reads):
        name, seq, qual = fastq_handle.next()
        sizes.append(len(seq))
    return 1.0 * sum(sizes) / len(sizes)


def wrapper(target, args):
    for base_path, dirs, files in os.walk(args.source_dir):
        if files:
            try:
                genome, institute, tissue = base_path.replace(args.source_dir, "").split("/")[1:]
            except ValueError:
                raise RuntimeError("Looks like the directory structure is not what we expected.")
            if args.include_tissues is not None and tissue not in args.include_tissues:
                continue
            elif args.include_genomes is not None and genome not in args.include_genomes:
                continue
            elif args.include_institutes is not None and institute not in args.include_institutes:
                continue
            bam_files = [x for x in files if x.endswith('.bam')]
            bam_map = {x.split('Aligned')[0]: os.path.join(base_path, x) for x in bam_files}
            if args.reference is not None:
                reference = args.reference
                ref_genome = os.path.basename(reference)  # hack
                assert ref_genome == 'C57B6J'
                for experiment, bam_path in bam_map.iteritems():
                    target.addChildTargetFn(extract_fastq, args=(genome, institute, tissue, reference, args.out_dir,
                                                                 experiment, bam_path, args.num_threads, ref_genome))
            else:
                reference = os.path.join(args.ref_dir, genome)
                for experiment, bam_path in bam_map.iteritems():
                    target.addChildTargetFn(extract_fastq, args=(genome, institute, tissue, reference, args.out_dir,
                                                                 experiment, bam_path, args.num_threads, None))


def extract_fastq(target, genome, institute, tissue, reference, out_dir, experiment, bam_path, num_threads, ref_genome):
    if is_paired_sequencing(bam_path):
        fwd_fastq_path = tmpFileGet(prefix=experiment, suffix='fwd.fastq')
        rev_fastq_path = tmpFileGet(prefix=experiment, suffix='rev.fastq')
        cmd = ['samtools', 'fastq', '-1', fwd_fastq_path, '-2', rev_fastq_path, bam_path]
        runProc(cmd)
        run_paired_star(target, genome, institute, tissue, reference, out_dir, experiment, fwd_fastq_path,
                        rev_fastq_path, num_threads, ref_genome)
    else:
        fastq_path = tmpFileGet(prefix=experiment)
        cmd = ['samtools', 'fastq', '-0', fastq_path, bam_path]
        runProc(cmd)
        run_single_star(target, genome, institute, tissue, reference, out_dir, experiment, fastq_path,
                        num_threads, ref_genome)


def run_paired_star(target, genome, institute, tissue, reference, out_dir, experiment, fwd_fastq_path, rev_fastq_path,
                    num_threads, ref_genome):
    # STAR wants a temp dir that doesn't exist, so we have to give it a fresh path because jobTree makes localTempDir()
    tmp_dir = os.path.join(target.getLocalTempDir(), "tmp_" + getRandomAlphaNumericString())
    out_path = build_out_dirs(out_dir, genome, institute, tissue, experiment) + "/"
    star_cmd = ['STAR', '--genomeDir', reference, '--readFilesIn', fwd_fastq_path, rev_fastq_path,
                '--outFileNamePrefix', out_path, '--outTmpDir', tmp_dir, '--runThreadN', num_threads]
    runProc(star_cmd + star_flags)
    os.remove(fwd_fastq_path)
    os.remove(rev_fastq_path)
    target.setFollowOnTargetFn(run_feature_counts, args=(genome, out_path, reference, ref_genome, num_threads, True))


def run_single_star(target, genome, institute, tissue, reference, out_dir, experiment, fastq_path, num_threads, ref_genome):
    tmp_dir = os.path.join(target.getLocalTempDir(), "tmp_" + getRandomAlphaNumericString())
    out_path = build_out_dirs(out_dir, genome, institute, tissue, experiment) + "/"
    star_cmd = ['STAR', '--genomeDir', reference, '--readFilesIn', fastq_path,
                '--outFileNamePrefix', out_path, '--outTmpDir', tmp_dir, '--runThreadN', num_threads]
    runProc(star_cmd + star_flags)
    os.remove(fastq_path)
    target.setFollowOnTargetFn(run_feature_counts, args=(genome, out_path, reference, ref_genome, num_threads, False))


def run_feature_counts(target, genome, out_path, reference, ref_genome, num_threads, is_paired=False):
    if ref_genome is None:
        ref = os.path.join(reference, genome + '.gtf')
    else:
        ref = os.path.join(reference, ref_genome + '.gtf')
    assert os.path.exists(ref), ref
    bam = os.path.join(out_path, 'Aligned.out.bam')
    assert os.path.exists(bam), bam
    out_counts = os.path.join(out_path, 'Aligned.out.bam.counts.cds')
    cmd = ['featureCounts', '-T', num_threads, '-t', 'CDS', '-g', 'gene_id', '-a', ref, '-o', out_counts, bam]
    if is_paired is True:
        cmd.append('-p')
    runProc(cmd)


def build_out_dirs(out_dir, genome, institute, tissue, experiment):
    out_path = os.path.join(out_dir, genome, institute, tissue, experiment)
    ensureDir(out_path)
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
    from bin.STAR_from_bam import *
    main()