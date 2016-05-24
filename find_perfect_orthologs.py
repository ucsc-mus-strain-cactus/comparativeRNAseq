"""
Find protein coding genes with perfect reciprocal 1-1 ortholog mappings of CDS sequence in a target genome.
"""
import sys
sys.path.extend(['/hive/users/ifiddes/ihategit/pipeline/', '/hive/users/ifiddes/ihategit/pipeline/submodules', '/hive/users/ifiddes/ihategit/pipeline/submodules/pycbio', '/hive/users/ifiddes/ihategit/pipeline/submodules/comparativeAnnotator'])
import argparse
import random
import tempfile
import logging
import pandas as pd
from pyfasta import Fasta
from collections import Counter, defaultdict
import pybedtools
from pybedtools.bedtool import BedTool
from pycbio.bio.bio import *
from pycbio.bio.transcripts import *
from pycbio.sys.fileOps import *
from pycbio.bio.intervals import *
from pycbio.sys.procOps import runProc, callProcLines



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('target_genome')
    parser.add_argument('--genome_dir', default='/hive/users/ifiddes/Comparative-Annotation-Toolkit/CAT/cat_work/genome_files')
    parser.add_argument('--ref-gp', default='/hive/users/ifiddes/Comparative-Annotation-Toolkit/CAT/cat_work/reference/Mus_musculus.GRCm38.83.gff3.gp')
    parser.add_argument('--attrs', default='/hive/users/ifiddes/ihategit/pipeline/gencode_vm8/C57B6J.tsv')
    parser.add_argument('--chain_dir', default='/hive/users/ifiddes/ihategit/pipeline/mouse_work_v4/chaining')
    parser.add_argument('--ref_genome', default='C57B6J')
    return parser.parse_args()


def strip_version(n):
    return n.split('.')[0]


def get_coding_genes(attrs):
    return {strip_version(x.split()[0]) for x in open(attrs) if x.split()[2] == 'protein_coding'}


def get_coding_transcripts(attrs):
    return {strip_version(x.split()[-2]) for x in open(attrs) if x.split()[-1] == 'protein_coding'}


def find_franken_splice(intervals):
    """
    Creates a franken-spliced transcript out of all exons for a gene. Does this by first merging all intervals
    as possible, then splicing them.
    """
    return gap_merge_intervals(intervals, 0)


def intervals_to_bed12(intervals, name, strand):
    """
    Converts a list of intervals to a single bed12 entry.
    """
    intervals = sorted(intervals, key=lambda x: x.start)
    start = int(intervals[0].start)
    stop = int(intervals[-1].stop)
    block_count = len(intervals)
    block_sizes = ','.join(map(str, [len(x) for x in intervals]))
    block_starts = ','.join(map(str, [x.start - start for x in intervals]))
    chrom = x.chromosome
    return '\t'.join(map(str, [chrom, start, stop, name, 0, strand, start, stop, 0, block_count, block_sizes, block_starts])) + '\n'


def load_tx_intervals(tx_dict, coding_genes, coding_transcripts):
    """
    Load exonic intervals instead of genomic
    """
    interval_map = defaultdict(list)
    for tx in tx_dict.itervalues():
        if tx.name2 not in coding_genes or tx.name not in coding_transcripts:
            continue
        cds_tx = Transcript(tx.get_bed(start_offset=tx.thick_start, stop_offset=tx.thick_stop))
        interval_map[strip_version(tx.name2)].extend(cds_tx.exon_intervals)
    r = {}
    for gene_name, intervals in interval_map.iteritems():
        chroms = Counter(x.chromosome for x in intervals)
        most_common_chrom = sorted(chroms.items())[-1][0]
        strands = Counter(x.strand for x in intervals)
        most_common_strand = sorted(strands.items())[-1][0]
        intervals = [x for x in intervals if x.chromosome == most_common_chrom and x.strand == most_common_strand]
        intervals = find_franken_splice(intervals)
        r[gene_name] = intervals_to_bed12(intervals, gene_name, convert_strand(most_common_strand))
    return r


def map_to_tgt(ref_fa, ref_sizes, tgt_fa, coding_transcripts, coding_genes, chain, tx_dict):
    bed = tmpFileGet()
    gp = tmpFileGet()
    fake_psl = tmpFileGet()
    fwd_unfiltered = tmpFileGet()
    fwd_filtered = tmpFileGet()
    r_intervals = load_tx_intervals(tx_dict, coding_genes, coding_transcripts)
    with open(bed, 'w') as outf:
        for i in r_intervals.itervalues():
            outf.write(i)
    cmd = ['bedToGenePred', bed, gp]
    runProc(cmd)
    cmd = ['genePredToFakePsl', '-chromSize={}'.format(ref_sizes), 'na', gp, fake_psl, '/dev/null']
    runProc(cmd)
    cmd = ['pslMap', '-chainMapFile', fake_psl, chain, fwd_unfiltered]
    runProc(cmd)
    cmd = [['sort', '-k10,10', fwd_unfiltered],
           ['pslCDnaFilter', '-localNearBest=0.05', '-filterWeirdOverlapped', '-decayMinCover', '/dev/stdin', fwd_filtered]]
    runProc(cmd)
    os.remove(bed)
    os.remove(gp)
    os.remove(fake_psl)
    os.remove(fwd_unfiltered)
    return fwd_filtered


def map_to_ref(fwd_filtered, chain):
    rev_unfiltered = tmpFileGet()
    cmd = ['pslMap', '-swapMap', '-chainMapFile', fwd_filtered, chain, rev_unfiltered]
    runProc(cmd)
    cmd = ['pslCDnaFilter', '-localNearBest=0.05', '-filterWeirdOverlapped', '-decayMinCover', rev_unfiltered, '/dev/stdout']
    return callProcLines(cmd)


def print_unique(fwd_filtered, rev_filtered):
    fwd_counts = Counter(x.split()[9] for x in open(fwd_filtered))
    rev_counts = Counter(x.split()[9] for x in rev_filtered)
    fwd_unique = set(x for x, y in fwd_counts.iteritems() if y == 1)
    rev_unique = set(x for x, y in rev_counts.iteritems() if y == 1)
    to_keep = fwd_unique & rev_unique
    assert len(to_keep) > 0
    for n in fwd_unique & rev_unique:
        print n
    os.remove(fwd_filtered)


def main():
    args = parse_args()
    chain = os.path.join(args.chain_dir, '{}_{}.chain.gz'.format(args.ref_genome, args.target_genome))
    coding_genes = get_coding_genes(args.attrs)
    coding_transcripts = get_coding_transcripts(args.attrs)
    tx_dict = get_transcript_dict(args.ref_gp)
    ref_fa = os.path.join(args.genome_dir, args.ref_genome + '.fa')
    ref_sizes = os.path.join(args.genome_dir, args.ref_genome + '.chrom.sizes')
    tgt_fa = os.path.join(args.genome_dir, args.target_genome + '.fa')
    fwd_filtered = map_to_tgt(ref_fa, ref_sizes, tgt_fa, coding_transcripts, coding_genes, chain, tx_dict)
    rev_filtered = map_to_ref(fwd_filtered, chain)
    print_unique(fwd_filtered, rev_filtered)


if __name__ == '__main__':
    main()