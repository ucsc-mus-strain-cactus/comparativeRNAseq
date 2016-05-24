"""
Can we improve our metrics by filtering for reads which align uniquely in BOTH datasets?
"""
import os
import sys
sys.path.extend(['/hive/users/ifiddes/ihategit/pipeline/', '/hive/users/ifiddes/ihategit/pipeline/submodules', '/hive/users/ifiddes/ihategit/pipeline/submodules/pycbio', '/hive/users/ifiddes/ihategit/pipeline/submodules/comparativeAnnotator'])
import pysam
from pycbio.sys.procOps import runProc
from collections import defaultdict, Counter
from glob import glob

genome = sys.argv[1]
experiment = sys.argv[2]
mm10_dir = 'mm10'
strain_dir = 'strain_specific'
strain_ref = os.path.join('STAR_references_v4', genome)
mm10_ref = 'mm10_reference/C57B6J/'
strain_gtf = os.path.join(strain_ref, genome + '.gtf')
mm10_gtf = os.path.join(mm10_ref, 'C57B6J.gtf')


def is_paired_sequencing(bamfile):
    # TODO: this is scary. Should check for unpaired being 0, and number paired == total number
    r = pysam.flagstat(bamfile)
    paired = int(r[5].split()[0])
    if paired != 0:
        return True
    else:
        return False


def find_not_unique_reads(bam1, bam2, is_paired):
    """
    Find the names of reads in both bams which have more than one mapping.
    """
    expected_value = 2 if is_paired else 1
    c1 = Counter(x.qname for x in pysam.Samfile(bam1))
    c2 = Counter(x.qname for x in pysam.Samfile(bam2))
    c1_unique = {x for x, y in c1.iteritems() if y == expected_value}
    c2_unique = {x for x, y in c2.iteritems() if y == expected_value}
    return c1_unique & c2_unique


def filter_bam(bam, out_bam, names):
    """
    Filter a BAM for reads in the read set
    """
    bam_handle = pysam.Samfile(bam)
    outf = pysam.Samfile(out_bam, 'wb', template=bam_handle)
    for x in bam_handle:
        if x.qname in names:
            outf.write(x)


def get_files():
    mm10_files = glob(os.path.join('star_alignments_v4/', mm10_dir, experiment, genome, '*', '*', '*', 'Aligned.out.bam'))
    strain_files = glob(os.path.join('star_alignments_v4/', strain_dir, experiment, genome, '*', '*', '*', 'Aligned.out.bam'))
    assert len(mm10_files) == len(strain_files)
    file_map = defaultdict(list)
    for p in [mm10_files, strain_files]:
        for f in p:
            e = os.path.basename(os.path.dirname(f))  # highest directory name
            file_map[e].append(f)
    assert all([len(x) == 2 for x in file_map.itervalues()])
    return file_map


def main():
    file_map = get_files()
    names = defaultdict(set)
    for e, (bam1, bam2) in file_map.iteritems():
        is_paired = is_paired_sequencing(bam1)
        unique_names = find_not_unique_reads(bam1, bam2, is_paired)
        for b in [bam1, bam2]:
            out_b = b + '.filtered.bam'
            filter_bam(b, out_b, unique_names)
            ref = mm10_gtf if 'mm10' in b else strain_gtf
            out_counts = out_b + '.counts.cds'
            cmd = ['featureCounts', '-T', '1', '-t', 'CDS', '-g', 'gene_id', '--primary', '--ignoreDup', '-Q', '30',
                   '-a', ref, '-o', out_counts, out_b]
            if is_paired is True:
                cmd.append('-p')
            runProc(cmd)


if __name__  == '__main__':
    main()