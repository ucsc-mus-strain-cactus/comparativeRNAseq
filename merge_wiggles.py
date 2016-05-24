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


def main():
    base_cmd = ['/cluster/home/ifiddes/WiggleTools/bin/wiggletools', 'median']
    for m in [mm10_dir, strain_dir]:
        files = glob(os.path.join('star_alignments_v4/', m, experiment, genome, '*', '*', '*', 'Aligned.sortedByCoord.out.bam.coverage.bw'))
        assert len(files) > 0
        cmd = base_cmd + files
        outf = os.path.join('star_alignments_v4', m, experiment, genome, 'median_counts.bw')
        runProc(cmd, stdout=outf)



if __name__  == '__main__':
    main()