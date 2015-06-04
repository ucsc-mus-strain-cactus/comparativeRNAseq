"""
Converts all of the default STAR wiggle output files to bigWigs
"""

import sys
import os
import argparse
import re
pattern = re.compile("")

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import setLoggingFromOptions, system, logger, fastqRead

target_folder = "/cluster/home/ifiddes/mus_strain_data/pipeline_data/rnaseq/STAR_output"
chrom_sizes = "/cluster/home/ifiddes/mus_strain_data/pipeline_data/assemblies/1504/C57B6J.chrom.sizes"

def wrapper(target, target_folder, chrom_sizes):
    for base_path, dirs, files in os.walk(target_folder):
        if files:
            wig_files = [os.path.join(base_path, x) for x in files if x.endswith(".wig")]
            try:
                genome, institute, tissue, experiment = base_path.replace(target_folder, "").split("/")[1:]
            except ValueError:
                raise RuntimeError("Looks like the directory structure is not what was expected.")
            for wig in wig_files:
                target.addChildTargetFn(rename_wig, args=(wig, chrom_sizes))  


def rename_wig(target, wig, chrom_sizes):
    new_wig = wig.replace(".wig", "") + ".renamed.wig"
    with open(new_wig, "w") as outf:
        for x in open(wig):
            if 'chrom' in x:
                step, chrom = x.split()
                chrom = chrom.split("=")[1]
                if chrom == "chrM":
                    chrom = "MT"
                elif "_" in chrom:
                    chrom = chrom.split("_")[1] + ".1"
                elif 'chr' in chrom:
                    chrom = chrom.replace('chr', '')
                else:
                    chrom = chrom.split(".")[0].split("_")[1]
                x = " ".join([step, 'chrom=' + chrom]) + "\n"
            outf.write(x)
    target.setFollowOnTargetFn(convert_wig, args=(new_wig, chrom_sizes,))


def convert_wig(target, wig, chrom_sizes):
    out_bw = wig.replace(".wig", ".bw")
    system("wigToBigWig {} {} {}".format(wig, chrom_sizes, out_bw))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    setLoggingFromOptions(args)

    i = Stack(Target.makeTargetFn(wrapper, args=(target_folder, chrom_sizes))).startJobTree(args)

    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == "__main__":
    from bin.convert_star_output import *
    main()
