"""
Munges all STAR output files into a target directory usable for making browser tracks.
"""

import sys, os, subprocess
from src.helperFunctions import mkdir_p, name_map


junction_map = {"1": "GT_AG", "2": "CT_AC", "3": "GC_AG", "4": "CT_GC", "5": "AT_AC", "6": "GT_AT"}

def rename_to_ensembl(chrom):
    if chrom == "chrM":
        chrom = "MT"
    elif "_" in chrom:
        chrom = chrom.split("_")[1] + ".1"
    elif 'chr' in chrom:
        chrom = chrom.replace('chr', '')
    else:
        chrom = chrom.split(".")[0].split("_")[1]
    return chrom


def splice_junction_line_to_bed(splice_junction_record, experiment, rename=None):
    """
    Specific case of turning an intron interval into the first and last two bases from a splice junction out file
    from STAR. This file uses 1 based positions.
    """
    chrom, start, stop, strand, motif, annot, uniq_map, multi_map, overhang = splice_junction_record.split()
    start = int(start) - 1
    stop = int(stop) - 1
    block_starts = "0,{}".format(stop - start - 2)
    if strand == "0":
        # undefined strand
        strand = "."
        color = "74,74,74"
    elif strand == "1":
        # plus strand
        strand = "+"
        color = "144,212,149"
    else:
        # minus strand
        assert strand == "2"
        strand = "-"
        color = "212,144,207"
    if rename is not None:
        chrom = rename(chrom)
    # hack to deal with weird entries where stop is within 2bp of start: bug in STAR?
    if stop - start - 2 < 2:
        return tuple(map(str, [chrom, start, stop, experiment, 0, strand, start, stop, color, "1", stop - start, "0"]))
    # TODO: map interval to actual gene name
    return tuple(map(str, [chrom, start, stop, experiment, 0, strand, start, stop, color, "2", "2,2", block_starts]))


def munge_files(files, target_path, rename=None):
    mkdir_p(os.path.dirname(target_path))
    records = set()
    for f in files:
        experiment = os.path.basename(f).replace("SJ.out.tab","")
        if len(experiment) == 0:
            experiment = os.path.dirname(f).split("/")[-1]
        for line in open(f):
            records.add(splice_junction_line_to_bed(line, experiment, rename))
    records = sorted(records, key = lambda x: [x[0], x[1]])
    with open(target_path, "w") as outf:
        for line in records:
            outf.write("\t".join(line) + "\n")

# the below code combines splice junction TAB files across experiments
# the STAR alignments were done here at UCSC and are against the reference genome

file_map = {}
for base_path, dirs, files in os.walk("/cluster/home/ifiddes/mus_strain_data/pipeline_data/rnaseq/STAR_output"):
    if files:
        genome, institute, tissue, experiment = base_path.split("/")[-4:]
        sj_path = os.path.join(base_path, "SJ.out.tab")
        if (genome, institute, tissue) not in file_map:
            file_map[(genome, institute, tissue)] = []
        file_map[(genome, institute, tissue)].append(sj_path)

for (genome, institute, tissue), files in file_map.iteritems():
    target_path = os.path.join("/cluster/home/ifiddes/mus_strain_data/pipeline_data/rnaseq/munged_STAR_data/GRCm38", genome, institute + "_" + tissue + ".sj.bed")
    munge_files(files, target_path, rename=rename_to_ensembl)

for (genome, institute, tissue), files in file_map.iteritems():
    target_path = os.path.join("/cluster/home/ifiddes/mus_strain_data/pipeline_data/rnaseq/munged_STAR_data/GRCm38", genome, institute + "_" + tissue + ".sj.bed")
    subprocess.Popen("bedSort {0} {0}".format(target_path), shell=True)

chrom_sizes = "/cluster/home/ifiddes/mus_strain_data/pipeline_data/assemblies/1504/C57B6J.chrom.sizes"
for (genome, institute, tissue), files in file_map.iteritems():
    bed_path = os.path.join("/cluster/home/ifiddes/mus_strain_data/pipeline_data/rnaseq/munged_STAR_data/GRCm38", genome, institute + "_" + tissue + ".sj.bed")
    bb_path = os.path.join("/cluster/home/ifiddes/mus_strain_data/pipeline_data/rnaseq/munged_STAR_data/GRCm38", genome, institute + "_" + tissue + ".sj.bb")
    subprocess.Popen("bedToBigBed -type=bed12 {} {} {}".format(bed_path, chrom_sizes, bb_path), shell=True)

# the below code combines splice junction TAB files across experiments
# the STAR alignments were done by Sanger and are against each indivdiual genome

file_map = {}
for base_path, dirs, files in os.walk("/cluster/home/ifiddes/mus_strain_data/data/assembly_rel_1505/splice_junction_data/REL-1504-chromosomes"):
    if files:
        genome, institute, tissue = base_path.split("/")[-3:]
        genome = name_map[genome]
        file_paths = [os.path.join(base_path, x) for x in files]
        if (genome, institute, tissue) not in file_map:
            file_map[(genome, institute, tissue)] = []
        file_map[(genome, institute, tissue)].extend(file_paths)


for (genome, institute, tissue), files in file_map.iteritems():
    target_path = os.path.join("/cluster/home/ifiddes/mus_strain_data/pipeline_data/rnaseq/munged_STAR_data/REL-1504-chromosomes", genome, institute + "_" + tissue + ".sj.bed")
    munge_files(files, target_path, rename=None)

for (genome, institute, tissue), files in file_map.iteritems():
    target_path = os.path.join("/cluster/home/ifiddes/mus_strain_data/pipeline_data/rnaseq/munged_STAR_data/REL-1504-chromosomes", genome, institute + "_" + tissue + ".sj.bed")
    subprocess.Popen("bedSort {0} {0}".format(target_path), shell=True)

chrom_size_base = "/cluster/home/ifiddes/mus_strain_data/pipeline_data/assemblies/1504/{}.chrom.sizes"
for (genome, institute, tissue), files in file_map.iteritems():
    bed_path = os.path.join("/cluster/home/ifiddes/mus_strain_data/pipeline_data/rnaseq/munged_STAR_data/REL-1504-chromosomes", genome, institute + "_" + tissue + ".sj.bed")
    bb_path = os.path.join("/cluster/home/ifiddes/mus_strain_data/pipeline_data/rnaseq/munged_STAR_data/REL-1504-chromosomes", genome, institute + "_" + tissue + ".sj.bb")
    chrom_sizes = chrom_size_base.format(genome)
    subprocess.Popen("bedToBigBed -type=bed12 {} {} {}".format(bed_path, chrom_sizes, bb_path), shell=True)
