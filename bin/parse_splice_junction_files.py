#!/usr/bin/env python
"""
Digs through the splice junction output from STAR and parses it to create BED files.
For biological replicates, the splice junctions will be combined into one set.
"""

import sys
import os
import argparse
import errno

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target_dir", help="Splice junction raw data dir", 
                        default="/hive/groups/recon/projs/mus_strain_cactus/data/assembly_rel_1505/"
                                "splice_junction_data")
    parser.add_argument("--out_dir", help="output dir", default="/cluster/home/ifiddes/mus_strain_data/pipeline_data/"
                                                                "rnaseq/splice_junction_bedfiles")
    return parser.parse_args()


def find_splice_junction_files(target_dir):
    file_map = {}
    for base_path, dirs, files in os.walk(target_dir):
        if files:
            assert not dirs
            assert all([x.endswith("tab") for x in files])
            target_genome, genome, institute, tissue = base_path.replace(target_dir, "").split("/")[1:]
            if target_genome not in file_map:
                file_map[target_genome] = {}
            if genome not in file_map[target_genome]:
                file_map[target_genome][genome] = {}
            if tissue not in file_map[target_genome][genome]:
                file_map[target_genome][genome][tissue] = []
            file_map[target_genome][genome][tissue] = [os.path.join(base_path, x) for x in files]
    return file_map


junction_map = {"1": "GT_AG", "2": "CT_AC", "3": "GC_AG", "4": "CT_GC", "5": "AT_AC", "6": "GT_AT"}

def splice_junction_line_to_bed(splice_junction_record):
    """
    Specific case of turning an intron interval into the first and last two bases from a splice junction out file
    from STAR. This file uses 1 based positions.
    """
    chrom, start, stop, strand, motif, annot, uniq_map, multi_map, overhang = splice_junction_record.split()
    start = int(start) - 1
    stop = int(stop) - 1
    block_starts = "0,{}".format(stop - start - 2)
    name = []
    if strand == "0":
        strand = "+"
        name.append("unknown_strand")
    elif strand == "1":
        strand = "+"
    else:
        strand = "-"
    if motif == "0":
        name.append("non_canonical")
    else:
        name.append(junction_map[motif])
    name = "/".join(name)
    # TODO: map interval to actual gene name
    return "\t".join(map(str, [chrom, start, stop, name, 0, strand, start, stop, 0, "2,2", block_starts])) + "\n"


def find_splice_junction_set(splice_junction_files):
    records = set()
    for infile in splice_junction_files:
        for line in open(infile):
            records.add(splice_junction_line_to_bed(line))
    sorted_records = sorted(records, key = lambda x: (x[0], x[1]))   # sort by chromosome then by start position
    return sorted_records


def get_record_map(file_map):
    record_map = {}
    for target_genome in file_map:
        if target_genome not in record_map:
            record_map[target_genome] = {}
        for genome in file_map[target_genome]:
            if genome not in record_map[target_genome]:
                record_map[target_genome][genome] = {}
            for tissue in file_map[target_genome][genome]:
                records = find_splice_junction_set(file_map[target_genome][genome][tissue])
                record_map[target_genome][genome][tissue] = records
    return record_map    


def main():
    args = parse_args()
    file_map = find_splice_junction_files(args.target_dir)
    record_map = get_record_map(file_map)
    for target_genome in record_map:
        for genome in record_map[target_genome]:
            for tissue in record_map[target_genome][genome]:
                out_path = os.path.join(args.out_dir, target_genome, genome, tissue)
                mkdir_p(out_path)
                with open(os.path.join(out_path, "splice_junctions.bed"), "w") as outf:
                    for x in record_map[target_genome][genome][tissue]:
                        outf.write(x)


if __name__ == "__main__":
    main()