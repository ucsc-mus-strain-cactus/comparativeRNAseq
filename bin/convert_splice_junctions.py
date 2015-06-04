"""
Changes all of the default STAR SJ.out.tab files into BED files.
"""

import sys
import os

target_folder = "/cluster/home/ifiddes/mus_strain_data/pipeline_data/rnaseq/STAR_output"

def splice_junction_line_to_bed(splice_junction_record):
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
    if chrom == "chrM":
        chrom = "MT"
    elif "_" in chrom:
        chrom = chrom.split("_")[1] + ".1"
    elif 'chr' in chrom:
        chrom = chrom.replace('chr', '')
    else:
        chrom = chrom.split(".")[0].split("_")[1]
    # hack to deal with weird entries where stop is within 2bp of start: bug in STAR?
    if stop - start - 2 < 2:
        return "\t".join(map(str, [chrom, start, stop, "SJ", 0, strand, start, stop, color, "0", stop - start, "0"])) + "\n"
    # TODO: map interval to actual gene name
    return "\t".join(map(str, [chrom, start, stop, "SJ", 0, strand, start, stop, color, "2", "2,2", block_starts])) + "\n"


def splice_junction_wrapper(sj_path):
    in_file_handle = open(sj_path)
    out_file_path = os.path.join(os.path.dirname(sj_path), "sj.bed")
    with open(out_file_path, "w") as out_file_handle:
        for x in in_file_handle:
            out_file_handle.write(splice_junction_line_to_bed(x))
    in_file_handle.close()


def main():
    for base_path, dirs, files in os.walk(target_folder):
        if files:
            sj_path = os.path.join(base_path, "SJ.out.tab")
            assert os.path.exists(sj_path)
            splice_junction_wrapper(sj_path)


if __name__ == "__main__":
    main()