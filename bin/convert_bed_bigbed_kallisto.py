import sys, os, subprocess

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

for base_dir, dirs, files in os.walk("/cluster/home/ifiddes/mus_strain_data/pipeline_data/rnaseq/log_expression_bedfiles"):
    if files:
        for f in files:
            with open(os.path.join(base_dir, f + ".renamed"), "w") as outf:
                for line in open(os.path.join(base_dir, f)):
                    line = line.split()
                    line[0] = rename_to_ensembl(line[0])
                    outf.write("\t".join(line) + "\n")

for base_dir, dirs, files in os.walk("/cluster/home/ifiddes/mus_strain_data/pipeline_data/rnaseq/log_expression_bedfiles"):
    if files:
        for f in files:
            if f.endswith("renamed"):
                f = os.path.join(base_dir, f)
                subprocess.Popen("bedSort {0} {0}".format(f), shell=True)



chrom_sizes = "~/mus_strain_data/pipeline_data/assemblies/1504/C57B6J.chrom.sizes"
for base_dir, dirs, files in os.walk("/cluster/home/ifiddes/mus_strain_data/pipeline_data/rnaseq/log_expression_bedfiles"):
    if files:
        for f in files:
            if f.endswith("renamed"):
                new = f + ".bb"
                f = os.path.join(base_dir, f)
                new_f = os.path.join(base_dir, new)
                subprocess.Popen("bedToBigBed {} {} {}".format(f, chrom_sizes, new_f), shell=True)
