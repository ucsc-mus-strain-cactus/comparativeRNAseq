"""
Merges expression wiggle files into one per tissue-institute.
Uses the java-genomics-toolkit for its awesome wigmath.Average function
"""
import sys, os, subprocess, multiprocessing, shutil
source_dir = "/hive/groups/recon/projs/mus_strain_cactus/data/rel-1509-rna-seq/ftp-mouse.sanger.ac.uk/REL-1509-Assembly-RNA-Seq"
target_dir = "/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/rnaseq/munged_STAR_data/REL-1509-chromosomes"
ref_chrom_sizes = "/cluster/home/ifiddes/mus_strain_data/pipeline_data/assemblies/1509_v2/{}.chrom.sizes"

wigmath_cmd = "./toolRunner.sh wigmath.Average -o {} {}"
convert_cmd = "wigToBigWig -clip {} {} {}"



name_map = {"CBA_J": "CBAJ", "BALB_cJ": "BALBcJ", "PWK_PhJ": "PWKPhJ", "NOD_ShiLtJ": "NODShiLtJ",
            "A_J": "AJ", "LP_J": "LPJ", "WSB_EiJ": "WSBEiJ", "Pahari_EiJ": "PAHARIEiJ",
            "C57BL_6NJ": "C57B6NJ", "AKR_J": "AKRJ", "DBA_2J": "DBA2J", "NZO_HlLtJ": "NZOHlLtJ",
            "129S1_SvImJ": "129S1", "CAST_EiJ": "CASTEiJ", "CAROLI_EiJ": "CAROLIEiJ",
            "C3H_HeJ": "C3HHeJ", "SPRET_EiJ": "SPRETEiJ", "FVB_NJ": "FVBNJ"}


def build_map(file_map, path, genome, institute, tissue, experiment):
    if genome not in file_map:
        file_map[genome] = {}
    if tissue not in file_map[genome]:
        file_map[genome][tissue] = {}
    if institute not in file_map[genome][tissue]:
        file_map[genome][tissue][institute] = {}
    file_map[genome][tissue][institute][experiment] = path


file_map = {}
for base_path, dirs, files in os.walk(source_dir):
    if files:
        genome, institute, tissue = base_path.split("/")[-3:]
        if genome not in ["PWK_PhJ"]:
            continue
        #genome = name_map[genome]
        for f in files:
            if f.endswith(".wig"):
                experiment = f.split("Signal")[0]
                wig_path = os.path.join(base_path, f)
                build_map(file_map, wig_path, genome, institute, tissue, experiment)

cmds = []
singletons = []
for genome in file_map:
    for tissue in file_map[genome]:
        for institute in file_map[genome][tissue]:
            files = []
            for experiment in file_map[genome][tissue][institute]:
                files.append(file_map[genome][tissue][institute][experiment])
            out_wig_path = os.path.realpath(os.path.join(target_dir, genome, institute + "_" + tissue + ".expression.wig"))
            if os.path.exists(out_wig_path):
                continue
            if len(files) < 2:
                singletons.append([files[0], out_wig_path])
            else:
                file_paths = " ".join(files)
                cmds.append(wigmath_cmd.format(out_wig_path, file_paths))


def run_cmd(c):
    subprocess.call(c, shell=True)


pool = multiprocessing.Pool(processes=5)

r = pool.map(run_cmd, cmds)
for inf, outf in singletons:
    shutil.copy(inf, outf)


for g in ["PWK_PhJ"]:
    for f in os.listdir(os.path.join(target_dir, g)):
        if not f.endswith("wig"):
            continue
        outf = os.path.join(target_dir, g, f.replace(".wig", ".bw"))
        f = os.path.join(target_dir, g, f)
        r = subprocess.Popen(convert_cmd.format(f, ref_chrom_sizes.format(g), outf), shell=True)
