"""
Makes symlinks for all BAMs for browser track building
"""

import os
from lib.general_lib import mkdir_p
# make symlinks for all alignments done against reference by Ian

file_map = {}
for base_path, dirs, files in os.walk("/cluster/home/ifiddes/mus_strain_data/pipeline_data/rnaseq/STAR_output"):
    if files:
        genome, institute, tissue, experiment = base_path.split("/")[-4:]
        bam_path = os.path.join(base_path, "Aligned.sortedByCoord.out.bam")
        bai_path = os.path.join(base_path, "Aligned.sortedByCoord.out.bam.bai")
        if (genome, institute, tissue) not in file_map:
            file_map[(genome, institute, tissue)] = []
        file_map[(genome, institute, tissue)].append([experiment, bam_path, bai_path])

target_folder = "/cluster/home/ifiddes/mus_strain_data/pipeline_data/rnaseq/munged_STAR_data/GRCm38"
for (genome, institute, tissue), files in file_map.iteritems():
    for experiment, bam_path, bai_path in files:
        bam = institute + "_" + tissue + "_" + experiment + ".sortedByCoord.bam"
        os.symlink(bam_path, os.path.join(target_folder, genome, bam))
        os.symlink(bai_path, os.path.join(target_folder, genome, bam + ".bai"))

# make symlinks for all alignments done against target assemblies by Sanger

file_map = {}
for base_path, dirs, files in os.walk("/hive/groups/recon/projs/mus_strain_cactus/data/rel-1509-rna-seq/ftp-mouse.sanger.ac.uk/REL-1509-Assembly-RNA-Seq"):
    if files:
        genome, institute, tissue = base_path.split("/")[-3:]
        #genome = name_map[genome]
        for f in files:
            if f.endswith(".bam"):
                experiment = f.split("Aligned")[0]
                bam_path = os.path.join(base_path, f)
                bai_path = bam_path + ".bai"
                if (genome, institute, tissue) not in file_map:
                    file_map[(genome, institute, tissue)] = []
                file_map[(genome, institute, tissue)].append([experiment, bam_path, bai_path])


target_folder = "/cluster/home/ifiddes/mus_strain_data/pipeline_data/rnaseq/munged_STAR_data/REL-1509-chromosomes"
for (genome, institute, tissue), files in file_map.iteritems():
    for experiment, bam_path, bai_path in files:
        bam = institute + "_" + tissue + "_" + experiment + ".sortedByCoord.bam"
        mkdir_p(os.path.join(target_folder, genome))
        tgt_bam = os.path.join(target_folder, genome, bam)
        tgt_bai = os.path.join(target_folder, genome, bam + ".bai")
        if not os.path.exists(tgt_bam):
            os.symlink(bam_path, tgt_bam)
            print "linked {}".format(tgt_bam)
        if not os.path.exists(tgt_bai):
            os.symlink(bai_path, tgt_bai)