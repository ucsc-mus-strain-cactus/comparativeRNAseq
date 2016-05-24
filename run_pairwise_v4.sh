#!/bin/bash
set -e

# Aligning the following sets of data:

# 1. sanger/brain (all strains)
# 2. jax/liver (8 strains)
# 3. ebi/heart (wild strains)
# 4. ebi/kidney (wild strains)
# 5. ebi/liver (wild strains)

export PYTHONPATH=./:/hive/users/ifiddes/ihategit/pipeline/submodules/:/hive/users/ifiddes/ihategit/pipeline/submodules/pycbio/

base_dir=star_alignments_v4
ref_dir=${base_dir}/mm10
tgt_dir=${base_dir}/strain_specific
jobtree_base_dir=star_alignments_v4_jobTrees

strain_assembly_refs=STAR_references_v4/
mm10_ref=mm10_reference/C57B6J


run_tgt_align() {
    tissue=$1
    institute=$2
    combo=${tissue}_${institute}
    jobtree_dir=${jobtree_base_dir}/tgt_${combo}
    mkdir -p `dirname ${jobtree_dir}`
    out_dir=${tgt_dir}/${combo}
    mkdir -p ${out_dir}
    python bin/STAR_from_bam.py  --ref_dir ${strain_assembly_refs} --out_dir ${out_dir} \
--jobTree ${jobtree_dir} --batchSystem parasol --parasolCommand ~/bin/remparasol --include_tissues ${tissue} \
--include_institutes ${institute}
}

run_ref_align() {
    tissue=$1
    institute=$2
    combo=${tissue}_${institute}
    jobtree_dir=${jobtree_base_dir}/ref_${combo}
    mkdir -p `dirname ${jobtree_dir}`
    out_dir=${ref_dir}/${combo}
    mkdir -p ${out_dir}
    python bin/STAR_from_bam.py  --reference ${mm10_ref} --out_dir ${out_dir} \
--jobTree ${jobtree_dir} --batchSystem parasol --parasolCommand ~/bin/remparasol --include_tissues ${tissue} \
--include_institutes ${institute} --maxThreads 40
}


tissue=brain
institute=sanger
run_tgt_align ${tissue} ${institute}
run_ref_align ${tissue} ${institute}

tissue=liver
institute=jax
run_tgt_align ${tissue} ${institute}
run_ref_align ${tissue} ${institute}

tissue=heart
institute=ebi
run_tgt_align ${tissue} ${institute}
run_ref_align ${tissue} ${institute}

tissue=kidney
institute=ebi
run_tgt_align ${tissue} ${institute}
run_ref_align ${tissue} ${institute}

tissue=liver
institute=ebi
run_tgt_align ${tissue} ${institute}
run_ref_align ${tissue} ${institute}
