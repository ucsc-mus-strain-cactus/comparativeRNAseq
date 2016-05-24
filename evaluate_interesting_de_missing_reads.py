# evaluating an interesting gene
import pysam
region = ['9', 61913284, 61914538]
names = set()
sh = pysam.Samfile("/hive/users/ifiddes/comparativeRNAseq/aligned_reference/C57BL_6NJ/sanger/brain/sanger_brain_ERR033015.sortedByCoord.bam")
for rec in sh.fetch(*region):
    names.add(rec.qname)


# find these names in the aligned to NJ file
sh_nj = pysam.Samfile("/hive/users/ifiddes/comparativeRNAseq/aligned_new_reference/C57BL_6NJ/sanger/brain/ERR033015/Aligned.sortedByCoord.out.bam")
recs = []
for rec in sh_nj:
    if rec.qname in names:
        recs.append(rec)

# how many alternate mappings did the reference have?
ref_recs = []
sh = pysam.Samfile("/hive/users/ifiddes/comparativeRNAseq/aligned_reference/C57BL_6NJ/sanger/brain/sanger_brain_ERR033015.sortedByCoord.bam")
for rec in sh:
    if rec.qname in names:
        ref_recs.append(rec)

outf = pysam.Samfile('nj_reads.bam', 'wb', template=sh_nj)
for x in recs:
    outf.write(x)


>>> len(names), len(recs), len(ref_recs)
(2050, 4112, 4616)


>>> len(names), len([x for x in recs if not x.is_secondary and not x.is_supplementary]), len([x for x in ref_recs if not x.is_secondary and not x.is_supplementary])
(2050, 4100, 4100)

# ref has 400 more secondary mappings
# none of this really makes sense... why is expression so different?

# htseq count raw data:
# ENSMUSG00000007892.7    2529    1811    77      71
