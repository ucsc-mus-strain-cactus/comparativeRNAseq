"""
Run HTSeq-count followed by DESeq2
"""

find aligned_reference -type f -name "*Coord*bam" | xargs -n 1 -P 20 -I{} sh -c  "htseq-count -m union -r pos -i gene_id -a 10 --stranded=no -f bam '{}' /hive/users/ifiddes/ihategit/pipeline/gencode_vm8/C57B6J.gtf > '{}'.htseq.counts"

find aligned_new_reference -type f -name "*Coord*bam" | xargs -n 1 -P 20 -I{} sh -c  'n="{}"; g=`dirname $n | cut -d "/" -f 2`; ref=STAR_references/${g}/${g}.gtf; htseq-count -m union -r pos -i gene_id -a 10 --stranded=no -f bam "{}" ${ref} > "{}".htseq.counts'


# use this to generate new references
find . -mindepth 1 -maxdepth 1 -type d | xargs -n 1 -P 5 -I {} sh -c 'n={}; g=${n:2}; cd ${g} && STAR --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles ${g}.fa --sjdbGTFfile ${g}.gtf'

