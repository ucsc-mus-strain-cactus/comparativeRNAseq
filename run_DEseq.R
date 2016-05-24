#!/usr/bin/env Rscript
library("DESeq2")
args <- commandArgs(TRUE)
genome <- args[1]
out_path <- args[2]
counts_path <- args[3]
data_path <- paste0(counts_path, '/', genome, '.counts.txt')
conds_path <- paste0(counts_path, '/', genome, '.conds.txt')
conds <- read.table(conds_path, sep="\t", header=T, row.names=1)
data <- read.table(data_path, sep="\t", header=T, row.names=1)
dds <- DESeqDataSetFromMatrix(countData=data, colData=conds, design=~condition)
dds$condition <- relevel(dds$condition, ref='ref')
dds <- DESeq(dds)
res <- results(dds)
pdf(paste0(out_path, '/', genome, '.maplot.pdf'))
plotMA(res, main=genome)
dev.off()
p <- paste0(out_path, '/', genome, '.results.rData')
save(dds, file=p)
# if none are significant, this will fail
baseMeanPerLvl <- try(sapply( levels(dds$condition), function(lvl) rowMeans( counts(dds,normalized=TRUE)[,dds$condition == lvl] ) ) )
if (class(baseMeanPerLvl) == 'try-error'){
    baseMeanPerLvl <- counts(dds, normalized=TRUE)
}
p2 <- paste0(out_path, '/', genome, '.data_table.tsv')
q <- cbind(baseMeanPerLvl, res$padj)
colnames(q) <- c('ref', 'tgt', '')
write.table(q, p2)

