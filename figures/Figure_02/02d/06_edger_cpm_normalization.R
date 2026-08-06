#!/usr/bin/env Rscript
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
library(edgeR)

args = commandArgs(trailingOnly=TRUE)
data <- read.table(args[1], sep=",", row.names=1, header=T)
RG <- DGEList(counts=data, group=rep(1,ncol(data)))
RG <- calcNormFactors(RG)
RG <- cpm(RG)
write.table(RG, file=args[2], sep=",", quote=F)
