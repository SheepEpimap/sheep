#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args)<4) stop("Usage: 17_fdr_adjust.R SAMPLE INPUT OUTPUT_DIR FDR_DETAIL_DIR")
sample_id<-args[[1]]; count_file<-args[[2]]; output_dir<-args[[3]]; fdr_detail_dir<-args[[4]]
suppressPackageStartupMessages(library(data.table))
if (!file.exists(count_file) || file.info(count_file)$size == 0) {
  warning(paste(sample_id, "has an empty ASE count table")); quit(status=0)
}
dir.create(output_dir,showWarnings=FALSE,recursive=TRUE); dir.create(fdr_detail_dir,showWarnings=FALSE,recursive=TRUE)
d<-fread(count_file,sep="\t",header=FALSE,check.names=FALSE)
if(ncol(d)!=8L) stop(sprintf("Expected 8 columns; observed %d",ncol(d)))
setnames(d,c("contig","start","end","ref","alt","refCount","altCount","totalCount"))
d[,c("start","end","refCount","altCount","totalCount"):=lapply(.SD,as.numeric),.SDcols=c("start","end","refCount","altCount","totalCount")]
d<-d[is.finite(totalCount)&totalCount>0&refCount>=3&altCount>=3]
d[,ratio:=pmin(refCount,altCount)/totalCount]; d<-d[ratio>=0.1]
if(nrow(d)==0L){warning(paste(sample_id,"has no sites after QC"));quit(status=0)}
d[,position:=start+1]; d[,variantID:=paste(contig,position,sep="_")]
d[,p_value:=mapply(function(x,n)binom.test(x=x,n=n,p=0.5,alternative="two.sided")$p.value,refCount,totalCount)]
d[,padj_bh:=p.adjust(p_value,"BH")]; d[,padj_holm:=p.adjust(p_value,"holm")]; d[,padj_hochberg:=p.adjust(p_value,"hochberg")]
cutoffs<-c(0.05,0.1,0.2,0.3,0.4,0.5)
sm<-rbindlist(lapply(cutoffs,function(z)data.table(target_fdr=z,sig_bh=sum(d$padj_bh<=z),sig_holm=sum(d$padj_holm<=z),sig_hochberg=sum(d$padj_hochberg<=z))))
res<-d[,.(contig,position,end,ref,alt,refCount,altCount,totalCount,variantID,ratio,p_value,padj_bh,padj_holm,padj_hochberg)]
fwrite(res,file.path(output_dir,paste0(sample_id,"_binomial_results.txt")),sep="\t"); fwrite(sm,file.path(fdr_detail_dir,paste0(sample_id,"_fdr_summary.txt")),sep="\t")
cat(sprintf("%s: after QC=%d; BH FDR<=0.1=%d\n",sample_id,nrow(res),sum(res$padj_bh<=0.1)))
