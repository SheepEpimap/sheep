#!/usr/bin/env Rscript
# Optional slow simulation retained from the source Markdown.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)<3)stop("Usage: 18_fdr_simulation.R SAMPLE INPUT OUTPUT_DIR [SIMULATIONS] [CORES]")
sample_id<-args[[1]]; count_file<-args[[2]]; output_dir<-args[[3]]; n_sim<-if(length(args)>=4)as.integer(args[[4]])else 1000L; n_cores<-if(length(args)>=5)as.integer(args[[5]])else 6L
suppressPackageStartupMessages(library(data.table))
if (!file.exists(count_file) || file.info(count_file)$size == 0) {
  warning(paste(sample_id, "has an empty ASE count table")); quit(status=0)
}; suppressPackageStartupMessages(library(parallel)); dir.create(output_dir,showWarnings=FALSE,recursive=TRUE)
d<-fread(count_file,header=FALSE,sep="\t"); if(ncol(d)!=8L)stop("Input must contain 8 columns"); setnames(d,c("contig","start","end","ref","alt","refCount","altCount","totalCount")); d[,c("refCount","altCount","totalCount"):=lapply(.SD,as.numeric),.SDcols=c("refCount","altCount","totalCount")]; d<-d[totalCount>0&refCount>=3&altCount>=3]; d[,ratio:=pmin(refCount,altCount)/totalCount]; d<-d[ratio>=0.1]; if(!nrow(d))quit(status=0)
actual_p<-sort(mapply(function(x,n)binom.test(x,n,p=0.5)$p.value,d$refCount,d$totalCount)); thresholds<-c(0.05,0.1,0.2,0.3,0.4,0.5)
worker<-function(target){vals<-replicate(n_sim,{sim_ref<-rbinom(nrow(d),size=d$totalCount,prob=0.5);sim_p<-sort(mapply(function(x,n)binom.test(x,n,p=0.5)$p.value,sim_ref,d$totalCount));cand<-rev(actual_p);ok<-vapply(cand,function(p){obs<-sum(actual_p<=p);false<-sum(sim_p<=p);obs>0&&false/obs<=target},logical(1));if(any(ok))cand[which(ok)[1]]else 0});cutoff<-mean(vals,na.rm=TRUE);data.table(target_fdr=target,sim_p_threshold=cutoff,significant_sites=sum(actual_p<=cutoff))}
cores<-max(1L,min(n_cores,detectCores()-1L,length(thresholds))); res<-if(cores>1L)rbindlist(mclapply(thresholds,worker,mc.cores=cores))else rbindlist(lapply(thresholds,worker)); fwrite(res,file.path(output_dir,paste0(sample_id,"_simulation_fdr_results.txt")),sep="\t")
