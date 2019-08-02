# University Paul Sabatier, France
# 11/02/2019


# This R script computes different 3DR ratios from the human genome: 
# - the standard 3DR when only motifs are provided; 
# - the 3DRp when predicted ChIP-seq peaks are provided;
# - the 3DRc when conservation scores are provided;



# Replace with the path of the Github folder.
setwd("PhyloCTCFLooping")


source("script/comp3DR.R")


library(GenomicRanges)
library(circlize)
library(bootstrap)
library(rtracklayer)


# Correspondance between assembly and specie's name
corresp=read.table("data/corresp_species_assemb.csv",sep='\t',header=T)

# Parameters
threshold=0.7 # Quantile threshold on motif score. As used in the article.
dist=100 # Minimum distance between consecutive motifs. As used in the article.

# Assembly
assembly="hg38"
data=read.table(paste0("data/CTCF_motif/",assembly,".tsv"),sep='\t',header=T)[,-c(1:2)]

# Compute 3DR for all CTCF motifs, R=1.28
CDF=comp3DR(data,assembly,threshold,dist,peak.GR=NULL,region.GR=NULL,testPalindromic=F)
print(CDF$RVec)

# Compute 3DRp for CTCF motifs + predicted CTCF peaks, Rp=1.44
CTCFpred=read.table(paste0("data/CTCF_deepbind/",assembly,"_deepbindpred.bed"),sep='\t',header=F)
CTCFpred=CTCFpred[CTCFpred[,7]>quantile(CTCFpred[,7],0.7,na.rm=T),]
CTCFpred=CTCFpred[!is.na(CTCFpred[,3]),]
CTCFpred.GR=GRanges(CTCFpred[,2],IRanges(CTCFpred[,3],CTCFpred[,4]))
CDFpred=comp3DR(data,assembly,threshold,dist,peak.GR=CTCFpred.GR)
print(CDFpred$RVec)

# Compute 3DRc for CTCF motifs + conservation, Rc=1.64
cons.GR=import.bed(paste0("data/CTCF_cons/",assembly,"_cons_50b.bed"))
cons.GR=cons.GR[cons.GR$score>quantile(cons.GR$score,0.7)] 
CDFcons=comp3DR(data,assembly,threshold,dist,peak.GR=cons.GR)
print(CDFcons$RVec)















D
