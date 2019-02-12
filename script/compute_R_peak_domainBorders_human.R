# University Paul Sabatier, France
# 11/02/2019


# This R script computes R at CTCF peaks or 3D domain borders.


# Replace with the path of the Github folder.
setwd("PhyloCTCFLooping")


source("script/compDistFun.R")

# Load library
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
assembly="hg19"
data=read.table(paste0("data/CTCF_motif/",assembly,".tsv"),sep='\t',header=T)[,-c(1:2)]

# Compute R for all CTCF motifs, R=1.28
CDF=compDistFun(data,assembly,threshold,dist,peak.GR=NULL,region.GR=NULL,testPalindromic=F)
print(CDF$RVec)

# Compute R for CTCF motifs + CTCF peaks, R=1.71 
CTCF.GR=import.bed("data/CTCF_peak/CTCF_Sydh_GM12878_hg19.bed")
CDFpeak=compDistFun(data,assembly,threshold,dist,peak.GR=CTCF.GR,region.GR=NULL,testPalindromic=F)
print(CDFpeak$RVec)

# Compute R for CTCF motifs located within Arrowhead domain borders, R=5.67 for 40kb borders (best R for 40kb)
ArrowDomborder.GR=import.bed("data/3DDomainBorder/ArrowheadDomBorder_GM12878_40kb_hg19.bed")
CDFborder=compDistFun(data,assembly,threshold,dist,peak.GR=ArrowDomborder.GR,region.GR=NULL,testPalindromic=F)
print(CDFborder$RVec)


