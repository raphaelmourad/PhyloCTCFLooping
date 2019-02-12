# University Paul Sabatier, France
# 11/02/2019


# This R script computes R for different chromosome regions.


# Replace with the path of the Github folder.
setwd("PhyloCTCFLooping")


source("script/compDistFun.R")


library(GenomicRanges)
library(circlize)
library(bootstrap)
library(rtracklayer)


# Correspondance between assembly and specie's name
corresp=read.table("data/corresp_species_assemb.csv",sep='\t',header=T)

# Parameters
threshold=0.7 # Quantile threshold on motif score. As used in the article.
dist=100 # Minimum distance between consecutive motifs. As used in the article.

# A/B compartments
compA.GR=import.bed("data/region/compAB/compA_GM12878_hg38_liftover.bed")
compB.GR=import.bed("data/region/compAB/compB_GM12878_hg38_liftover.bed")

# Subcomp
subcompAB=read.table("data/region/subcomp/subcomp_GM12878_hg38_liftover.bed", header=F)
subcompAB.GR=GRanges(subcompAB[,1],IRanges(subcompAB[,2],subcompAB[,3]),type=subcompAB[,4])
subcompAB.GR=subcompAB.GR[!is.na(subcompAB.GR$type)]
subcompA1.GR=subcompAB.GR[subcompAB.GR$type=="A1"]
subcompA2.GR=subcompAB.GR[subcompAB.GR$type=="A2"]
subcompB1.GR=subcompAB.GR[subcompAB.GR$type=="B1"]
subcompB2.GR=subcompAB.GR[subcompAB.GR$type=="B2"]
subcompB3.GR=subcompAB.GR[subcompAB.GR$type=="B3"]
subcompB4.GR=subcompAB.GR[subcompAB.GR$type=="B4"]

# Isochores
isochores.GR=import.bed("data/region/isochore/isochores_isoSegmenter_hg38.bed")
isochoresH1.GR=isochores.GR[isochores.GR$name=="H1"]
isochoresH2.GR=isochores.GR[isochores.GR$name=="H2"]
isochoresH3.GR=isochores.GR[isochores.GR$name=="H3"]
isochoresL1.GR=isochores.GR[isochores.GR$name=="L1"]
isochoresL2.GR=isochores.GR[isochores.GR$name=="L2"]


# Make list
regionList=list(compA=compA.GR,compB=compB.GR,
	subcompA1=subcompA1.GR,subcompA2=subcompA2.GR,subcompB1=subcompB1.GR,
	subcompB2=subcompB2.GR,subcompB3=subcompB3.GR,subcompB4=subcompB4.GR,
	isochoresL1=isochoresL1.GR,isochoresL2=isochoresL2.GR,isochoresH1=isochoresH1.GR,
	isochoresH2=isochoresH2.GR,isochoresH3=isochoresH3.GR)
regions=c("compA","compB","subcompA1","subcompA2","subcompB1",
	"subcompB2","subcompB3","subcompB4","isochoresL1","isochoresL2",
	"isochoresH1","isochoresH2","isochoresH3")


# Compute R for different chromosome regions
data=read.table(paste0("data/CTCF_motif/hg38.tsv"),sep='\t',header=T)[,-c(1:2)]
matPvalPM=NULL
for(i in 1:length(regionList)){
 CDFi=compDistFun(data,assembly,threshold,dist,region.GR=regionList[[i]],testPalindromic=F)
 if(length(CDFi)>1){
  matPvalPM=rbind(matPvalPM,data.frame(Region=regions[i],t(CDFi$RVec)))
 }
}


# Print results
print(matPvalPM)


