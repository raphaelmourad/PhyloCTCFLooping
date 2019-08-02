# University Paul Sabatier, France
# 11/02/2019


# This R script computes 3DR for different species.


# Replace with the path of the Github folder.
setwd("PhyloCTCFLooping")

# Library to load.
library(GenomicRanges)
library(circlize)
library(bootstrap)

# Function to compute R.
source("script/comp3DR.R")

# Parameters
threshold=0.7 # Quantile threshold on motif score.
dist=100 # Minimum distance between consecutive motifs.

# List motif tsv files. Motif positions were computed using FIMO with default parameters.
lf=list.files("data/CTCF_motif")

# Compute 3DR.
# Store results in matrices "matPvalPM" and "matPM".
matPvalPM=NULL
matPM=NULL
for(i in 1:length(lf)){
 # Assembly name
 assemblyi=strsplit(lf[i],'[.]')[[1]][1]

 # Load motif binding
 datai=read.table(paste0("data/CTCF_motif/",lf[i]),sep='\t',header=T)[,-c(1:2)]

 # Compute R
 CDFi=comp3DR(datai,assemblyi,threshold,dist)

 # Store results
 if(length(CDFi)>1){
  matPvalPM=rbind(matPvalPM,data.frame(Assembly=assemblyi,t(CDFi$RVec))) 
 }

 print(assemblyi)
}

# Compute FDR q value
matPvalPM$X3DR.qvalue=p.adjust(matPvalPM$X3DR.pvalue,"fdr")

# Species names
corresp=read.table("data/corresp_species_assemb.csv",sep='\t',header=T)
matPvalPM=merge(matPvalPM,corresp,by.x="Assembly",by.y="Assembly")


# Print results
print(matPvalPM)














