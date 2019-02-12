# University Paul Sabatier, France
# 11/02/2019


# This R script computes ancestral R reconstruction using precomputed R from different species and a phylogenetic tree.


# Replace with the path of the Github folder.
setwd("PhyloCTCFLooping")

# Library to load.
library(phytools)

# Parameters
mode="_cons" # "" "_deepbind" "_cons"
spec="mammals" # "mammals" or "vertebrates" 

# Load R for different species (precomputed)
file_matPvalPM=paste0("results/matPvalPM/matPvalPM_q0.7_dist100",mode,".csv")
matPvalPM=read.table(file=file_matPvalPM,sep='\t',header=T)

# Filters genomes with enough motif pairs
matPvalPM=matPvalPM[matPvalPM$Npairs>8000,]

# Select species
if(spec=="mammals"){
 mammals=as.character(read.table("data/mammals.csv",sep='\t',header=T)[,1])
 matPvalPM=matPvalPM[matPvalPM$Name%in%mammals,]
}else if(spec=="vertebrates"){
 vertebrates=as.character(read.table("data/vertebrates.csv",sep='\t',header=T)[,1])
 matPvalPM=matPvalPM[matPvalPM$Name%in%vertebrates,]
}
assnonum=gsub("[[:digit:]]","",as.character(matPvalPM[,1]))

# Load tree and parse it
treeall=read.newick("data/tree/hg38.100way.nh")
label=treeall$tip.label
labelnonum=gsub("[[:digit:]]","",label)
tree=drop.tip(treeall,tip=treeall$tip.label[!labelnonum%in%assnonum])
tree$node.label=(length(tree$tip.label)+1):(length(tree$tip.label)*2-1)
tree$tip.label=as.character(matPvalPM$Assembly[match(labelnonum,assnonum)][!is.na(matPvalPM$Assembly[match(labelnonum,assnonum)])])

# Match species in "matPvalPM" with species in the tree
matchTreeMat=matPvalPM[match(labelnonum,assnonum),]
matchTreeMat=matchTreeMat[!is.na(matchTreeMat[,1]),]

# Ancestral R reconstruction 
stat=matchTreeMat$R
names(stat)=tree$tip.label
ancML=anc.ML(tree,x=stat)

# Plot ancestral character reconstruction
dir.create("results/phylo/")
file_plottree_FC=paste0("results/phylo/ancestral_R_reconstruction_",spec,mode,".pdf")
pdf(file_plottree_FC,10,12)
plotANC=contMap(tree,x=stat,method="anc.ML",col="blue",type="phylogram")
dev.off()

# Mantel test to assess if R is evolutionary conserved
distSpecies=dist.nodes(tree)[1:length(tree$tip.label),1:length(tree$tip.label)]
statDiff=sapply(stat,function(x){abs(x-stat)})
mantel.test(distSpecies,statDiff,nperm=1e5)




