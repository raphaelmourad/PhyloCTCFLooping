# University Paul Sabatier, France
# 11/02/2019


# This R script assesses difference of 3DR between two species genomes.


# Replace with the path of the Github folder.
setwd("PhyloCTCFLooping")

# Function to compute the null distribution of the difference of 3DR values between two species
# This function uses permutations to compute the null distribution
compDiff3DRPerm<-function(ass1,ass2){
 matPMt=matPM[matPM[,3]==ass1|matPM[,3]==ass2,]
 idxS1=sample(1:nrow(matPMt),sum(matPM[,3]==ass1))
 idxS2=setdiff(1:nrow(matPMt),idxS1)
 m1=matPMt[idxS1,]
 m2=matPMt[idxS2,]
 R1=median(m1[m1[,2]=="pm",1])/median(m1[m1[,2]=="mp",1])
 R2=median(m2[m2[,2]=="pm",1])/median(m2[m2[,2]=="mp",1])
 return(R2-R1)
}


# Load precomputed 3DR for different species
file_matPvalPM=paste0("results/matPvalPM/matPvalPM_q0.7_dist100.csv")
matPvalPM=read.table(file=file_matPvalPM,sep='\t',header=T)

# Load precomputed distances between consecutive motifs depending on orientation
file_matPM=paste0("results/matPM/matPM_q0.7_dist100.csv.gz")
matPM=read.csv(file=gzfile(file_matPM),sep='\t',header=T)

# Compute pairwise test between 2 species (here human and mouse genomes)
ass1="hg38" # Human genome
ass2="mm10" # Mouse genome
numPerm=1e4 # number of permutations (the higher, the more accurate the p-value)
diff3DRexp=sapply(1:numPerm,function(x){compDiff3DRPerm(ass1,ass2)})
diff3DRobs=matPvalPM[matPvalPM[,1]==ass2,2]-matPvalPM[matPvalPM[,1]==ass1,2]
pvaldiff3DR=sum(abs(diff3DRexp)>abs(diff3DRobs))/length(diff3DRexp)

# Print results
print(pvaldiff3DR)



