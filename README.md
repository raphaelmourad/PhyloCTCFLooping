# Studying 3D genome evolution using genomic sequence

![alt text](https://github.com/morphos30/PhyloCTCFLooping/blob/master/approach.png)

**Overview**

We propose a novel approach to study the 3D genome evolution in vertebrates using the genomic sequence only, e.g. without the need for Hi-C data. The approach is simple and relies on comparing the distances between convergent and divergent CTCF motifs (ratio 3DR).

**Systems Requirements**

The scripts were written in R language. 

To run the scripts, you need several R packages. To install the packages:
`install.packages(c("circlize","bootstrap","phytools"))` \
`source("https://bioconductor.org/biocLite.R")` \
`biocLite("GenomicRanges")` 

**Usage**

Before computing ratio 3DR, one must use FIMO to scan for CTCF motifs. FIMO outputs a "tsv" file that is used by the "comp3DR.R" function to compute 3DR. You can use FIMO online (http://meme-suite.org/tools/fimo), select the genome assembly and upload the CTCF MEME file "data/CTCF_meme/MA0139.1.meme".

You can also use precomputed CTCF motifs from FIMO for selected assemblies: hg19, hg38, bosTau8, ce11, dm6, mm10, rn6 and xenTro7. These motifs are available in the folder "data/CTCF_motif".

In this package, there are three main folders: 
- The folder "data" contains: the CTCF motif PWM in MEME format for FIMO (subfolder "CTCF_meme"), CTCF motifs called by FIMO (subfolder "CTCF_motif"), CTCF motif DeepBind scores (subfolder "CTCF_deepbind"), CTCF motif conservation scores (subfolder "CTCF_cons"), CTCF ChIP-seq peaks from GM12878 ENCODE cells (subfolder "CTCF_peak"), 3D domain borders (subfolder "3DDomainBorder"), chromosome regions (subfolder "region"), phylogenetic tree (subfolder "tree").
- The folder "script" contains seven R scripts: "comp3DR.R" is the function to compute 3DR, "compute_3DR_species.R" to compute 3DR in different species, "compute_3DR_3DRp_3DRc_human.R" to compute different 3DR's in human, "compute_3DR_regions_human.R" to compute 3DR for different chromosome regions, "compute_3DR_peak_domainBorders_human.R" to compute 3DR at CTCF peaks and 3D domain borders, "test_difference3DR_2species.R" to test the difference of 3DR values between two species and "ancestral_3DR_reconstruction.R" to use ancestral 3DR reconstruction. 
- The folder "results" contains three subfolders: precomputed 3DR for different species (subfolder "matPvalPM"), DNA motif prediction results (subfolder "predMotif"), precomputed distances between consecutive motifs depending on orienation (subfolder "matPM") and ancestral 3DR reconstruction results (subfolder "phylo"). 

**References**
Raphael Mourad. Studying 3D genome evolution using genomic sequence. BioRxiv, May 23, 2019.
https://www.biorxiv.org/content/10.1101/646851v1.abstract

**Contact**:
raphael.mourad@ibcg.biotoul.fr
raphael.mourad@univ-tlse3.fr
