# University Paul Sabatier, France
# 11/02/2019


# Function to compute R ratio with p-values.
compDistFun<-function(dat,assembly,threshold,dist,peak.GR=NULL,region.GR=NULL,testPalindromic=F){

 if(is.null(assembly)){assembly=NA}

 # Filter and sort dat
 dat=dat[apply(dat,1,function(x){sum(is.na(x))==0}),]
 dat=dat[order(as.character(dat[,1])),]

 # Chromosome info
 if(!is.na(assembly)){ # Read from UCSC
  chrominfoi=read.chromInfo(chromInfo = system.file(package = "circlize","extdata", "chromInfo.txt"), species = assembly, chromosome.index = NULL, sort.chr = TRUE)
  chrleni=chrominfoi$chr.len[names(chrominfoi$chr.len)%in%unique(as.character(dat[,1]))]
  chrleni=chrleni[order(names(chrleni))]
  seqinfoi=Seqinfo(names(chrleni), seqlengths=chrleni, isCircular=NA, genome=assembly)
  asslen=sum(seqlengths(seqinfoi))
 }else{
  asslen=NA
 }

 # Motif GRanges
 motif.GR=GRanges(as.character(dat[,1]),IRanges(dat[,2],dat[,3]),strand=as.character(dat[,4]),score=dat[,5])
 motif.GR=sort(motif.GR,by = ~ seqnames + start)
 if(!is.na(assembly)){
  seqinfo(motif.GR)=seqinfoi
  motif.GR=motif.GR[motif.GR$score>quantile(motif.GR$score,threshold)]
  #motif.GR=motif.GR[motif.GR$score>15]
  motif.GR=motif.GR[seqnames(motif.GR)%in%seqnames(seqinfo(motif.GR))[seqlengths(seqinfo(motif.GR))>10000]] # remove small contigs (length<10000pb)
 }else{
  motif.GR=motif.GR[motif.GR$score>quantile(motif.GR$score,threshold)]
  seqinfoi=NULL
 }

 # Test palindrome
 if(testPalindromic){ 
  motifns.GR=motif.GR
  strand(motifns.GR)='*'
  palindromicValue=(length(unique(motif.GR))/length(unique(motifns.GR)))-1
 }else{palindromicValue=NA}

 # Overlap with peaks
 if(length(peak.GR)>0){
  olPeak=findOverlaps(motif.GR,peak.GR)
  motif.GR=motif.GR[unique(queryHits(olPeak))]
 }

 # Distance motif matrix
 distMotif=data.frame(chr1=seqnames(motif.GR[1:(length(motif.GR)-1)]),chr2=seqnames(motif.GR[2:length(motif.GR)]),pos1=start(motif.GR[1:(length(motif.GR)-1)]),pos2=start(motif.GR[2:length(motif.GR)]),strand1=strand(motif.GR[1:(length(motif.GR)-1)]),strand2=strand(motif.GR[2:length(motif.GR)]),score1=motif.GR[1:(length(motif.GR)-1)]$score,score2=motif.GR[2:length(motif.GR)]$score)
 distMotif$weight=distMotif$score1*distMotif$score2
 distMotif=distMotif[distMotif$chr1==distMotif$chr2,]
 distMotif=distMotif[(distMotif$pos2-distMotif$pos1)>dist,]

 # Analysis depending on genomic regions
 if(length(region.GR)>0){ 
  distMotif.GR=GRanges(distMotif[,1],IRanges(distMotif[,3],distMotif[,4]),strand1=distMotif[,5],strand2=distMotif[,6])
  olRegion=findOverlaps(distMotif.GR,region.GR,type="within")
  distMotifRegion.GR=distMotif.GR[queryHits(olRegion)]
  distMotif=data.frame(chr1=seqnames(distMotifRegion.GR),chr2=seqnames(distMotifRegion.GR),
	pos1=start(distMotifRegion.GR),pos2=end(distMotifRegion.GR),
	strand1=distMotifRegion.GR$strand1,strand2=distMotifRegion.GR$strand2)
 }

 # Distance depending of motif orientation
 distMotifpp=distMotif[distMotif$strand1=='+'&distMotif$strand2=='+',]
 distMotifmm=distMotif[distMotif$strand1=='-'&distMotif$strand2=='-',]
 distMotifpm=distMotif[distMotif$strand1=='+'&distMotif$strand2=='-',]
 distMotifmp=distMotif[distMotif$strand1=='-'&distMotif$strand2=='+',]
 distpp=distMotifpp$pos2-distMotifpp$pos1
 distmm=distMotifmm$pos2-distMotifmm$pos1
 distpm=distMotifpm$pos2-distMotifpm$pos1
 distmp=distMotifmp$pos2-distMotifmp$pos1

 # Assess distance differences
 if(length(distpp)>20 & length(distpm)>20 & length(distmp)>20 & length(distmm)>20){

  dataDist=data.frame(dir=factor(c(rep("pm",length(distpm)),rep("mp",length(distmp)),
	rep("pp",length(distpp)),rep("mm",length(distmm))),levels=c("pm","mp","pp","mm")),
	dist=c(distpm,distmp,distpp,distmm))
  distby=round(as.vector(by(dataDist$dist,dataDist$dir,median)))

  # Wilcoxon tests
  ttpmmp=wilcox.test(distpm,distmp)
  ttppmm=wilcox.test(distpp,distmm)

  # Store distances depending on orientation
  distpmmp=c(distpm,distmp,distpp,distmm)
  chrpmmp=c(as.character(distMotifpm$chr1),as.character(distMotifmp$chr1),
	as.character(distMotifpp$chr1),as.character(distMotifmm$chr1))
  dirpmmp=c(rep("pm",length(distpm)),rep("mp",length(distmp)),
	rep("pp",length(distpp)),rep("mm",length(distmm)))
  distanceMat=data.frame(Distance=distpmmp,Orientation=dirpmmp,
	Assembly=rep(assembly,length(distpmmp)),Chr=chrpmmp)

  # Bootstrap estimation of 95% confidence interval for R
  Rboot=bootstrap(1:length(c(distpm,distmp)),100,func<-function(x,xdata){xd=xdata[x,];median(xd[xd[,2]=="pm",1])/median(xd[xd[,2]=="mp",1])},xdata=data.frame(d=c(distpm,distmp),pm=c(rep("pm",length(distpm)),rep("mp",length(distmp)))))
  Rse=sd(Rboot$thetastar)
  Rvalue=distby[1]/distby[2]
  Rlower=Rvalue-(Rse*1.96)
  Rupper=Rvalue+(Rse*1.96)

  # Output
  RVec=c(Rvalue,Rlower,Rupper,Cvalue=distby[3]/distby[4],R.pvalue=ttpmmp$p.value,C.pvalue=ttppmm$p.value,
	Npairs=nrow(distMotif),Nmotifs=length(unique(c(distMotif$pos1,distMotif$pos2))),
	MedianDistPM=distby[1],MedianDistMP=distby[2],MedianDistPP=distby[3],MedianDistMM=distby[4],
	AssemblyLength=asslen,palindromicValue)
  names(RVec)=c("R","Rlower","Rupper","C","R.pvalue","C.pvalue","Npairs","Nmotifs","MedianDistPM","MedianDistMP",
	"MedianDistPP","MedianDistMM","AssemblyLength","PalindromeValue")
 
  res=list(RVec=RVec,distanceMat=distanceMat)
  return(res)
 }else{
  return("Not enough motifs to compute R.")
 }
}# END OF compDistFun()
