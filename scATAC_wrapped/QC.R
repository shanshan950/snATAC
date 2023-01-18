#! /usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
name=args[1]
genome=args[2]
### To generate mouse/human scatter plot
library(dplyr)
library(scales)
library(ggplot2)
library(gridExtra)
##################################################################################### To generate histogram of fragments counts
x<-read.table(paste("./Demit.SPLIT.",genome,"/bedsummary",sep=""))
thetitle<-paste(name,genome,":",length(which(x$V1>1000)),"Cells (",length(x$V1),"BC")
data.frame(value=log10(x$V1)) %>% ggplot(.)+aes(value)+geom_histogram(binwidth=0.05,color="black")+xlim(1,5)+theme_classic()+ggtitle(thetitle)+geom_vline(xintercept=3,linetype=2)+labs(x="log10 Fragments")+theme(plot.title = element_text( face="bold",hjust = 0.5))+ggtitle(thetitle) -> p.fragmentsN



####################################################################################### To generate fragment size distribution in human and mouse

fragment<-read.table(paste(name,".",genome,".PE.monoclonal.bed",sep=""))
cbind(fragment,size=fragment[,3]-fragment[,2]) %>% ggplot(.)+aes(size)+geom_density()+theme_classic()+ggtitle(paste(name,genome,"frament size"))+theme(plot.title = element_text( face="bold",hjust = 0.5)) ->p.fragment.size
#

######################################################################################### To generate TSS and peak fraction histogram

# Mouse
TSSPeak<-read.table(paste("./Demit.SPLIT.",genome,"/Summary.TSSPeak",sep=""),header=T)
data.frame(TSSPeak,TSSfraction=TSSPeak$FragTSS/TSSPeak$Fragments,PEAKfraction=TSSPeak$FragPeak/TSSPeak$Fragments)-> TSSPeak.ratio
# To generate histogram
p.TSS<-ggplot(TSSPeak.ratio)+aes(TSSfraction)+geom_histogram(color="black",binwidth=0.01)+xlim(0,1)+ggtitle(paste(name,genome,": \nFraction of fragment on TSS per cell"))+theme_classic()+theme(plot.title = element_text( face="bold",hjust = 0.5,size=10))
p.Peak<-ggplot(TSSPeak.ratio)+aes(PEAKfraction)+geom_histogram(color="black",binwidth=0.01)+xlim(0,1)+ggtitle(paste(name,genome,": \nFraction of fragment on Peak per cell"))+theme_classic()+theme(plot.title = element_text( face="bold",hjust = 0.5,size=10))


# Mouse  To generate histogram of fragments on peaks counts in Mouse
thetitle<-paste(name,genome,":",length(which(TSSPeak$FragPeak>1000)),"Cells with 1000 onPEAK")
data.frame(value=log10(TSSPeak$FragPeak)) %>% ggplot(.)+aes(value)+geom_histogram(binwidth=0.05,color="black")+xlim(1,5)+theme_classic()+ggtitle(thetitle)+geom_vline(xintercept=3,linetype=2)+labs(x="log10 Fragments")+theme(plot.title = element_text( face="bold",hjust = 0.5))+ggtitle(thetitle) -> p.fragmentsNonPEAK





#######################################################################Reads QC
RawQC<-read.table("RawQC",header=T)
RawQC<-as.data.frame(t(RawQC))
RawQC<-data.frame(Step=row.names(RawQC),Reads=RawQC$V1)
RawQC$Step<-factor(RawQC$Step,levels=c("Rawreads",paste("Uniqmap",genome,sep=""),paste("Demit",genome,sep="")))
p.rawmap<-ggplot(RawQC)+aes(Step,Reads,label=comma(Reads),fill=Step)+geom_bar(stat="identity",color="black")+geom_text(size=6,vjust=1.5)+ggtitle(paste(name," Raw reads mapping"))+theme_classic()+theme(plot.title = element_text( face="bold",hjust = 0.5,size=15),axis.text=element_text(size=14),axis.title=element_blank(),legend.text=element_text(size=12),legend.title=element_blank())+scale_fill_brewer(palette="Set2")

###mousesummary
Summary<-read.table(paste("./Demit.SPLIT.",genome,"/Summary",sep=""))

barcoderate<-data.frame(Step=c(paste("Demit",genome,sep=""),"Barcoded"),Reads=c(RawQC[2,2],sum(Summary$V2)))
barcoderate$Step<-factor(barcoderate$Step,levels=c(paste("Demit",genome,sep=""),"Barcoded"))
p.barcoderate<-ggplot(barcoderate)+aes(Step,Reads,label=comma(Reads),fill=Step)+geom_bar(stat="identity",color="black")+geom_text(size=6,vjust=1.5)+ggtitle(paste(name," Uniq maped barcoded",genome))+theme_classic()+theme(plot.title = element_text( face="bold",hjust = 0.5,size=15),axis.text=element_text(size=14),axis.title=element_blank(),legend.text=element_text(size=12),legend.title=element_blank())+scale_fill_brewer(palette="Set2")


DedupDemito<-data.frame(Step=c("mean.total","mean.Dedup"),Reads=c(round(mean(Summary$V2),digits=0),round(mean(Summary$V3),digits=0)))
DedupDemito$Step<-factor(DedupDemito$Step,levels=c("mean.total","mean.Dedup"))
thetitle<-paste(nrow(Summary),"Detected barcode combination")

p.DedupDemito<-ggplot(DedupDemito)+aes(Step,Reads,label=comma(Reads),fill=Step)+geom_bar(stat="identity",color="black")+geom_text(size=6,vjust=1.5)+ggtitle(thetitle)+theme_classic()+theme(plot.title = element_text( face="bold",hjust = 0.5,size=15),axis.text=element_text(size=14),axis.title=element_blank(),legend.text=element_text(size=12),legend.title=element_blank())+scale_fill_brewer(palette="Set1")







pdf(paste(name,".QC.pdf",sep=""))
#print(p.species)
print(p.fragmentsN)
print(p.fragmentsNonPEAK)
print(p.fragment.size)
grid.arrange(p.TSS,p.Peak)
dev.off()


pdf(paste(name,".Reads.QC.pdf",sep=""),width=8,height=8)
print(p.rawmap)
grid.arrange(p.barcoderate,p.DedupDemito)
dev.off()
