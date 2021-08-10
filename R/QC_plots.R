#!/usr/bin/env Rscript

################  LOAD OPTPARSE
library("optparse")
library("data.table")
library("ggplot2")
library("RColorBrewer")
library("ggrepel")

################ DEFINE THE INPUT PARAMETERS
option_list = list(
  make_option(c("-d", "--dir"), type="character", default=NULL,help="Path to the samples' 'raw_coverage_per_smMIP.txt' and 'filtered_read_counts.txt' files (Need to be in a single folder).", metavar="character"),
  make_option(c("-c", "--code"), type="character", default=getwd(),help="Path to smMIP tools source functions file. If not supplied, it is assumed that the file (smMIPs_Function.R) is located in your working directory", metavar="character"),
  make_option(c("-t", "--type"), type="character", default="smMIP",help="Options are smMIP (default) or gene. If the smMIPs IDs format is Gene_001..., then -t 'gene' can be used. This will color smMIPs targeted the same gene with the same color", metavar="character"))
  
#Takes user input parameters
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser);
if (is.null(opt$dir)){ 
  print_help(opt_parser)
  stop("the -d parameter must be supplied", call.=FALSE)
}

if (is.null(opt$code)){
  source(paste0(getwd(),"/smMIPs_Functions.R"))
} else {
  source(paste0(opt$code,"/smMIPs_Functions.R"))
}

############### LOADING THE DATA AND GENERATE PLOTS
defaultW <- getOption("warn") 
options(warn = -1) 

tab=load_qc_coverage_per_smMIP()   
sm=data.table(melt(tab,id.vars = "smMIPs"))
sm=sm[,value:=mean(value),by=smMIPs]
sm=unique(sm[,variable:=NULL])

if(opt$t=="gene"){
  sm$gene=gsub("_.*","",sm$smMIPs)
  sm=sm[order(gene)]
  n <- length(unique(sm$gene))
  col.par = function(n) sample(seq(0.3, 1, length.out=250),n); 
  cols = rainbow(n, s=col.par(n), v=col.par(n))[sample(1:n,n)]
  
  p1=ggplot(data=sm, aes(x=smMIPs, y=value, fill=gene)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_classic()   + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position="bottom") + theme(legend.title = element_blank()) +guides(fill=guide_legend(nrow=3,byrow=TRUE)) +
  ylab("Mean Coverage") + scale_fill_manual(values = cols)
} else {
  sm=sm[order(smMIPs)]
  p1=ggplot(data=sm, aes(x=smMIPs, y=value)) +
    geom_bar(stat="identity", color="black", position=position_dodge())+
    theme_classic()   + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("Mean Coverage")
}

ggsave(paste0(opt$dir,"/Average_Coverage_per_smMIP.pdf"),p1, width=nrow(sm)/7, height=10, units="in")
write.table(tab,file=paste0(opt$dir,"/Cohort_Coverage_per_smMIP.txt"),col.names = T,row.names = F,quote = F,sep = '\t')
#######

tab=load_qc_filtered.reads()
names(tab)=c("sample_ID","Total Reads", "Low Mapping Score", "Bad Sam Flag", "Unexpected Insert Length", "Hard Clip", "Removed due to Filtered Mate", "Off Target", "Good Reads")
write.table(tab,file=paste0(opt$dir,"/Cohort_Filter_Reason.txt"),col.names = T,row.names = F,quote = F,sep = '\t')

names(tab)=c("sample_ID","Total Reads", "Low Mapping\nScore", "Bad Sam\nFlag", "Unexpected\nInsert Length", "Hard Clip", "Removed due to\nFiltered Mate", "Off Target", "Good Reads")
sm=data.table(melt(as.data.table(tab[,3:8])))
p2=ggplot(sm, aes(x=variable, y=value)) + 
  geom_boxplot(notch=T,outlier.colour="black", outlier.shape=16,
  outlier.size=2) + xlab("Filter") + ylab("Percentage of Reads")  +
  theme_classic()   + theme(axis.text.x = element_text(size=10, angle = 45, vjust = 1, hjust=1, face="bold"), axis.text.y=element_text(size=12,face="bold"),axis.title=element_text(size=14,face="bold")) 

ggsave(paste0(opt$dir,"/Filter_Reason.pdf"),p2, width=6, height=10, units="in")

sm=data.table(melt(as.data.table(tab[,c(2,9)])))
p3=ggplot(sm, aes(x=variable, y=value)) + 
  geom_boxplot(notch=T,outlier.colour="black", outlier.shape=16,
               outlier.size=2) + xlab("") + ylab("No. of Read")  + 
  theme_classic() + theme(axis.text.x = element_text(size=10, angle = 45, vjust = 0.5, face="bold"), axis.text.y=element_text(size=12,face="bold"),axis.title=element_text(size=14,face="bold")) 

ggsave(paste0(opt$dir,"/Number_of_Reads.pdf"),p3, width=6, height=10, units="in")

tab$Percentage=100*tab[,9]/tab[,2]
sm=tab[,c(1,10)]
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}
n=sm$sample_ID
n[!is_outlier(sm$Percentage)]=""
sm$outlier=n
p4=ggplot(sm, aes(x=0, y=Percentage)) + 
  geom_boxplot(notch=T,outlier.colour="black", outlier.shape=16,
               outlier.size=2) + xlab("Good Reads") + ylab("Percentage of Total")  + 
  theme_classic() + theme(axis.ticks.x=element_blank(), axis.text.y=element_text(size=12,face="bold"),axis.title=element_text(size=14,face="bold")) +
  coord_cartesian(ylim = c(0,100)) + 
    ggrepel::geom_text_repel(data = sm, aes(label = outlier),nudge_x=2,segment.size  = 0.2, size=2)+
    scale_x_discrete(labels=c(0,1))
    
ggsave(paste0(opt$dir,"/Outliers.pdf"),p4, width=6, height=10, units="in")

options(warn = defaultW)
