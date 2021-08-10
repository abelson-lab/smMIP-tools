#!/usr/bin/env Rscript
################  LOAD OPTPARSE
library("optparse")

################ DEFINE THE INPUT PARAMETERS
option_list = list(
  make_option(c("-b", "--bam.file"), type="character", default=NULL,help="Path to the filtered bam file output by map_smMIPs_extract_UMIs.R (sample_clean.bam) [MENDATORY]", metavar="character"),
  make_option(c("-p", "--panel.file"), type="character", default=NULL,help="Path to the smMIP design file [MENDATORY]", metavar="character"),
  make_option(c("-s", "--sample.name"), type="character", default=NULL,help="Sample ID that will be used to name the output file(s). If not provided, the name of the folder containing the bam file is assumed to be the sample name", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,help="Path for the output pileup file(s). If not supplied, the output will be saved within the folder that contain the bam file", metavar="character"),
  make_option(c("-O", "--tmp.output"), type="character", default=NULL,help="Path for temporary files. If not supplied, the tmp.output will be the same as defined by the -o option. When working on HPC, supplying a local folder to write and read temporary files can increase speed", metavar="character"),
  make_option(c("-c", "--code"), type="character", default=getwd(),help="Path to smMIP tools source functions file. If not supplied, it is assumed that the file (smMIPs_Function.R) is located in your working directory", metavar="character"),
  make_option(c("-d", "--mnd"), type="integer", default=1,help="Minimum depth to consider in the pileup. ", metavar="character"),
  make_option(c("-m", "--mmq"), type="integer", default=50,help="Minimum mapping quality to consider in the pileup. ", metavar="character"),
  make_option(c("-q", "--mbq"), type="integer", default=10,help="Minimum base quality to consider in the pileup. ", metavar="character"),
  make_option(c("-r", "--rank"), type="character", default="F",help="Options are 'F' or 'T'. Reporting all the possible alleles (A,C,T,G,-,+) if they observed. When 'T', only the allele with the most read support (other than the reference) at each genomic postion will be reported. This option can help to remove many errors from the data yet it comes with the risk of missing real mutations", metavar="character"),
  make_option(c("-u", "--umi"), type="character", default="T",help="Options are 'T' or 'F'. If 'T', the UMI sequences that are associated with each allele will be reported. This information is required to estimate read-to-sample misassignment", metavar="character"),
  make_option(c("-f", "--family.size"), type="integer", default=0,help="The minimum read-family size to consider in single strand consensus pileups", metavar="character"),
  make_option(c("-v", "--consensus.cutoff"), type="numeric", default=0.7,help="The percentage of reads that must have the same base to form a consensus", metavar="character"),
  make_option(c("-t", "--threads"), type="integer", default=1,help="The number of cores to use for parallel processing", metavar="character"))


#takes user input parameters
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser);
if (is.null(opt$bam.file) | is.null(opt$panel.file)){
  print_help(opt_parser)
  stop("The -b and -p parameters must be supplied", call.=FALSE)
}

if (is.null(opt$code)){ #load source functions
  source(paste0(getwd,"/smMIPs_Functions.R"))
} else {
  source(paste0(opt$code,"/smMIPs_Functions.R"))
}

if (is.null(opt$output)){ 
  opt$output=dirname(opt$bam.file)
} else if (!dir.exists(opt$output)){
      dir.create(opt$output)
  }

if (is.null(opt$tmp.output)){ 
  opt$tmp.output=opt$output
  del=0
} else {
  dir.create(opt$tmp.output)
  del=1
}

if (is.null(opt$sample.name)){ 
  opt$sample.name=basename(dirname(opt$bam.file))
}

if (opt$consensus.cutoff<0.51){ 
  print(paste0("The consensus cutoff must be at the range of 0.51-1, yet you entered ",opt$consensus.cutoff,". It is now set by default to 0.7"),quote = F)
  opt$consensus.cutoff=0.7
}

################  LOAD LIBRARIES
library("Rsamtools")
library("data.table")


#Print all the parameters
cat("###############################\n")
cat("        Run Parameters\n")
cat("###############################\n")
print(opt)
cat("###############################\n")
cat("            Running...\n")
cat("###############################\n")


################ CREATE A DATA OBJECT
data<-list()
cat("Loading bam\n")
data$samtable<-as.data.table(load.bam(opt$bam.file)) 
data$panel<-load.panel(opt$panel.file)
data$samtable$smmip=gsub(".*[||]","",data$samtable$qname)
data$samtable$umi=gsub(".*:","",data$samtable$smmip)
data$samtable$smmip=gsub(":.*","",data$samtable$smmip)

################ DEFINE THE PILEUP PARAMETERS 
p_param <- PileupParam(max_depth=10000000, 
                       min_nucleotide_depth=opt$mnd, 
                       distinguish_strand=T,
                       min_mapq=opt$mmq,
                       min_base_quality=opt$mbq, #Default is 10 since high values may bias variant allele frequencies.
                       distinguish_nucleotides = T,
                       ignore_query_Ns = F,
                       include_deletions = T,
                       include_insertions = T)


################ CREATING HEADER
u=unique(data$samtable$rname)
header=as.data.frame(matrix(nrow = length(u)+2,ncol=11))
header[1:(length(u)+2),]=""
colnames(header)=c("qname","flag","rname","pos","mapq","cigar","mrnm","mpos","isize","seq","qual")
header[,1]=c(rep("@SQ",length(u)),"@RG","@PG")
header[,2]=c(paste0("SN:",u),paste0("ID:",opt$sample.name),"ID:bwa")
header[1:length(u),3]=rep("LN:1",length(u))

############### TO FILTER THE PILEUPS SO THEY WILL NOT INCLUDE THE REFERENCE ALLELES
data$panel$target_seq[which(data$panel$probe_strand=="-")]=reverse(chartr("ATGC","TACG",data$panel$target_seq[which(data$panel$probe_strand=="-")])) #determine the reference alleles on the positive strand

tmp1=data.table("smMIP"=rep(data$panel$id,unlist(lapply(strsplit(data$panel$target_seq,""),function(x) length(x)))),
                "chr"=rep(data$panel$chr,unlist(lapply(strsplit(data$panel$target_seq,""),function(x) length(x)))),
                "pos"=unlist(apply(cbind(data$panel$target_start,data$panel$target_stop),1,function(x) seq(x[1],x[2],1))), #we are not interested to call mutations from the smMIPs arms
                "ref"=unlist(strsplit(data$panel$target_seq, "")))

################ SMMIP-LEVEL RAW PILEUP 
cat(paste("Creating smMIP-level raw pileups. Please notice, temporary files will be written in",opt$tmp.output))
cat("\n")

smmip_piles=split(data$samtable, by="smmip")
data$samtable=NULL

pile_raw=mclapply(1:length(smmip_piles), mc.cores = opt$threads, function(i) {
  pile=pileup_foreach_smmip(i)
  if(round(i/length(smmip_piles),2) %in% seq(0.01,0.99,0.01)) {
    system(paste0("printf '\\rCreating smMIP-level raw pileups :  ",round(100*i/length(smmip_piles)),"%%     '"))
  } else if (i==length(smmip_piles)){
    system(paste0("printf '\\rCreating smMIP-level raw pileups :  100%%     '"))
  }
  pile
})
system(paste0("printf '\\rCreating smMIP based pileups :  100%%     '"))
cat("\n")

pile_raw=rbindlist(pile_raw)
names(pile_raw)[1]="chr"
pile_raw=pile_raw[!(paste(smMIP,chr,pos,nucleotide) %in% paste(tmp1$smMIP,tmp1$chr,tmp1$pos,tmp1$ref))] #remove the reference alleles based on the nucleutides indicated in the smMIP design file

cat("Writing raw pileup file to disk\n")
write.table(pile_raw,file=paste0(opt$output,"/",opt$sample.name,"_raw_pileup.txt"),col.names = T,row.names = F,quote = F,sep = '\t')
rm(pile_raw)

################ SMMIP-LEVEL AND UMI-LEVEL PILEUP 
if(!((length(unique(data$panel$length.left.umi))==1 & 0 %in% unique(data$panel$length.left.umi)) &
     (length(unique(data$panel$length.right.umi))==1 & 0 %in% unique(data$panel$length.right.umi)))){
  
  cat(paste("Creating smMIP and UMI level pileups. Please notice, temporary files will be written in",opt$tmp.output))
  cat("\n")

  for(j in 1:length(smmip_piles)){
    smmip_umi_piles=split(rbindlist(smmip_piles[j]), by=c("smmip","umi"))
    idx=grep("flag",names(smmip_umi_piles))
    if(length(idx)>0){smmip_umi_piles=smmip_umi_piles[-grep("flag",names(smmip_umi_piles))]}
    if(length(smmip_umi_piles)==0){
      smmip_piles[[j]]=data.frame("seqnames"=NA, "pos"=NA, "strand"=NA, "nucleotide"=NA, "count"=NA, "coverage_at_position"=NA, "VAF"=NA, "smMIP"=NA, "UMI"=NA)
      smmip_piles[[j]]=smmip_piles[[j]][-1,]
    } else {
        smmip_piles[[j]]=rbindlist(mclapply(1:length(smmip_umi_piles), mc.cores = opt$threads, mc.cleanup=T,
               mc.silent=F, function(i) {
                 pile=pileup_foreach_smmip.umi(i)
                 pile
               }))
    }
    if(round(j/length(smmip_piles),2) %in% seq(0.01,0.99,0.01)) {
      system(paste0("printf '\\rCreating smMIP-UMI consensus pileups :  ",round(100*j/length(smmip_piles)),"%%     '"))
    } else if (j==length(smmip_piles)){
      system(paste0("printf '\\rCreating smMIP-UMI consensus pileups :  100%%     '"))
    }
  }
    
  system(paste0("printf '\\rCreating smMIP-UMI consensus pileups :  100%%     '"))
  cat("\n")
  
  pile_sscs=rbindlist(smmip_piles)
  rm(smmip_piles)
  pile_sscs[,c("a","b") := list(coverage_at_position,VAF)]
  pile_sscs[,coverage_at_position := sum(count),by=list(seqnames,pos,strand,smMIP)] #total SSCS coverage
  pile_sscs[,count := sum(count),by=list(seqnames,pos,strand,smMIP,nucleotide)] #total SSCS supporting reads for each allele
  if(opt$umi=="T"){
    pile_sscs[,c("family_sizes","VAF_in_families","UMI") := list(paste(a,collapse=","),paste(b,collapse=","),paste(UMI,collapse=",")),by=list(seqnames,pos,strand,smMIP,nucleotide)]
  } else if (opt$umi=="F"){
      pile_sscs[,c("family_sizes","VAF_in_families") := list(paste(a,collapse=","),paste(b,collapse=",")),by=list(seqnames,pos,strand,smMIP,nucleotide)]
    }
  pile_sscs=unique(pile_sscs[,c("VAF","a","b") := list(count/coverage_at_position,NULL,NULL)])
  names(pile_sscs)[1]="chr"
  pile_sscs=pile_sscs[!(paste(smMIP,chr,pos,nucleotide) %in% paste(tmp1$smMIP,tmp1$chr,tmp1$pos,tmp1$ref))] #remove the reference alleles based on the nucleutides indicated in the smMIP design file

  cat("Writing SSCS pileup file")
  write.table(pile_sscs,file=paste0(opt$output,"/",opt$sample.name,"_sscs_pileup.txt"),col.names = T,row.names = F,quote = F,sep = '\t')
}

if(del==1){
 unlink(opt$tmp.output, recursive=T)
}

cat("\n")
cat("###############################\n")
cat("             DONE\n")
cat("###############################\n")



