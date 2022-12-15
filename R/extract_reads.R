#!/usr/bin/env Rscript

################  LOAD OPTPARASE
library("optparse")

################ INPUT PARAMETERS
#define parameters
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,help="Path to the bam or sam file", metavar="character"),
  make_option(c("-p", "--position"), type="character", default=NA,help="Chromosom and genomic positions, e.g, chr13:28608217-28608354", metavar="character"),
  make_option(c("-i", "--identifier"), type="character", default=NA,help="Any character that can identify reads. e.g., accepted identifiers are chromosomes, full or partial CIGAR strings, full or partial read names and read sequences. Multiple features can be input in the following format: chr15,8M1D143M,AAGAGATGG. Comma corresponds to AND operator and vertical bar '|' to OR operator. Read help (Rscript extract_reads.R) to learn how to use this option properly", metavar="character"),
  make_option(c("-m", "--multiple"), type="character", default=NULL,help="Path to a file describing the location of multiple bam/sam files position and pattern. This should be a delimited text where the first line corresponds to file_path  position  pattern", metavar="character"),
  make_option(c("-s", "--split"), type="character", default="y",help="Options are 'y' or 'n'.One file per sample or not when -m is used", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,help="Path for output sam [MENDATORY]", metavar="character"))

#takes user input parameters
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser);

h=0
if(is.null(opt$output)){h=1}
if(is.null(opt$file) & is.null(opt$multiple)){h=1}
if(is.null(opt$position) & is.null(opt$identifier) & is.null(opt$multiple)){h=1}
if(!is.null(opt$file) & !is.null(opt$multiple)){h=1}

if(!is.null(opt$file) & (!is.null(opt$position) | !is.null(opt$identifier))) {o="option1"}
if(!is.null(opt$multiple)) {o="option2"}

if (h==1) {
  print_help(opt_parser)
  stop("Option1: -f, -o and at least one of -p or -i parameters must be supplied. Option2: use the parameters -m and -o.  
       This code can query a single position, or a single string pattern, or both a position and a pattern in a single sample.
       If you would like to extract reads from multiple samples provide your query as a tab delimited text file. Columns 1-3 shown be named 'file_path' 'position' 'pattern'.
       Write NA in case you want to search read based on position or pattern alone. Please notice that searching reads based on pattern alone might take long.  
       Example to run based on option1: Rscript extract_reads.R -f path_to_file -p chr13:28608217-28608354 -i '[0-9]D|[0-9]I' -o output_path
       Example to run based on option2: Rscript extract_reads.R -m path_to_file_with_multiple_queries -o output_path
       ", call.=FALSE)
}


################ LOAD LIBRARIES
library("Rsamtools")
library("data.table")

print(opt,quote = F) 

################ MAIN CODE
dts=list() #
if (o=="option1"){
  opt$multiple=data.table("file_path"=opt$file,"position"=opt$position, "pattern"=opt$identifier)
} else if (o=="option2"){ opt$multiple=fread(opt$multiple,header=T,sep="\t") }


for(i in 1:nrow(opt$multiple)){
   ################ Load the file
   #try to load a bam file with Rsamtools
   print("Loading the SAM/BAM file...",quote = F)
   param <- ScanBamParam(what="rname")
   a=try(as.data.table(scanBam(opt$multiple$file_path[i],param=param)[[1]]),silent = T)
   if(length(grep("Error",a[1]))>0){ #if error load sam
     dt=fread(opt$multiple$file_path[i],header=F,sep="\t")
     idx=grep("@", dt$V1[1:1000])
     if(length(idx)>0){dt=dt[-idx,]}
     colnames(dt)[1:11]=c("qname","flag","rname","pos","mapq","cigar","mrnm","mpos","isize","seq","qual")
   } else {# loading bam
     dt=as.data.table(do.call("DataFrame",scanBam(opt$multiple$file_path[i])[[1]]))
     dt=dt[,c("strand","qwidth") := NULL]
   }
   
  ################ Position based selection
  dts[[i]]=data.table()
   if(!is.na(opt$multiple$position[i])){ 
     dts[[i]]=dt[rname==strsplit(opt$multiple$position[i],":")[[1]][1] & pos>=as.numeric(strsplit(opt$multiple$position[i],":|-")[[1]][2]) & pos<as.numeric(strsplit(opt$multiple$position[i],":|-")[[1]][3])]
   } 
  ################ Pattren based selection
  if(!is.na(opt$multiple$pattern[i])){  
    if(nrow(dts[[i]])==0){dts[[i]]=dt}
    patterns=unlist(strsplit(opt$multiple$pattern[i],"[,]")[[1]])
    for(j in 1:length(patterns)){
      backup=dts[[i]]
      dts[[i]]=dts[[i]][qname==patterns[j] | rname==patterns[j] | pos==patterns[j] | cigar==patterns[j] | seq==patterns[j]]
      if(nrow(dts[[i]])==0){
        dts[[i]]=backup
        dts[[i]]=dts[[i]][apply(dts[[i]][,-1],1,function(x) grepl(patterns[j],paste(x,collapse = " ")))] #This could take long if the search was not set properly, especially with big files
      }
    }
  }
  dts[[i]]$qname=paste0(dts[[i]]$qname,"|",basename(opt$multiple$file_path[i]))
}	    

if(is.na(opt$multiple$position[i]) & is.na(opt$multiple$pattern[i])){ 
  print("No genomic locations or patterns to search were input.",quote = F)
  }

################ Merging all the reads and write the sam
if(opt$split=="n"){
 dt=unique(rbindlist(dts))
 dt=dt[-which(is.na(dt$mpos))]
 u=unique(dt$rname)
 header=as.data.frame(matrix(nrow = length(u)+2,ncol=11))
 header[1:(length(u)+2),]=""
 colnames(header)=colnames(dt)[1:11]
 header[,1]=c(rep("@SQ",length(u)),"@RG","@PG")
 header[,2]=c(paste0("SN:",u),"ID:extracted_reads","ID:bwa")
 header[1:length(u),3]=rep("LN:1",length(u))
 write.table(rbind(header,as.data.frame(dt[order(dt$rname,dt$pos),])),file=paste0(opt$output,"/","selected_reads.sam"),col.names = F,row.names = F,quote = F,sep = '\t') #Can be load to IGV
 print(paste("Wrote",dim(dt)[1],"reads to selected_reads.sam"),quote = F)
} else if(opt$split=="y"){
   for (i in 1:length(dts)){
     idx=which(is.na(dts[[i]]$mpos))
     if(length(idx)>0){dts[[i]]=dts[[i]][-idx]}   
     n=basename(opt$multiple$file_path[i])
     n=paste0(gsub(".bam","",n),"_",opt$multiple$position[i],"_",opt$multiple$pattern[i])
     n=gsub("_NA","",n)
     if(nrow(dts[[i]])!=0){
       u=unique(dts[[i]]$rname)
       header=as.data.frame(matrix(nrow = length(u)+2,ncol=11))
       header[1:(length(u)+2),]=""
       colnames(header)=colnames(dts[[i]])[1:11]
       header[,1]=c(rep("@SQ",length(u)),"@RG","@PG")
       header[,2]=c(paste0("SN:",u),"ID:extracted_reads","ID:bwa")
       header[1:length(u),3]=rep("LN:1",length(u))   
       write.table(rbind(header,as.data.frame(dts[[i]][order(dts[[i]]$rname,dts[[i]]$pos),])),file=paste0(opt$output,"/",n,".sam"),col.names = F,row.names = F,quote = F,sep = '\t')
       print(paste("Wrote",nrow(dts)[[i]], "reads to",paste0(opt$output,"/",n,".sam")))
     } else if(nrow(dts[[i]])==0){
       print(paste("The search resulted with no reads for",n)) 
     }
     }
  }
