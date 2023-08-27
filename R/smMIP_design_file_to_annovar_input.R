#!/usr/bin/env Rscript

library("optparse")

################ INPUT PARAMETERS
#define parameters
option_list = list(
  make_option(c("-p", "--panel.file"), type="character", default=NULL,help="Path to smMIP design file [MENDATORY if -i was not supplied]", metavar="character"),
  make_option(c("-c", "--code"), type="character", default=getwd(),help="Path to smMIP tools source functions, smMIPs_Function.R file. If not supplied, it assumes that you are executing this code (smMIP_design_file_to_annovar_input.R) from the same folder where the smMIPs_Functions.R file is" , metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,help="Path for the output new formated file If not supplied, the output will be saved within the folder that contain the smMIP design file/input file", metavar="character"))

#takes user input parameters
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser);
if ( is.null(opt$panel.file)) {
  print_help(opt_parser)
  stop("The -p parameter must be supplied", call.=FALSE)
}

if (is.null(opt$output)) {
  opt$output=dirname(opt$panel.file)
}

if (is.null(opt$code)){ #load source functions
  source(paste0(getwd,"/smMIPs_Functions.R"))
} else {
  source(paste0(opt$code,"/smMIPs_Functions.R"))
}

file.name=basename(opt$panel.file)
################  LOAD LIBRARIES
library("data.table")
library("IRanges")


#Print all the parameters
cat("###############################\n")
cat("        Run Parameters\n")
cat("###############################\n")
print(opt)

################  CREATING A TABLE WITH ALL THE POSSIBLE OBSERVED ALLELES USING THE TARGETED GENOMIC LOCI. ALTERNATIVELY, LOADING THE USER INPUT FILE

panel<-load.panel(opt$panel.file)
panel$target_seq[which(panel$probe_strand=="-")]=reverse(chartr("ATGC","TACG",panel$target_seq[which(panel$probe_strand=="-")])) #determine the reference alleles on the positive strand

tmp1=data.table("chr"=rep(panel$chr,unlist(lapply(strsplit(panel$target_seq,""),function(x) length(x)))),
                  "pos1"=unlist(as.data.frame(apply(cbind(panel$target_start,panel$target_stop),1,function(x) seq(x[1],x[2],1)))), #we are not interested to call mutations from the smMIPs arms
                  "pos2"=unlist(as.data.frame(apply(cbind(panel$target_start,panel$target_stop),1,function(x) seq(x[1],x[2],1)))),
                  "ref"=unlist(strsplit(panel$target_seq, "")),
                  "alt"="A",
                  "smMIP"=rep(panel$id,unlist(lapply(strsplit(panel$target_seq,""),function(x) length(x)))))
tmp2=tmp1
tmp2$alt="T"
tmp3=tmp1
tmp3$alt="C"
tmp4=tmp1
tmp4$alt="G"

genomicBases=rbind(tmp1,tmp2,tmp3,tmp4)
rm(tmp1,tmp2,tmp3,tmp4)

genomicBases=genomicBases[genomicBases$alt!=genomicBases$ref,]
if(length(grep("chr",head(genomicBases$chr))>0)){genomicBases$chr=gsub("chr","",genomicBases$chr)}
write.table(genomicBases,file=paste0(opt$output,"/input_for_annovar_",file.name),col.names = F,row.names = F,quote = F,sep = '\t')

cat("\n")
cat("###############################\n")
cat("             DONE\n")
cat("###############################\n")

print(paste0("Wrote the input_for_annovar file in"," ",opt$output))

