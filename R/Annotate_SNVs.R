#!/usr/bin/env Rscript

library("optparse")

################ INPUT PARAMETERS
#define parameters
option_list = list(
  make_option(c("-p", "--panel.file"), type="character", default=NULL,help="Path to smMIP design file [MENDATORY if -i was not supplied]", metavar="character"),
  make_option(c("-i", "--input"), type="character", default=NULL,help="Path to a tab delimited table with the following columns: chr, pos, ref, alt. [MENDATORY if -p was not supplied]", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,help="Path for the output annotated panel. If not supplied, the output will be saved within the folder that contain the smMIP design file/input file", metavar="character"),
  make_option(c("-c", "--code"), type="character", default=getwd(),help="Path to smMIP tools source functions, smMIPs_Function.R file. If not supplied, it assume the code share the same folder as this code folder with this code (Annotate_SNVs.R)", metavar="character"),
  make_option(c("-t", "--threads"), type="integer", default=1,help="Specify the number of threads to use for parallel processing", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default="GRCh37",help="Specify the genome assembly to be queried. For all the possible options please refer to the cellbaseR package documentation", metavar="character"),
  make_option(c("-s", "--species"), type="character", default="hsapiens",help="Specify the species to be queried. For all the possible options please refer to the cellbaseR package documentation", metavar="character"))

#takes user input parameters
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser);
if ( is.null(opt$panel.file) & is.null(opt$input) ) {
  print_help(opt_parser)
  stop("Either the -p or -i parameter must be supplied", call.=FALSE)
}

if ( !is.null(opt$panel.file) & !is.null(opt$input) ) {
  print_help(opt_parser)
  stop("Either the -p or -i parameter must be supplied. Not both", call.=FALSE)
}

if ( is.null(opt$panel.file) & !is.null(opt$input)) {
  opt$output=dirname(opt$input)
  file.name=basename(opt$input)
}

if ( !is.null(opt$panel.file) & is.null(opt$input)) {
  opt$output=dirname(opt$panel.file)
  file.name=basename(opt$panel.file)
}


if (is.null(opt$code)){ #load source functions
  source(paste0(getwd,"/smMIPs_Functions.R"))
} else {
  source(paste0(opt$code,"/smMIPs_Functions.R"))
}


################  LOAD LIBRARIES
library("data.table")
library("cellbaseR")
library("IRanges")


#Print all the parameters
cat("###############################\n")
cat("        Run Parameters\n")
cat("###############################\n")
print(opt)
cat("###############################\n")
cat("            Running...\n")
cat("###############################\n")


################  CREATING A TABLE WITH ALL THE POSSIBLE OBSERVED ALLELES USING THE TARGETED GENOMIC LOCI. ALTERNATIVELY, LOADING THE USER INPUT FILE
if(!is.null(opt$panel.file)){
  panel<-load.panel(opt$panel.file)
  panel$target_seq[which(panel$probe_strand=="-")]=reverse(chartr("ATGC","TACG",panel$target_seq[which(panel$probe_strand=="-")])) #determine the reference alleles on the positive strand

  tmp1=data.table("smMIP"=rep(panel$id,unlist(lapply(strsplit(panel$target_seq,""),function(x) length(x)))),
                  "chr"=rep(panel$chr,unlist(lapply(strsplit(panel$target_seq,""),function(x) length(x)))),
                  "pos"=unlist(apply(cbind(panel$target_start,panel$target_stop),1,function(x) seq(x[1],x[2],1))), #we are not interested to call mutations from the smMIPs arms
                  "ref"=unlist(strsplit(panel$target_seq, "")),
                  "alt"="A")
  tmp2=tmp1
  tmp2$alt="T"
  tmp3=tmp1
  tmp3$alt="C"
  tmp4=tmp1
  tmp4$alt="G"

  genomicBases=rbind(tmp1,tmp2,tmp3,tmp4)
  colnames(genomicBases)[3]="pos"
  rm(tmp1,tmp2,tmp3,tmp4)
  idx=1:nrow(genomicBases)
} else if (!is.null(opt$input)){
  genomicBases=fread(opt$input,header=T,sep="\t")
  if(length(grep("annotated",names(genomicBases)))>0){
    idx=which(genomicBases$annotated=="FALSE")
  } else {
    idx=1:nrow(genomicBases)
  }
}

################  POPULATE THE TABLE
cb <- CellBaseR(species = opt$species)
cbparam <- CellBaseParam(assembly = opt$genome) #define genome

print("Annotating. Please wait...",quote = F)
a=mclapply(idx, mc.cores = opt$threads, mc.cleanup=T, mc.silent=F ,function(i) {
  variant_annotation_using_cellbaseR(i)
})

x=sapply(rbindlist(a),class)
idx2=which(!(names(x) %in% names(genomicBases)))
if(length(idx2)>0){
  for(j in names(x)[idx2]){
    genomicBases[,j]=""
  }
}
genomicBases[] <- lapply(genomicBases, as.character)
genomicBases[idx,c("gene","protein","cosmic","maf","variant_type","cadd_scaled","annotated")]=rbindlist(a)


###############  RENAME AMINO ACIDS WITH THEIR SINGLE LETTER ABBREVIATIONS
f=c("STOP","ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR")
r=c("*","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")

for(i in 1:length(f)){
  genomicBases$protein=gsub(f[i],r[i],genomicBases$protein)
}

system(paste0("printf '\\rDownloading SNV annotations :  100%%     '"))
cat("\n")
cat("Writing the annotated table to disk\n")

u=unique(genomicBases$annotated)
if(length(u)==1 & "TRUE" %in% u){genomicBases=genomicBases[,-"annotated"]}

write.table(genomicBases,file=paste0(opt$output,"/annotated_",file.name),col.names = T,row.names = F,quote = F,sep = '\t')

cat("\n")
cat("###############################\n")
cat("             DONE\n")
cat("###############################\n")
