#!/usr/bin/env Rscript

library("optparse")

################ INPUT PARAMETERS
#define parameters
option_list = list(
  make_option(c("-i", "--input.file"), type="character", default=NULL,help="Path to the raw annovar annotated panel file [MENDATORY]", metavar="character"))

#takes user input parameters
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser);
if ( is.null(opt$input.file)) {
  print_help(opt_parser)
  stop("The -i parameter must be supplied", call.=FALSE)
}

#/u/sabelson/ANNOVAR/annovar/table_annovar.pl input_for_annovar_Myeloid_Panel_Targets.txt /.mounts/labs/abelsonlab/private/annovar/datasets/humandb -buildver hg19 -out test_output -protocol refGene,cosmic96_coding,gnomad_genome,cadd13 -operation g,f,f,f -nastring NA --remove --otherinfo


################  LOAD LIBRARIES
library("data.table")

#Print all the parameters
cat("###############################\n")
cat("        Run Parameters\n")
cat("###############################\n")
print(opt)

input=fread(opt$input.file,header=T,sep="\t",showProgress = FALSE)
V1=input$ExonicFunc.refGene
V1[which(is.na(V1))]=input$Func.refGene[which(is.na(V1))]
input$V1=V1
input=input[,c("Otherinfo1","Chr","Start","Ref","Alt","Gene.refGene","AAChange.refGene","cosmic96_coding","gnomAD_genome_ALL","V1","CADD13_PHRED")]

colnames(input)=c("smMIP","chr","pos","ref","alt","gene","protein","cosmic","maf","variant_type","cadd_scaled")
if(length(grep("chr",head(input$chr)))==0){input$chr=paste0("chr",input$chr)}

write.table(input,file=paste0(dirname(opt$input),"/annotated_Target_MIPgen.txt"),col.names = T,row.names = F,quote = F,sep = '\t')

cat("\n")
cat("###############################\n")
cat("             DONE\n")
cat("###############################\n")

print(paste0("Wrote the final annotation file for smMIP-tools in"," ",dirname(opt$input)))



