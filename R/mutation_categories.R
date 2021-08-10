#!/usr/bin/env Rscript

library("optparse")

################ INPUT PARAMETERS
#define parameters
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,help="Path to the output from calling_mutations.R. i.e., called_mutations.txt [MENDATORY]", metavar="character"),
  make_option(c("-c", "--code"), type="character", default=getwd(),help="Path to smMIP tools source functions file. If not supplied, it is assumed that the file (smMIPs_Function.R) is located in your working directory", metavar="character"))

#takes user input parameters
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser);
if ( is.null(opt$panel.file) & is.null(opt$input) ) {
  print_help(opt_parser)
  stop("The -i parameter must be supplied", call.=FALSE)
}

if (is.null(opt$code)){ #load source functions
  source(paste0(getwd,"/smMIPs_Functions.R"))
} else {
  source(paste0(opt$code,"/smMIPs_Functions.R"))
}

################  LOAD LIBRARIES
library("data.table")

calls=fread(opt$input,header=T,sep="\t",showProgress = FALSE)

categorize()
