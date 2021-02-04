#!/bin/bash -l
#$ -cwd
#$ -S /bin/bash

############################################       R E A D  M E      ######################################################
## Built on UGE 8.6.4. The code may need editing if you are working with a different grid engine
## This code will run map_smMIPs_extract_UMIs.R for every bam file in the input folder.
## Only the mandatory parameters are used. If you want to change the defaults, you must add your changes to the code below.
## This code file needs to be in the same folder as map_smMIPs_extract_UMIs.R to work properly.
###########################################################################################################################

#Enter the path to the directory containing the bam files
BAM=/.mounts/labs/abelsonlab/private/smMIP_ARCH/CODE/example_github/Example/bams/  #CHANGE HERE

#Enter the path to the smMIP-design file
TARGET=/.mounts/labs/abelsonlab/private/smMIP_ARCH/CODE/example_github/Example/supplemental_files/Target_MIPgen.txt  #CHANGE HERE


for file in $( ls $BAM | grep bam); do
	
	### GENERATING THE SCRIPT 
	SAMPLE="$(echo $file | sed 's/.bam.*//')" #Please note the sample name and the file name assumed to be the same (everything after .bam will be removed). 
	echo -e "#!/bin/bash -l\n#$ -S /bin/bash\n#$ -cwd\n\n" > automatic_script_for_map_smMIPs_extract_UMIs_$SAMPLE.sh
	echo -e "module load rstats/3.6\n" >> automatic_script_for_map_smMIPs_extract_UMIs_$SAMPLE.sh #Loading R. Delete this line if R is rooted. Change the module name to your R module if it isnt.
	echo -e "Rscript map_smMIPs_extract_UMIs.R -b $BAM/$file -p $TARGET -s $SAMPLE" >> automatic_script_for_map_smMIPs_extract_UMIs_$SAMPLE.sh 
	
	### RUNNING THE SCRIPT
	qsub -P abelsonlab -l h_vmem=16g automatic_script_for_map_smMIPs_extract_UMIs_$SAMPLE.sh	

done

