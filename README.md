# smMIP-tools: a computational toolset for read processing and mutation detection from molecular-inversion-probes derived data

Here we present an example of how to use smMIP-tools main two features:
The first is a read processing tool that performs quality control steps, generates read-smMIP linkages and retrieval of molecular tags.
The second is an error-aware variant caller capable of detecting both single nucleotide variants and short insertion and deletions.

## Quick Start
smMIP-tools can be executed from the terminal. There is no need for installation. Copy the [code](https://www.) to your folder of choice.

## Dependencies
This pipeline requires the following software and packages:
| Program | Packages                                       |
| ------------------------|------------------------- |
| R (https://www.r-project.org) | optparse, data.table, parallel, Rsamtools, dplyr |               

   

## Read Processing
### Description
map_smMIPs_extract_UMIs.R takes as an input a paired-end read alignment bam file and a smMIP design file containing information about each probe and its targeted sequenced. It applys a set of filters on the input bam file to discard hard clipped reads, reads with low mapping quality, paired reads with an unexpected insert size or improper alignment orientations. To validate the proper structure of reads, and to identify corrupted UMI sequences read-smMIP linkages are being conducted. The final output contains a couple of quality control summary files concerning raw and consensus reads as well as a BAM file containing only high-quality reads. UMIs sequences and smMIPs-of-origin will be included in each read’s header.

### Configuration 
1) Make sure to have a [smMIP-design file](https://www.). Such files can easly be generated with [MIPgen](http://shendurelab.github.io/MIPGEN) executables.
2) Please generate a single folder and copy into it all the bam files that you want to analyse. We provide [bam files](http://) that can be used with this manual.  
map_smMIPs_extract_UMIs.R was built to process one bam file at a time. A simple shell script for parallel processing is provided [here](https:/www)

### Running the code
Run map_smMIPs_extract_UMIs.R with the required input parameters:
```
-b, --bam.file, Path to bam file. [MANDATORY]
-p, --panel.file, Path to smMIP design file. [MANDATORY]
-s, --sample.name, Sample ID that will be used for the output bam [MANDATORY]
-o, --output, Path for output files.  If not supplied, a new folder which is named based on the -s parameter will be generated within the folder that contain bam file.
-c, --code, Path to smMIPs_Function.R file. If not supplied, it assume the code share the same folder as this code folder with this code (map_smMIPs_extract_UMIs.R)
-f, --filtered.reads,  default="y", options="y" or "n". Output a sam file that contain the filtered reads. 
-t, --threads, default=1. Specify the number of threads to use. 
-O, --OVERLAP, default=0.95. Fine-tuning the overlap between reads and smMIPs . Used in the map.smip_to_site function
-M, --MAPQ, default=50. MAPQ cut-off. Used in the filter.on.mappingscore function
  ```





Who to contact: sagi.abelson@gmail.com
