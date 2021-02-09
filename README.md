# smMIP-tools: a computational toolset for processing and analysis of single molecule molecular-inversion-probes derived data

Here we present an example of how to use smMIP-tools to process and analyze sequencing data obtained from targeted smMIP panels. 
The pipeline covers the following topics:
1) Read processing and quality control steps, construction of read-smMIP linkages and retrieval of molecular tags.
2) smMIP-level base calls summary.
3) smMIP panel annotation
4) Mutation calling using a error-aware variant caller capable of detecting both single nucleotide variants and short insertion and deletions.

## Quick Start
The supplied example data contain sequencing information obtained by a single smMIP that covers a recurrent mutation of which we will attemp to identify.
smMIP-tools code can be executed from the terminal. There is no need for installation. Copy the [code files](https://github.com/BioSoft/smMIP-tools/tree/main/R) to your folder of choice.

## Dependencies
This pipeline requires the following software and packages:
| Program | Packages                                       |
| ------------------------|------------------------- |
| R (https://www.r-project.org) | optparse, data.table, parallel, Rsamtools, dplyr, IRanges, cellbaseR |               

*You will be prompt to install/update other packages when installing the above packages.
   

## Read Processing
### Description
map_smMIPs_extract_UMIs.R takes as an input a paired-end read alignment bam file and a smMIP design file containing information about each probe and its targeted sequenced. It applys a set of filters on the input bam file to discard hard clipped reads, reads with low mapping quality, paired reads with an unexpected insert size or improper alignment orientations. To validate the proper structure of reads, and to identify corrupted UMI sequences read-smMIP linkages are being conducted. The final output contains a couple of quality control summary files concerning raw and consensus reads as well as a BAM file containing only high-quality reads. UMIs sequences and smMIPs-of-origin will be included in each read’s header.

### Configuration 
1) Make sure to have a [smMIP-design file](https://github.com/BioSoft/smMIP-tools/blob/main/Example/supplemental_files/Target_MIPgen.txt). We used [MIPgen](http://shendurelab.github.io/MIPGEN) to design the smMIP used in this example. If you designed your smMIPs differently make sure that your file contain the following columns (case sensitve) and their relavant information: chr, ext_probe_start, ext_probe_stop, ext_probe_sequence, lig_probe_start, lig_probe_stop, lig_probe_sequence, mip_scan_start_position, mip_scan_stop_position, scan_target_sequence, mip_sequence ('N' represent UMIs), probe_strand, mip_name
2) Please generate a single folder and copy into it all the bam files that you want to analyse. We provide [bam files](https://github.com/BioSoft/smMIP-tools/tree/main/Example/bams) that can be used with this manual.  

### Running the code
Run map_smMIPs_extract_UMIs.R with the required input parameters:
```
-b, --bam.file. Path to bam file. [MANDATORY]
-p, --panel.file. Path to smMIP design file. [MANDATORY]
-s, --sample.name. Sample ID that will be used to name the output bam. [MANDATORY]
-o, --output, Path for output files.  Path for output files. A new folder which is named based on the -s parameter will be generated by default. If not supplied, the sample folder will be generated within the folder that contain the bam files.
-c, --code, Path to smMIP tools source functions, smMIPs_Function.R file. If not supplied, it assume the code share the same folder as this code folder with this code (map_smMIPs_extract_UMIs.R).
-f, --filtered.reads,  default="n", options="y" or "n". Output a sam file that contain the filtered reads. A sam file for the non-filtered reads will also be generated.  
-t, --threads, default=1. Specify the number of threads to use for parallel processing. 
-O, --OVERLAP, default=0.95. Fine-tuning the overlap between reads and smMIPs. Used in the map.smip_to_site function
-M, --MAPQ, default=50. MAPQ cut-off. Used in the filter.on.mappingscore function
  ```
\*map_smMIPs_extract_UMIs.R was built to process one bam file at a time. An example shell script to assign jobs to multiple HPC cluster cores is provided [here](https://github.com/BioSoft/smMIP-tools/blob/main/R/map_smMIPs_extract_UMIs_MULTIPLE_BAMS.sh)

```
Rscript /.mounts/example_github/R/map_smMIPs_extract_UMIs.R -b /.mounts/example_github/Example/bams/control1.bam -p /.mounts/example_github/Example/supplemental_files/Target_MIPgen.txt -s control1
  
Loading required package: GenomeInfoDb
Loading required package: BiocGenerics
Loading required package: parallel

Attaching package: ‘BiocGenerics’

The following objects are masked from ‘package:parallel’:

    clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    clusterExport, clusterMap, parApply, parCapply, parLapply,
    parLapplyLB, parRapply, parSapply, parSapplyLB

The following objects are masked from ‘package:stats’:

    IQR, mad, sd, var, xtabs

The following objects are masked from ‘package:base’:

    anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    union, unique, unsplit, which, which.max, which.min

Loading required package: S4Vectors
Loading required package: stats4

Attaching package: ‘S4Vectors’

The following object is masked from ‘package:base’:

    expand.grid

Loading required package: IRanges
Loading required package: GenomicRanges
Loading required package: Biostrings
Loading required package: XVector

Attaching package: ‘Biostrings’

The following object is masked from ‘package:base’:

    strsplit


Attaching package: ‘dplyr’

The following objects are masked from ‘package:Biostrings’:

    collapse, intersect, setdiff, setequal, union

The following object is masked from ‘package:XVector’:

    slice

The following objects are masked from ‘package:GenomicRanges’:

    intersect, setdiff, union

The following object is masked from ‘package:GenomeInfoDb’:

    intersect

The following objects are masked from ‘package:IRanges’:

    collapse, desc, intersect, setdiff, slice, union

The following objects are masked from ‘package:S4Vectors’:

    first, intersect, rename, setdiff, setequal, union

The following objects are masked from ‘package:BiocGenerics’:

    combine, intersect, setdiff, union

The following objects are masked from ‘package:stats’:

    filter, lag

The following objects are masked from ‘package:base’:

    intersect, setdiff, setequal, union


Attaching package: ‘data.table’

The following objects are masked from ‘package:dplyr’:

    between, first, last

The following object is masked from ‘package:GenomicRanges’:

    shift

The following object is masked from ‘package:IRanges’:

    shift

The following objects are masked from ‘package:S4Vectors’:

    first, second

###############################
        Run Parameters
###############################
$bam.file
[1] "/.mounts/example_github/Example/bams/control1.bam"

$panel.file
[1] "/.mounts/example_github/Example/supplemental_files/Target_MIPgen.txt"

$sample.name
[1] "control1"

$code
[1] "/.mounts/example_github/R"

$filtered.reads
[1] "n"

$threads
[1] 1

$OVERLAP
[1] 0.95

$MAPQ
[1] 50

$help
[1] FALSE

$output
[1] "/.mounts/example_github/Example/bams"

###############################
            Running...
###############################
Loading bam
Filtering reads based on low mapping score
Filtering reads based on bad sam flag
Filtering reads based on mate distance
Filtering hard clipped reads
Filtering reads with missing mates
Mapping smMIPs to reads. Considering overlap and distance measurements :  100%      
Verifying mapping by local smMIP arms alignment : 100%%     
Extracting UMIs
Writing coverage per smMIP and filtered reads summary files
Adding smMIPs names and UMI sequences to the reads names
Writing UMIs per smmip summary
Writing the clean bam file
[1] "/.mounts/example_github/Example/bams/control1/control1_clean.bam"

###############################
             DONE
###############################
  ```

## smMIP Level Base Calls Summary
### Description
smMIP_level_raw_and_consensus_pileups.R takes as an input a clean bam file generated by map_smMIPs_extract_UMIs.R and a smMIP design file to generate a smMIP-level base calls summaries. smMIP-level base calls summaries provides an additional way to identify errors and differentiate them from real mutations. If your smMIPs include unique molecular identifers, consensus base calls summaries will also be generated. If you are not familiar with the concept of consensus base calls you may read our Nucleic Acid Research paper ([PMID: 31127310](https://academic.oup.com/nar/article/47/15/e87/5498633)).

### Configuration 
Make sure to have a proper smMIP design file. It is being used to infer automatically the length of UMIs, remove reference alleles' summaries to reduce file size and to focus the query only on the targeted loci of each smMIP. The porpose of the latter is to withhold analyis of genomic loci associated with smMIP's arms hybridization. 

### Running the code
Run smMIP_level_raw_and_consensus_pileups.R with the required input parameters:

```
-b, --bam.file. Path to the filtered bam file output by map_smMIPs_extract_UMIs.R (sample_clean.bam) [MENDATORY]
-p, --panel.file. Path to smMIP design file [MENDATORY]
-s, --sample.name. Sample ID that will be used to name the output file(s). If not provided the name of the folder containing the bam file is assumed to be the sample name
-o, --output. Path for the output pileup file(s). If not supplied, the output will be saved within the folder that contain the bam file
-c, --code. Path to smMIP tools source functions, smMIPs_Function.R file. If not supplied, it assume the code share the same folder as this code folder with this code (smMIP_specific_raw_and_consensus_pileups.R)"
-d, --mnd, default=1. Minimum depth to consider in the pileup
-m, --mmq, default=50. Minimum mapping quality to consider in the pileup
-q, --mbq"),default=10. Minimum base quality to consider in the pileup
-r, --rank, default="F". Options are 'F' or 'T'. Reporting all the possible alleles (A,C,T,G,-,+) if they observed. When 'T', only the allele with the most read support (other than the reference) at each genomic postion will be reported. This option can help to remove many errors from the data yet it comes with the risk of missing real mutations.
-u, --umi, default="T". Options are 'T' or 'F'. Reporting the UMI sequences that are associated with each allele will be reported. Provides with the data needed to estimate read-to-sample misassignment
-f, --family.size, default=0. The minimum read-family size to consider in single strand consensus pileups
-v, --consensus.cutoff", default=0.7. The percentage of reads that must have the same base to form a consensus
-t, --threads, default=1. Specify the number of threads to use"
  ```

\*smMIP_level_raw_and_consensus_pileups.R was built to process one bam file at a time. An example shell script to assign jobs to multiple HPC cluster cores is provided [here](https://github.com/BioSoft/smMIP-tools/blob/main/R/smMIP_level_raw_and_consensus_pileups_MULTIPLE_BAMS.sh)

```
Rscript smMIP_level_raw_and_consensus_pileups.R -b /.mounts/example_github/Example/bams/control1/control1_clean.bam -p /.mounts/example_github/Example/supplemental_files/Target_MIPgen.txt
  
Loading required package: GenomeInfoDb
Loading required package: BiocGenerics
Loading required package: parallel

Attaching package: ‘BiocGenerics’

The following objects are masked from ‘package:parallel’:

    clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    clusterExport, clusterMap, parApply, parCapply, parLapply,
    parLapplyLB, parRapply, parSapply, parSapplyLB

The following objects are masked from ‘package:stats’:

    IQR, mad, sd, var, xtabs

The following objects are masked from ‘package:base’:

    anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    union, unique, unsplit, which, which.max, which.min

Loading required package: S4Vectors
Loading required package: stats4

Attaching package: ‘S4Vectors’

The following object is masked from ‘package:base’:

    expand.grid

Loading required package: IRanges
Loading required package: GenomicRanges
Loading required package: Biostrings
Loading required package: XVector

Attaching package: ‘Biostrings’

The following object is masked from ‘package:base’:

    strsplit


Attaching package: ‘data.table’

The following object is masked from ‘package:GenomicRanges’:

    shift

The following object is masked from ‘package:IRanges’:

    shift

The following objects are masked from ‘package:S4Vectors’:

    first, second

###############################
        Run Parameters
###############################
$bam.file
[1] "/.mounts/example_github/Example/bams/control1/control1_clean.bam"

$panel.file
[1] "/.mounts/example_github/Example/supplemental_files/Target_MIPgen.txt"

$code
[1] "/.mounts/example_github/R"

$mnd
[1] 1

$mmq
[1] 50

$mbq
[1] 10

$rank
[1] "F"

$umi
[1] "T"

$family.size
[1] 0

$consensus.cutoff
[1] 0.7

$threads
[1] 1

$help
[1] FALSE

$output
[1] "/.mounts/example_github/Example/bams/control1"

$sample.name
[1] "control1"

###############################
            Running...
###############################
Loading bam
Creating smMIP-level raw pileups. Please notice, temporary files will be written in /.mounts/example_github/Example/bams/control1
Creating smMIP based pileups :  100%         
Writing raw pileup file to disk
Creating smMIP and UMI level pileups. Please notice, temporary files will be written in /.mounts/example_github/Example/bams/control1
Creating smMIP-UMI consensus pileups :  100%     
Writing SSCS pileup file
###############################
             DONE
###############################
  ```

The output:

```
head /.mounts/example_github/Example/bams/control1/control1_raw_pileup.txt
  
chr	pos	strand	nucleotide	count	coverage_at_position	VAF	smMIP
chr9	5073732	+	A	9	7001	0.00128553063848022	JAK2_001
chr9	5073732	+	G	2	7001	0.000285673475217826	JAK2_001
chr9	5073733	+	A	10	7001	0.00142836737608913	JAK2_001
chr9	5073733	+	C	1	7001	0.000142836737608913	JAK2_001
chr9	5073733	+	G	3	7001	0.000428510212826739	JAK2_001
chr9	5073734	+	A	12	7001	0.00171404085130696	JAK2_001
chr9	5073734	+	C	1	7001	0.000142836737608913	JAK2_001
chr9	5073735	+	A	4	7001	0.000571346950435652	JAK2_001
chr9	5073735	-	A	20	6936	0.00288350634371396	JAK2_001
  ```
```
head /.mounts/example_github/Example/bams/control1/control1_sscs_pileup.txt
  
chr	pos	strand	nucleotide	count	coverage_at_position	VAF	smMIP	UMI	family_sizes	VAF_in_families
chr9	5073844	+	G	1	318	0.00314465408805031	JAK2_001	TTAA_AATC	1	1
chr9	5073844	-	G	1	320	0.003125	JAK2_001	TTAA_AATC	1	1
chr9	5073775	+	G	1	320	0.003125	JAK2_001	TCAT_AAGT	1	1
chr9	5073789	+	C	1	320	0.003125	JAK2_001	TCAT_AAGT	1	1
chr9	5073845	+	T	10	316	0.0316455696202532	JAK2_001	TTGT_CCTC,CAAT_ATTC,TTAT_AACA,AATT_TTAT,TTAT_TAAT,TATC_TTAT,TTTT_TCAA,AATG_TGTT,TTTT_TACT,CTAT_ACAT	1,1,1,1,1,1,1,1,1,1	1,1,1,1,1,1,1,1,1,1
chr9	5073741	-	C	1	317	0.00315457413249211	JAK2_001	TTTA_GTAC	1	1
chr9	5073753	-	A	2	320	0.00625	JAK2_001	TTTA_GTAC,TCAA_ATTT	1,1	1,1
chr9	5073759	-	C	2	319	0.00626959247648903	JAK2_001	TTTA_GTAC,TCAA_TATT	1,1	1,1
chr9	5073828	+	G	2	318	0.00628930817610063	JAK2_001	TTTA_GTAC,GTTT_AATA	1,1	1,1
  ```

## smMIP Target Panel Annotation
### Description
Annotate_SNVs.R can take two types of files as an input: 1) A [smMIP-design file](https://github.com/BioSoft/smMIP-tools/blob/main/Example/supplemental_files/Target_MIPgen.txt) or 2) a tab delimited file with the columns "chr", "pos", "ref", "alt" . The code uses the  "cellbaseR" package to summarize gene ID, amino acid changes, Cosmic, Minor allele Frequency (MAF), Variant type (synonymous, missense, etc..) and deleteriousness (CADD score) annotations. If the input is a smMIP-design file, all the possible SNVs that might be detected using each smMIP in your target smMIP panel will be annotated. 

### Configuration 
cellbaseR makes use of a web service API. Make sure to run it on a processor(s) that is connected to the internet. If the connection was interupted the code will continue to run until it reaches the last variant. A column named "annotated" will be generated to mark variant that were successfully annotated and those that were not. If this column exists in your output you can try to run Annotate_SNVs.R once more using this output and the -i parameter. If the "annotated" column does not exits it means that all the variants were successfully annotated.

### Running the code
Run Annotate_SNVs.R with the required input parameters:

```
-p, --panel.file. "Path to smMIP design file [MENDATORY if -i was not supplied]
-i, --input. Path to a tab delimited table with the following columns: chr, pos, ref, alt. [MENDATORY if -p was not supplied]", metavar="character"),
-o, --output. Path for the output annotated panel. If not supplied, the output will be saved within the folder that contain the smMIP design file/input file.
-c, --code. Path to smMIP tools source functions, smMIPs_Function.R file. If not supplied, it assume the code share the same folder as this code folder with this code (Annotate_SNVs.R)
-t, --threads, default=1. Specify the number of threads to use for parallel processing. 
-g, --genome, default="GRCh37". Specify the genome assembly to be queried. For all the possible options please refer to the cellbaseR package documentation. 
-s, --species, default="hsapiens". Specify the species to be queried. For all the possible options please refer to the cellbaseR package documentation.
 ```
    
```
Rscript Annotate_SNVs.R -p /.mounts/example_github/Example/supplemental_files/Target_MIPgen.txt
  
Loading required package: BiocGenerics
Loading required package: parallel

Attaching package: ‘BiocGenerics’

The following objects are masked from ‘package:parallel’:

    clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    clusterExport, clusterMap, parApply, parCapply, parLapply,
    parLapplyLB, parRapply, parSapply, parSapplyLB

The following objects are masked from ‘package:stats’:

    IQR, mad, sd, var, xtabs

The following objects are masked from ‘package:base’:

    anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    union, unique, unsplit, which, which.max, which.min

Loading required package: S4Vectors
Loading required package: stats4

Attaching package: ‘S4Vectors’

The following objects are masked from ‘package:data.table’:

    first, second

The following object is masked from ‘package:base’:

    expand.grid


Attaching package: ‘IRanges’

The following object is masked from ‘package:data.table’:

    shift

###############################
        Run Parameters
###############################
$panel.file
[1] "/.mounts/example_github/Example/supplemental_files/Target_MIPgen.txt"

$code
[1] "/.mounts/example_github/R"

$threads
[1] 1

$genome
[1] "GRCh37"

$species
[1] "hsapiens"

$help
[1] FALSE

$output
[1] "/.mounts/example_github/Example/supplemental_files"

###############################
            Running...
###############################
[1] Annotating. Please wait...
Downloading SNV annotations :  100%     
Writing the annotated table to disk

###############################
             DONE
###############################
  ```
  
The output:

```
head /.mounts/example_github/Example/supplemental_files/annotated_Target_MIPgen.txt
  
smMIP	chr	pos	ref	alt	gene	protein	cosmic	maf	variant_type	cadd_scaled
JAK2_001	chr9	5073732	T	A	JAK2	ENST00000539801:L604H; ENST00000381652:L604H; ENST00000544510:L455H	NA	NA	missense_variant	31
JAK2_001	chr9	5073733	T	A	JAK2	ENST00000539801:L604L; ENST00000381652:L604L; ENST00000544510:L455L	NA	NA	synonymous_variant	9.98999977111816
JAK2_001	chr9	5073734	T	A	JAK2	ENST00000539801:S605T; ENST00000381652:S605T; ENST00000544510:S456T	NA	NA	missense_variant	23.7999992370605
JAK2_001	chr9	5073735	C	A	JAK2	ENST00000539801:S605Y; ENST00000381652:S605Y; ENST00000544510:S456Y	NA	NA	missense_variant	33
JAK2_001	chr9	5073736	T	A	JAK2	ENST00000539801:S605S; ENST00000381652:S605S; ENST00000544510:S456S	NA	NA	synonymous_variant	8.86999988555908
JAK2_001	chr9	5073737	C	A	JAK2	ENST00000539801:H606N; ENST00000381652:H606N; ENST00000544510:H457N	NA	NA	missense_variant	27.2999992370605
JAK2_001	chr9	5073738	A	A	JAK2	ENST00000539801:H606H; ENST00000381652:H606H; ENST00000544510:H457H	NA	NA	synonymous_variant	0
JAK2_001	chr9	5073739	C	A	JAK2	ENST00000539801:H606Q; ENST00000381652:H606Q; ENST00000544510:H457Q	COSM29116,haematopoietic and lymphoid tissue	0.000144591296666667	missense_variant	26.8999996185303
JAK2_001	chr9	5073740	A	A	JAK2	ENST00000539801:K607K; ENST00000381652:K607K; ENST00000544510:K458K	NA	NA	synonymous_variant	0
  ```




Who to contact: sagi.abelson@gmail.com
