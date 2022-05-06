# Pseudo-nitzschia2bRAD
Pseudo-nitzschia spp. produces a neurotoxin, domoic acid (DA), which can cause illness in humans, commonly known as amnesic shellfish poisoning, mass-mortality of marine animals, and closure of commercial shellfish fisheries during bloom events. Understanding and forecasting these harmful algal blooms is a primary management goal; however, accurately predicting the onset and severity of Pseudo-nitzschia bloom events remains difficult, in part because the underlying drivers of bloom formation have not been fully resolved. Recent work suggests that the genetic composition of a Pseudo-nitzschia bloom may be a better predictor of toxicity than prevailing environmental conditions. Here, we test the ability of reduced representation sequencing (2b-RAD) to delineate inter and intra-species level genetic variation of Pseudo-nitzschia. We generated 2bRAD libraries from pure cultures of Pseudo-nitzschia spp. and mock community mixes. We also generated sequencing libraries from longitudinal field samples from the San Pedro Shelf. This repository contains the bioinformatic and statistical scripts necessary to re-create analyses, as well as raw input data used to generate figures. FASTQ reads for both the Mock Community Mix (SUB9673624) and Natural Sample libraries (SUB10064804) can be obtained from NCBI’s SRA under BioProject PRJNA749297. The final Pseudo-nitzschia spp. 2bRAD reference library is hosted at https://dornsife.usc.edu/labs/carlslab/data/.

Files in this repository 
-----------

1. SeagrantCleanProcessingPipeline.txt: Annotated bioinformatic workflow for producing 2bRAD reference libraries (including read QC, assembly, and contaminant filtering). Custom scripts used in this analysis are listed below. To determine what each of these scripts does run them without arguments. Note that this pipeline was written for a cluster which uses a SLURM scheduler. 
	- Perl Scripts
	  - BcgI_Extract.pl
	  - TruncateFastq.pl
	  - uniquerOne.pl
	  - mergeUniq.pl
	  - BcgI_Extract_FullTag.pl
	  - splitFasta.pl
	  - taxfiles.pl
	  - TaxaOriginByBlast2.pl
	  - expression_compiler.pl
	  - trim2bRAD_1barcode_dedup.pl
	  - CombineExpression.pl
	- Shell Scripts
	  - organize.sh
	  - organizeRan.sh
	  - batch.sh
	  - batch2.sh
	  - SortContams.sh
	  - TagCount.sh
	
2. Analyses_ForGithub.R: Annotated R script for conducting statistical analyses of community and population-level data
	- Input file: Psub_bcgIN1sam.fasta.br.all.stripped.tab
	- Input file: AllContaminantCounts.csv
	- Input file: FPrate.csv
	- Input file: Accuracy.csv
	- Input file: NatSamLibs.csv
	- Input file: ARISA.csv
	- Input file: PpunCoverage.csv
	- Input file: CMH_analysis_ByYear.csv
	
	- Output file: Psubcontams.csv

