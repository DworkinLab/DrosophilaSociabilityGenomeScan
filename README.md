# Data and scripts repository for Rosenbaum, Torabi-Marashi et al. 2025.

This repository contains all data (except raw sequencing files) and scripts associated with the study:

Rosenbaum, Torabi-Marashi et al. 2025. Evolution of sociability: Genome scans and gene validation. *Evolution* Accepted.

Raw Sequence data for this study is available on [NCBI SRA bioproject  PRJNA1311514](https://www.ncbi.nlm.nih.gov/sra/PRJNA1311514).

A static version of this repository (frozen at time of final acceptance) is available on [Figshare](https://doi.org/10.6084/m9.figshare.30500915.v1). **Please note:** figshare has a complete copy of all of the `.sync files` (file sizes are too big for github).

**Please note**: If you wish to use this pipeline for the genomic analysis, and you are new to such pipelines, before using the scripts before, you will likely find the README for the data and scripts associated with [Audet et al. 2024](https://github.com/DworkinLab/Audet_etal_Evolution_2024) a simple overview of the pipeline. Then the scripts below will likely be more useful for replicating the specific analyses from this study.


### Data

#### Phenotypic Data

`JR_SociabilityScores_AllCrosses_formatted.csv`  Behavioural data (sociability) for the RNAi experiments to test candidate genes. These are analyzed in the script `JR_SA_SociabilityRNAiCrosses_March2025_ID.Rmd`

Day: Day of experiment within an experimental block 
Time: Time (in 15 minute intervals) of phenotyping. So 1 corresponds to the first measure at 15 minutes into the behavioural assessment, 2 is 30 minutes etc...

Arena: The arena (unit of sampling) unique to a trial within a day.

Treatment: RNAi or control treatment

Soc_index: Measure of social aggregation

ExperimentalBlock: Blocking variable for the entire experiment (across all RNAi strains)

StrainID: Bloomington Drosophila Stock Center identifer


Additional meta-data for the experiment is in file `Construct_meta_data_JR_correct_EXPblocks.csv`.

stock_number: Bloomington Drosophila Stock Center identifer

TestDates: Beginning date for experiment (corresponding to Day 1)

TRiP: Whether the RNAi strain was from the TRiP or VDRC collection (for this study all were TRiP).

vector: Which vector was used to make this particular strain

genes: gene name

FBgn: Flybase gene identifier

Insertion location: Location of transgene insertion

Insertion Chromosom: Chromosomal location of transgene insertion

full_strain_genotype: Formal genotypic information

Notes: Any other notes about the strains either from Bloomington or our own work.

#### Allele counts data

`.sync` files of allele counts and coverage by samples. Because of file size limitations only the .sync files for the X chromosomes are part of this repo. Sync files for autosomes need to be regenerated from the analysis pipeline.See below for more information on how they are generated (mostly via the shell scripts) and the Rscripts for analysis. Alternatively, go to the static version of the repository on [figshare](https://doi.org/10.6084/m9.figshare.30500915.v1) which has all of the sync files.

#### qRT-PCR

All `.xlsx` files beginning with "JR_" as well as `rg_JRR_5122025_data.xls` are the  output from the qPCR runs. These are analyzed in the script `SociabilityPopGen_qPCR_results_JRR.R`

### Intermediate Data

#### snp_eff_cmh_fst_overlaps

SnpEff outputs and bed files from CMH and FST outputs (to be used for SnpEff) in all contrasts.

#### AlleleF_plotting.txt

Plotting text file to pull relevant gene names. Note that anywhere in the R scripts where alleleF_plotting.txt is being used, that is just snp_eff.txt for a given contrast, just with the snpeff header removed.



### Output

#### ACER and FST outputs for low versus high contrast

#### gene_lists

gene lists from all analyses as well as overlapping with gene lists found in Torabi-Marashi et al. 2025

### Scripts

#### R_Analyses

`JR_SA_SociabilityRNAiCrosses_March2025_ID.Rmd` - Analysis of behavioural data for RNAi experiments

`alleleFrequencyPlottingLoop__Sept8_2025.Rmd` - Allele Frequency plotting loop to plot identified SNPs into single PDF

`cmh_updated.R` - CMH test using the ACER package in R (split between autosomes and X chromosome)

`dokb_snps_AF_plotting.Rmd` - alternate Allele Frequency plotting for SNPs in *dokb* gene

`Figure1_allele_frequency_plotting.R` - alternate Code used to help generate allele frequency plot in figure 1

`fst_functions.R` - helper functions to plot FST

`gene_curations.R` - overlaps with Torabi-Marashi et al 2025 and between FST and CMH

`geneCurations_popGen_PostSubmit.R` - overlapping genes with other studies

`overlapSimulation.R` - simulation to test how likely the number of observed genes would overlap by chance

`read_sync_SpaceDelimited.R` - helper function to read sync files with space delimitation

`revisions_paper_FST_Plot.R` - Plotting FST

`updated_all_points_dokb_snps_AF_plotting.Rmd` - Allele Frequency plotting for SNPs in *dokb* gene

`updated_All_points_Figure1_allele_frequency_plotting.R` - Code used to help generate allele frequency plot in figure 1

#### Unix

Unix files to generate from BAM to sync files. All scripts are light and their titles are self explanatory to their function



### Intermediate Data

#### snp_eff_cmh_fst_overlaps - SnpEff outputs and bed files from CMH and FST outputs (to be used for SnpEff) in all contrasts

#### AlleleF_plotting.txt - plotting text file to pull relevant gene names
Note that anywhere in the R scripts where alleleF_plotting.txt is being used, that is just snp_eff.txt for a given contrast, just with the snpeff header removed



### Unix

Unix files for genomics analysis going from raw `.fastq` files to generate BAM to sync files. Titles are self explanatory to their function.

The relative order of the scripts (from `.fastq` to `.bam` to `.sync`) is:

#### QC & trimming adaptor and low quality sequence from sequence data

`fastqc_dna.sh` 

`trim_rename.sh` 

`trim_fastqc.sh` 

#### Index and map reads
`index_bwa.sh` 

`bwa_map.sh` 

#### Samtools processing

`sam_to_bam.sh`

`filter.sh`

`bam_fastqc.sh`

`bam_quality_check.sh`

`extractCoreGenomes.sh`

`sort_byName.sh`

`sort_byCoord.sh`

`fixmated.sh`

`markdup.sh`

#### Picard and GATK
`picard_RG.sh`

`gatk_mark_indels.sh`

`gatk_realign.sh`


#### Merge replicates and create mpileup
`ANC_merge.sh`

`U_merge.sh`

`C_merge.sh`

`D_merge.sh`

`make_mpileup.sh`

#### SNP Calling

`run_PoolSNP.sh`

`repeat_masker.sh`

`detectIndels_Kapun.sh`


#### Mpileup to Sync
`mp2sync.sh`

`subset_syncByVCF.sh`

#### FST
`run_fst.sh`

#### Annotation
`local_snp_eff.sh`

`snpEff.sh`

#### Misc helper function
`subset_vcf_byBed.sh`


### R_Analyses

`cmh_updated.R` - CMH test using the ACER package in R (split between autosomes and X chromosome)

`revisions_paper_FST_Plot.R` - Plotting FST

`fst_functions.R` - helper functions to plot FST

`alleleFrequencyPlottingLoop__Sept8_2025.Rmd` - Allele Frequency plotting loop to plot identified SNPs into single PDF

`gene_curations.R` - overlaps with Torabi-Marashi et al 2025 and between FST and CMH

`geneCurations_popGen_PostSubmit.R` - overlapping genes with other studies

`overlapSimulation.R` - simulation to test how likely the number of observed genes would overlap by chance

`read_sync_SpaceDelimited.R` - helper function to read sync files with space delimitation

`updated_all_points_dokb_snps_AF_plotting.Rmd` - Allele Frequency plotting for SNPs in dokb gene

`updated_All_points_Figure1_allele_frequency_plotting.R` - Code used to help generate allele frequency plot in figure 1

`dokb_snps_AF_plotting.Rmd` - alternate Allele Frequency plotting for SNPs in dokb gene

`Figure1_allele_frequency_plotting.R` - alternate Code used to help generate allele frequency plot in figure 1
