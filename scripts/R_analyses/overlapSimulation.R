## Updated simulation for overlap between our PG results and other sociability genomics studies.


# updated by ID, April 21th, 2025
# This includes an adjustment for the proportion of genes that have orthologs between species under consideration.

# Bratlen et al 2021
# 56 unique genes from their GWAS. 40 of which with orthologs in Dm
# 19,427 protein coding genes in humans

# For our results from pop gen it is 169 unique genes that have human orthologs according to 
# 15094 protein coding genes Drosophila


# Only considering common genes
# This is probably a more sensible check.
# This uses the total gene number from:
# https://flybase.org/docs/releasenotes.tx
# https://www.ncbi.nlm.nih.gov/genome/47?genome_assembly_id=204923
#  Genome R6.48, with 17896 annotated gene models

# Of the 15094 protein coding genes, there are 9096 with human orthologs according to dmel_human_orthologs_disease_fb_2025_02.tsv (downoaded 21/04/25)
# https://s3ftp.flybase.org/releases/FB2025_02/precomputed_files/orthologs/dmel_human_orthologs_disease_fb_2025_02.tsv.gz

propOrth <- 9096/15094


# note I am using the total number of annotated genes as we did evaluate all genes in Drosophila, coding and non-coding. So an assumption that the amount of orthology is shared for protein and other genes.

sharedGenes <- function( totalGenes_species1 = 19427, sigGenes_species1 = 40, totalGenes_species2 = 17896, sigGenes_species2 = 169, propOrthologs = propOrth) {
  
  totalGenes_species1 <- round(totalGenes_species1*propOrthologs)
  sigGenes_species2 <- round(sigGenes_species2*propOrthologs)
  
  sites1 <- sample(totalGenes_species1, size = sigGenes_species1, replace = F)
  sites2 <- sample(totalGenes_species2, size = sigGenes_species2, replace = F)
  intersection_length <- length(intersect(sites1, sites2))
  return(intersection_length)
}



# total SNPs is all dros genes we are checking? and sig sites


# This is the formula for doing this, straight from Ians code

HowManySharedSites <- replicate(10^5, sharedGenes())

max(HowManySharedSites) # the maximum number of shared sites among all the simulations

quantile(HowManySharedSites, 
         probs = c(0.5, 0.9, 0.95, 0.99, 1))

mean(HowManySharedSites > 0) # what proportion of simulations have any sites in common.


# Shpigler 10317 total genes (from supplemental table 1 of Shpigler). 1,057 DGE (main body of text). Of which there were 6651 orthologs with Drosophila melanogaster (Supp table 1 Shpigler).

propOrth2 <- 6651/10317

sharedGenes2 <- function( totalGenes_species1 = 10317, sigGenes_species1 = 1057, totalGenes_species2 = 17896, sigGenes_species2 = 169, propOrthologs = propOrth2) {
  
  totalGenes_species2 <- round(totalGenes_species2*propOrthologs)
  sigGenes_species2 <- round(sigGenes_species2*propOrthologs)
  
  sites1 <- sample(totalGenes_species1, size = sigGenes_species1, replace = F)
  sites2 <- sample(totalGenes_species2, size = sigGenes_species2, replace = F)
  intersection_length <- length(intersect(sites1, sites2))
  return(intersection_length)
}

# This is the formula for doing this, straight from Ians code

HowManySharedSites <- replicate(10^5, sharedGenes2())

max(HowManySharedSites) # the maximum number of shared sites among all the simulations

quantile(HowManySharedSites, 
         probs = c(0.5, 0.9, 0.95, 0.99, 1))

mean(HowManySharedSites > 0) # what proportion of simulations have any sites in common.



## Comparison for candidate genes identified by DGE and pop gen in D.mel

sharedGenes3 <- function( totalGenes_species1 = 17896, sigGenes_species1 = 328, totalGenes_species2 = 17896, sigGenes_species2 = 169) {
  
  sites1 <- sample(totalGenes_species1, size = sigGenes_species1, replace = F)
  sites2 <- sample(totalGenes_species2, size = sigGenes_species2, replace = F)
  intersection_length <- length(intersect(sites1, sites2))
  return(intersection_length)
}

HowManySharedSites <- replicate(10^5, sharedGenes3())

max(HowManySharedSites) # the maximum number of shared sites among all the simulations

quantile(HowManySharedSites, 
         probs = c(0.5, 0.9, 0.95, 0.99, 1))

mean(HowManySharedSites > 0) # what proportion of simulations have any sites in common.
