## Comp the FST/CMH genes - and then maybe just the FST and CMH back to the human and bee studies


## Grabbing these genes out of the gene_curations.R script and then going from here

## Path to gene_curations.R  is "/Users/arteen/Desktop/School/Thesis_Writing/scripts/populationGenomics/r/gene_curations.R"


sharedSNPs <- function( totalSNPs = 5099, sigSitesMine = 327, sigSitesOtherPaper = 56) {
  sites1 <- sample(totalSNPs, size = sigSitesMine, replace = F)
  sites2 <- sample(totalSNPs, size = sigSitesOtherPaper, replace = F)
  intersection_length <- length(intersect(sites1, sites2))
  return(intersection_length)
}

###################################################################################################


lowVHigh_Comp_df <- data.frame(matrix(nrow = 4,
                                         ncol = 3))
row.names(lowVHigh_Comp_df) <- c("Wang", "Shpigler",
                                    "Woodard", "Bralten")

colnames(lowVHigh_Comp_df) <- c("Overlapping", "MaxSimOverlapping",
                                   "Sim_95_percentile")


# LOW V HIGH - DOWN V UP
length(unique(pop_gen_UVD_genes$GeneId)) # 333

# Wang et al
length(wang_gene_curation_list)

LH_wang_olap <- pop_gen_UVD_genes[pop_gen_UVD_genes$GeneId %in% wang_gene_curation_list, ]

write.table(x = unique(LH_wang_olap$GeneId), file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Writing/Paper/Post_defence/popGen_WangOverlap.txt",
            sep = "/t", col.names = F, row.names = F, quote = F)
length(unique(LH_wang_olap$GeneId)) #5

wangSites <- replicate(10^4, sharedSNPs(totalSNPs = 6306, sigSitesMine = 333, sigSitesOtherPaper = 115))
max(wangSites) # 17
min(wangSites) # 0
quantile(wangSites, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))



lowVHigh_Comp_df[1,] <- as.integer(c(5,15,10))


# Shpigler et al
shp_gene_list

LH_shpig_olap <- pop_gen_UVD_genes[pop_gen_UVD_genes$GeneId %in% shp_gene_list, ]

length(unique(LH_shpig_olap$GeneId)) #28

write.table(x = unique(LH_shpig_olap$GeneId), file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Writing/Paper/Post_defence/popGen_ShpiglerOverlap.txt",
            sep = "/t", col.names = F, row.names = F, quote = F)


shpigSites <- replicate(10^4, sharedSNPs(totalSNPs = 6306, sigSitesMine = 333, sigSitesOtherPaper = 1057))
max(shpigSites) # 81
min(shpigSites) # 28
quantile(shpigSites, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))

lowVHigh_Comp_df[2,] <- as.integer(c(28,83,67))

# woodard et al
woodard_gene_list

LH_woodard_olap <- pop_gen_UVD_genes[pop_gen_UVD_genes$GeneId %in% woodard_gene_list, ]

length(unique(LH_woodard_olap$GeneId)) #2

write.table(x = unique(LH_woodard_olap$GeneId), file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Writing/Paper/Post_defence/popGen_WoodardOverlap.txt",
            sep = "/t", col.names = F, row.names = F, quote = F)

woodardSites <- replicate(10^4, sharedSNPs(totalSNPs = 6306, sigSitesMine = 333, sigSitesOtherPaper = 212))
max(woodardSites) # 25
min(woodardSites) # 0
quantile(woodardSites, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))

lowVHigh_Comp_df[3,] <- as.integer(c(2,24,17))

#Bralten et al

braltenGeneList <- bralten_csv$drosID

LH_bralten_olap <- pop_gen_UVD_genes[pop_gen_UVD_genes$GeneId %in% braltenGeneList, ]

length(unique(LH_bralten_olap$GeneId)) #2

write.table(x = unique(LH_bralten_olap$GeneId), file = "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Writing/Paper/Post_defence/popGen_BraltenOverlap.txt",
            sep = "/t", col.names = F, row.names = F, quote = F)


braltenSites <- replicate(10^4, sharedSNPs(totalSNPs = 5099, sigSitesMine = 333))
max(braltenSites) # 15
min(braltenSites) # 0
quantile(braltenSites, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))

lowVHigh_Comp_df[4,] <- as.integer(c(2,12,7))

write.csv(lowVHigh_Comp_df,
          "/Users/arteen/Desktop/School/Projects/SociabilityRNA/Writing/Paper/Post_defence/lowVHigh_popGen_compToOtherStudiesUPDATED.csv",
          quote = F,
          row.names = TRUE)
##########################################################################################
ControlVLow_Comp_df <- data.frame(matrix(nrow = 4,
                                         ncol = 3))
row.names(ControlVLow_Comp_df) <- c("Wang", "Shpigler",
                                    "Woodard", "Bralten")

colnames(ControlVLow_Comp_df) <- c("Overlapping", "MaxSimOverlapping",
                                   "Sim_95_percentile")

controlVhigh_Comp_df <- data.frame(matrix(nrow = 4,
                                         ncol = 3))
row.names(controlVhigh_Comp_df) <- c("Wang", "Shpigler",
                                    "Woodard", "Bralten")

colnames(controlVhigh_Comp_df) <- c("Overlapping", "MaxSimOverlapping",
                                   "Sim_95_percentile")

ancestorVLow_Comp_df <- data.frame(matrix(nrow = 4,
                                         ncol = 3))
row.names(ancestorVLow_Comp_df) <- c("Wang", "Shpigler",
                                    "Woodard", "Bralten")

colnames(ancestorVLow_Comp_df) <- c("Overlapping", "MaxSimOverlapping",
                                   "Sim_95_percentile")

ancestorVhigh_Comp_df <- data.frame(matrix(nrow = 4,
                                         ncol = 3))
row.names(ancestorVhigh_Comp_df) <- c("Wang", "Shpigler",
                                    "Woodard", "Bralten")

colnames(ancestorVhigh_Comp_df) <- c("Overlapping", "MaxSimOverlapping",
                                   "Sim_95_percentile")

# CONTROL VS LOW - DOWN V CONTROL
length(unique(pop_gen_CVD_genes$GeneId))

# Wang et al
length(wang_gene_curation_list)

CL_wang_olap <- pop_gen_CVD_genes[pop_gen_CVD_genes$GeneId %in% wang_gene_curation_list, ]

wangSites <- replicate(10^4, sharedSNPs(totalSNPs = 6306, sigSitesMine = 244, sigSitesOtherPaper = 115))
max(wangSites) # 17
min(wangSites) # 0
quantile(wangSites, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))



ControlVLow_Comp_df[1,] <- as.integer(c(13,15,8))


# Shpigler et al
shp_gene_list

CL_shpig_olap <- pop_gen_CVD_genes[pop_gen_CVD_genes$GeneId %in% shp_gene_list, ]

shpigSites <- replicate(10^4, sharedSNPs(totalSNPs = 6306, sigSitesMine = 244, sigSitesOtherPaper = 1057))
max(shpigSites) # 81
min(shpigSites) # 28
quantile(shpigSites, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))

# woodard et al
woodard_gene_list

woodardSites <- replicate(10^4, sharedSNPs(totalSNPs = 6306, sigSitesMine = 244, sigSitesOtherPaper = 212))
max(woodardSites) # 25
min(woodardSites) # 0
quantile(woodardSites, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))

#Bralten et al

bralten_csv

braltenSites <- replicate(10^4, sharedSNPs(totalSNPs = 5099, sigSitesMine = 244))
max(braltenSites) # 15
min(braltenSites) # 0
quantile(braltenSites, 
         probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))





## CONTROL VS HIGH - UP V CONTROL
length(unique(pop_gen_CVU_genes$GeneId))



## ANCESTOR VS LOW / DOWN
length(unique(pop_gen_AVL_genes$GeneId)) 

## ANCESTOR VS HIGH / UP
length(unique(pop_gen_AVH_genes$GeneId)) # 184 genes 


# for each of these 5 lists, want to see the overlap with the 4 different studies papers


# SO each of the 4 studies gene lists are 

# Wang et al
length(wang_gene_curation_list)

# Shpigler et al
shp_gene_list

# woodard et al
woodard_gene_list

#Bralten et al

bralten_csv




woodard_gene_subset <- all_selection_genes[all_selection_genes$geneID %in% woodard_gene_list, ]


## Basically


# Subset <-  popgengeneList[popgengeneList %in% paperGeneList,  ] - and then just the number

# and then after we'll run some simulations to see the number of overlapping by chance - and going to ignore the +/- 10% AT THE MOMENT





## Or just stick them all together into a table and make my life easier !
