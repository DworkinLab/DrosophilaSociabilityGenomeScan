library(UpSetR)
`%notin%` <- Negate(`%in%`)
## Pop Gen

## Vs this listdownVcontrol_subset

pop_gen_CVD_genes <- read.table(file = "/Users/arteen/Desktop/pop_gen_out/revisions/lowVcontrol/snpEff_genes.txt",
                                header = TRUE) #336 genes

length(unique(pop_gen_CVD_genes$GeneId)) # 245 genes 

head(downVcontrol_subset)

popgen_overlap <- downVcontrol_subset[downVcontrol_subset$geneID %in% pop_gen_CVD_genes$GeneId,]

popgen_overlap <- pop_gen_CVD_genes[pop_gen_CVD_genes$GeneId %in% downVcontrol_subset$geneID,]

CVL_DGE_overlap <- data.frame(geneID = unique(popgen_overlap$GeneId),
                              geneName = unique(popgen_overlap$GeneName))

# write.table(CVL_DGE_overlap,
#             file = "/Users/arteen/Desktop/pop_gen_out/revisions/gene_lists/dge_olap/CVL_geneList.csv",
#             quote = FALSE,
#             sep = ",",
#             row.names = FALSE,
#             col.names = TRUE)

high_and_moderate_index <- c(grep(1, pop_gen_CVD_genes$variants_impact_MODERATE),
                             grep(2, pop_gen_CVD_genes$variants_impact_MODERATE),
                             grep(3, pop_gen_CVD_genes$variants_impact_MODERATE))

high_and_moderate <- unique(pop_gen_CVD_genes$GeneId[high_and_moderate_index])

high_and_moderate_overlap <- downVcontrol_subset[downVcontrol_subset$geneID %in% high_and_moderate,]


#### Now for low vs High


pop_gen_UVD_genes <- read.table(file = "/Users/arteen/Desktop/pop_gen_out/revisions/lowVhigh/snpEff_genes.txt",
                                header = TRUE)

length(unique(pop_gen_UVD_genes$GeneId))


head(upVdown_subset)

popgen_overlap_UVD <- upVdown_subset[upVdown_subset$geneID %in% pop_gen_UVD_genes$GeneId,]

popgen_overlap_UVD <- pop_gen_UVD_genes[pop_gen_UVD_genes$GeneId %in% upVdown_subset$geneID,]

LVH_DGE_overlap <- data.frame(geneID = unique(popgen_overlap_UVD$GeneId),
                              geneName = unique(popgen_overlap_UVD$GeneName))
# 
# write.table(LVH_DGE_overlap,
#             file = "/Users/arteen/Desktop/pop_gen_out/revisions/gene_lists/LVH_geneList.csv",
#             quote = FALSE,
#             sep = ",",
#             row.names = FALSE,
#             col.names = TRUE)

high_and_moderate_index <- c(grep(1, pop_gen_UVD_genes$variants_impact_HIGH),
                             grep(1, pop_gen_UVD_genes$variants_impact_MODERATE),
                             grep(2, pop_gen_UVD_genes$variants_impact_MODERATE))

high_and_moderate <- pop_gen_UVD_genes$GeneId[high_and_moderate_index]

high_and_moderate_overlap <- upVdown_subset[upVdown_subset$geneID %in% high_and_moderate,]

#### Now for Control vs High


pop_gen_CVU_genes <- read.table(file = "/Users/arteen/Desktop/pop_gen_out/revisions/controlVhigh/snpEff_genes.txt",
                                header = TRUE) 

length(unique(pop_gen_CVU_genes$GeneId)) # 184 genes 




popgen_overlap_CVU <- controlVup_subset[controlVup_subset$geneID %in% pop_gen_CVU_genes$GeneId,]


popgen_overlap_CVU <- pop_gen_CVU_genes[pop_gen_CVU_genes$GeneId %in% controlVup_subset$geneID,]

CVH_DGE_overlap <- data.frame(geneID = unique(popgen_overlap_CVU$GeneId),
                              geneName = unique(popgen_overlap_CVU$GeneName))
# 
# write.table(CVH_DGE_overlap,
#             file = "/Users/arteen/Desktop/pop_gen_out/revisions/gene_lists/dge_olap/CVH_geneList.csv",
#             quote = FALSE,
#             sep = ",",
#             row.names = FALSE,
#             col.names = TRUE)

high_and_moderate_index <- c(grep(1, pop_gen_CVU_genes$variants_impact_HIGH),
                             grep(1, pop_gen_CVU_genes$variants_impact_MODERATE),
                             grep(2, pop_gen_CVU_genes$variants_impact_MODERATE),
                             grep(3, pop_gen_CVU_genes$variants_impact_MODERATE),
                             grep(5, pop_gen_CVU_genes$variants_impact_MODERATE))

high_and_moderate <- unique(pop_gen_CVU_genes$GeneId[high_and_moderate_index])

high_and_moderate_overlap <- controlVup_subset[controlVup_subset$geneID %in% high_and_moderate,]




## Ancestor vs Low
pop_gen_AVL_genes <- read.table(file = "/Users/arteen/Desktop/pop_gen_out/revisions/ancestorVlow/snpEff_genes.txt",
                                header = TRUE) 

length(unique(pop_gen_AVL_genes$GeneId)) 


downVcontrol_subset

popgen_overlap_CVU <- downVcontrol_subset[downVcontrol_subset$geneID %in% pop_gen_AVL_genes$GeneId,]


popgen_overlap_AVL <- pop_gen_AVL_genes[pop_gen_AVL_genes$GeneId %in% downVcontrol_subset$geneID,]

AVL_DGE_overlap <- data.frame(geneID = unique(popgen_overlap_AVL$GeneId),
                              geneName = unique(popgen_overlap_AVL$GeneName))

# write.table(AVL_DGE_overlap,
#             file = "/Users/arteen/Desktop/pop_gen_out/revisions/gene_lists/dge_olap/AVL_geneList.csv",
#             quote = FALSE,
#             sep = ",",
#             row.names = FALSE,
#             col.names = TRUE)



high_and_moderate_index <- c(grep(1, pop_gen_AVL_genes$variants_impact_HIGH),
                             grep(1, pop_gen_AVL_genes$variants_impact_MODERATE),
                             grep(2, pop_gen_AVL_genes$variants_impact_MODERATE),
                             grep(3, pop_gen_AVL_genes$variants_impact_MODERATE),
                             grep(4, pop_gen_AVL_genes$variants_impact_MODERATE))

high_and_moderate <- unique(pop_gen_AVL_genes$GeneId[high_and_moderate_index])

high_and_moderate_overlap <- controlVup_subset[controlVup_subset$geneID %in% high_and_moderate,]


popgen_overlap_AVL_DTU <- downVcontrol_dtu_subset[downVcontrol_dtu_subset$gene %in% pop_gen_AVL_genes$GeneId,]



## Ancestor vs High
pop_gen_AVH_genes <- read.table(file = "/Users/arteen/Desktop/pop_gen_out/revisions/ancestorVhigh/snpEff_genes.txt",
                                header = TRUE) 

length(unique(pop_gen_AVH_genes$GeneId)) # 184 genes 


downVcontrol_subset

popgen_overlap_CVU <- controlVup_subset[controlVup_subset$geneID %in% pop_gen_AVH_genes$GeneId,]

popgen_overlap_AVH <- pop_gen_AVH_genes[pop_gen_AVH_genes$GeneId %in% controlVup_subset$geneID,]

AVH_DGE_overlap <- data.frame(geneID = unique(popgen_overlap_AVH$GeneId),
                              geneName = unique(popgen_overlap_AVH$GeneName))

# write.table(AVH_DGE_overlap,
#             file = "/Users/arteen/Desktop/pop_gen_out/revisions/gene_lists/dge_olap/AVH_geneList.csv",
#             quote = FALSE,
#             sep = ",",
#             row.names = FALSE,
#             col.names = TRUE)




high_and_moderate_index <- c(grep(1, pop_gen_AVH_genes$variants_impact_MODERATE),
                             grep(2, pop_gen_AVH_genes$variants_impact_MODERATE),
                             grep(3, pop_gen_AVH_genes$variants_impact_MODERATE))

high_and_moderate <- unique(pop_gen_AVH_genes$GeneId[high_and_moderate_index])

high_and_moderate_overlap <- controlVup_subset[controlVup_subset$geneID %in% high_and_moderate,]


popgen_overlap_AVH_DTU <- controlVup_dtu_subset[controlVup_dtu_subset$gene %in% pop_gen_AVH_genes$GeneId,]






# pop gen genes overlapping with DTU genes

geneDictionary <- read.delim(file = "/Users/arteen/Downloads/FlyGeneDictionary.txt")
geneDictionary$FBgnID <- geneDictionary$validated_id

DTU_popgen_overlap_UVD <- upVdown_dtu_subset[upVdown_dtu_subset$transcript %in% pop_gen_UVD_genes$TranscriptId,]

DTU_popgen_overlap_CVD <- downVcontrol_dtu_subset[downVcontrol_dtu_subset$transcript %in% pop_gen_CVD_genes$TranscriptId,]

DTU_popgen_overlap_CVU <- controlVup_dtu_subset[controlVup_dtu_subset$transcript %in% pop_gen_CVU_genes$TranscriptId,]

DTU_popgen_overlap_AVL <- downVcontrol_dtu_subset[downVcontrol_dtu_subset$transcript %in% pop_gen_AVL_genes$TranscriptId,]

DTU_popgen_overlap_AVH <- controlVup_dtu_subset[controlVup_dtu_subset$transcript %in% pop_gen_AVH_genes$TranscriptId,]



####
LVH_DTU_overlap <- data.frame(transcript = DTU_popgen_overlap_UVD$transcript,
                              geneID = DTU_popgen_overlap_UVD$gene)
gene_name <- list()
for (i in 1:nrow(LVH_DTU_overlap)) {
  gene <- LVH_DTU_overlap$geneID[i]
  rowIndex <- which(geneDictionary$validated_id == gene)
  geneName_toAdd <- geneDictionary[rowIndex, 3]
  gene_name <- append(gene_name, geneName_toAdd)
}
geneName <- unlist(gene_name)
LVH_DTU_overlap <- cbind(LVH_DTU_overlap, geneName)
write.table(LVH_DTU_overlap,
            file = "/Users/arteen/Desktop/pop_gen_out/revisions/gene_lists/dtu_olap/LVH_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####

####
CVL_DTU_overlap <- data.frame(transcript = DTU_popgen_overlap_CVD$transcript,
                              geneID = DTU_popgen_overlap_CVD$gene)
gene_name <- list()
for (i in 1:nrow(CVL_DTU_overlap)) {
  gene <- CVL_DTU_overlap$geneID[i]
  rowIndex <- which(geneDictionary$validated_id == gene)
  geneName_toAdd <- geneDictionary[rowIndex, 3]
  if (length(geneName_toAdd) == 0){ # if the gene name isnt in the dictionary
    geneName_toAdd <- gene
  }
  gene_name <- append(gene_name, geneName_toAdd)
}
geneName <- unlist(gene_name)
CVL_DTU_overlap <- cbind(CVL_DTU_overlap, geneName)
write.table(CVL_DTU_overlap,
            file = "/Users/arteen/Desktop/pop_gen_out/revisions/gene_lists/dtu_olap/CVL_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####

####
CVH_DTU_overlap <- data.frame(transcript = DTU_popgen_overlap_CVU$transcript,
                              geneID = DTU_popgen_overlap_CVU$gene)
gene_name <- list()
for (i in 1:nrow(CVH_DTU_overlap)) {
  gene <- CVH_DTU_overlap$geneID[i]
  rowIndex <- which(geneDictionary$validated_id == gene)
  geneName_toAdd <- geneDictionary[rowIndex, 3]
  if (length(geneName_toAdd) == 0){ # if the gene name isnt in the dictionary
    geneName_toAdd <- gene
  }
  gene_name <- append(gene_name, geneName_toAdd)
}
geneName <- unlist(gene_name)
CVH_DTU_overlap <- cbind(CVH_DTU_overlap, geneName)
write.table(CVH_DTU_overlap,
            file = "/Users/arteen/Desktop/pop_gen_out/revisions/gene_lists/dtu_olap/CVH_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####

####
AVL_DTU_overlap <- data.frame(transcript = DTU_popgen_overlap_AVL$transcript,
                              geneID = DTU_popgen_overlap_AVL$gene)
gene_name <- list()
for (i in 1:nrow(AVL_DTU_overlap)) {
  gene <- AVL_DTU_overlap$geneID[i]
  rowIndex <- which(geneDictionary$validated_id == gene)
  geneName_toAdd <- geneDictionary[rowIndex, 3]
  if (length(geneName_toAdd) == 0){ # if the gene name isnt in the dictionary
    geneName_toAdd <- gene
  }
  gene_name <- append(gene_name, geneName_toAdd)
}
geneName <- unlist(gene_name)
AVL_DTU_overlap <- cbind(AVL_DTU_overlap, geneName)
write.table(AVL_DTU_overlap,
            file = "/Users/arteen/Desktop/pop_gen_out/revisions/gene_lists/dtu_olap/AVL_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####

####
AVH_DTU_overlap <- data.frame(transcript = DTU_popgen_overlap_AVH$transcript,
                              geneID = DTU_popgen_overlap_AVH$gene)
gene_name <- list()
for (i in 1:nrow(AVH_DTU_overlap)) {
  gene <- AVH_DTU_overlap$geneID[i]
  rowIndex <- which(geneDictionary$validated_id == gene)
  geneName_toAdd <- geneDictionary[rowIndex, 3]
  if (length(geneName_toAdd) == 0){ # if the gene name isnt in the dictionary
    geneName_toAdd <- gene
  }
  gene_name <- append(gene_name, geneName_toAdd)
}
geneName <- unlist(gene_name)
AVH_DTU_overlap <- cbind(AVH_DTU_overlap, geneName)
write.table(AVH_DTU_overlap,
            file = "/Users/arteen/Desktop/pop_gen_out/revisions/gene_lists/dtu_olap/AVH_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####



length(unique(DTU_popgen_overlap_UVD$transcript))
length(unique(DTU_popgen_overlap_CVD$transcript))
length(unique(DTU_popgen_overlap_CVU$transcript))






## Upset plot for the three pop gen lists

to_list <- list("Low versus High" = unique(pop_gen_UVD_genes$GeneId),
                "Control versus Low" = (unique(pop_gen_CVD_genes$GeneId)),
                "Control versus High" = unique(pop_gen_CVU_genes$GeneId),
                "Ancestor versus High" = unique(pop_gen_AVH_genes$GeneId),
                "Ancestor versus Low" = unique(pop_gen_AVL_genes$GeneId))



upset(data = fromList(to_list),
      nsets = 5,
      matrix.color = "red",
      main.bar.color = "black",
      mainbar.y.label = "Number of Genes",
      sets.bar.color = "Black",
      sets.x.label = "Total Set Size",
      point.size = 2.2,
      line.size = 0.7,
      order.by = c("freq"),
      decreasing = c(T),
      group.by = "degree")



############# Writing Gene Lists #############
## Now going to make csvs of gene lists
## going to include both gene name and FBgnID


## CMH First

AVL_CMH_genes <- read.table(file = "/Users/arteen/Desktop/pop_gen_out/revisions/ancestorVlow/cmh/snpEff_genes.txt",
                                header = TRUE) 

AVL_CMH_geneList <- data.frame(geneID = unique(AVL_CMH_genes$GeneId),
                               geneName = unique(AVL_CMH_genes$GeneName))


write.table(AVL_CMH_geneList,
            file = "/Users/arteen/Desktop/pop_gen_out/revisions/gene_lists/cmh/AVL_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####


AVH_CMH_genes <- read.table(file = "/Users/arteen/Desktop/pop_gen_out/revisions/ancestorVhigh/cmh/snpEff_genes.txt",
                            header = TRUE) 

AVH_CMH_geneList <- data.frame(geneID = unique(AVH_CMH_genes$GeneId),
                               geneName = unique(AVH_CMH_genes$GeneName))

head(AVH_CMH_geneList)
write.table(AVH_CMH_geneList,
            file = "/Users/arteen/Desktop/pop_gen_out/revisions/gene_lists/cmh/AVH_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####


LVH_CMH_genes <- read.table(file = "/Users/arteen/Desktop/pop_gen_out/revisions/lowVHigh/cmh/snpEff_genes.txt",
                            header = TRUE) 

LVH_CMH_geneList <- data.frame(geneID = unique(LVH_CMH_genes$GeneId),
                               geneName = unique(LVH_CMH_genes$GeneName))

head(LVH_CMH_geneList)
write.table(LVH_CMH_geneList,
            file = "/Users/arteen/Desktop/pop_gen_out/revisions/gene_lists/cmh/LVH_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####


CVL_CMH_genes <- read.table(file = "/Users/arteen/Desktop/pop_gen_out/revisions/lowVcontrol/cmh/snpEff_genes.txt",
                            header = TRUE) 

CVL_CMH_geneList <- data.frame(geneID = unique(CVL_CMH_genes$GeneId),
                               geneName = unique(CVL_CMH_genes$GeneName))

head(CVL_CMH_geneList)
write.table(CVL_CMH_geneList,
            file = "/Users/arteen/Desktop/pop_gen_out/revisions/gene_lists/cmh/CVL_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####


CVH_CMH_genes <- read.table(file = "/Users/arteen/Desktop/pop_gen_out/revisions/controlVhigh/cmh/snpEff_genes.txt",
                            header = TRUE) 

CVH_CMH_geneList <- data.frame(geneID = unique(CVH_CMH_genes$GeneId),
                               geneName = unique(CVH_CMH_genes$GeneName))

head(CVH_CMH_geneList)
write.table(CVH_CMH_geneList,
            file = "/Users/arteen/Desktop/pop_gen_out/revisions/gene_lists/cmh/CVH_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####




AVL_FST_genes <- read.table(file = "/Users/arteen/Desktop/pop_gen_out/revisions/ancestorVlow/fst/snpEff_genes.txt",
                            header = TRUE) 

AVL_FST_geneList <- data.frame(geneID = unique(AVL_FST_genes$GeneId),
                               geneName = unique(AVL_FST_genes$GeneName))

head(AVL_FST_geneList)
write.table(AVL_FST_geneList,
            file = "/Users/arteen/Desktop/pop_gen_out/revisions/gene_lists/fst/AVL_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####




AVH_FST_genes <- read.table(file = "/Users/arteen/Desktop/pop_gen_out/revisions/ancestorVhigh/fst/snpEff_genes.txt",
                            header = TRUE) 

AVH_FST_geneList <- data.frame(geneID = unique(AVH_FST_genes$GeneId),
                               geneName = unique(AVH_FST_genes$GeneName))

head(AVH_FST_geneList)
write.table(AVH_FST_geneList,
            file = "/Users/arteen/Desktop/pop_gen_out/revisions/gene_lists/fst/AVH_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####




LVH_FST_genes <- read.table(file = "/Users/arteen/Desktop/pop_gen_out/revisions/lowVHigh/fst/snpEff_genes.txt",
                            header = TRUE) 

LVH_FST_geneList <- data.frame(geneID = unique(LVH_FST_genes$GeneId),
                               geneName = unique(LVH_FST_genes$GeneName))

head(LVH_FST_geneList)
write.table(LVH_FST_geneList,
            file = "/Users/arteen/Desktop/pop_gen_out/revisions/gene_lists/fst/LVH_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####




CVL_FST_genes <- read.table(file = "/Users/arteen/Desktop/pop_gen_out/revisions/lowVcontrol/fst/snpEff_genes.txt",
                            header = TRUE) 

CVL_FST_geneList <- data.frame(geneID = unique(CVL_FST_genes$GeneId),
                               geneName = unique(CVL_FST_genes$GeneName))

head(CVL_FST_geneList)
write.table(CVL_FST_geneList,
            file = "/Users/arteen/Desktop/pop_gen_out/revisions/gene_lists/fst/CVL_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####




CVH_FST_genes <- read.table(file = "/Users/arteen/Desktop/pop_gen_out/revisions/controlVhigh/fst/snpEff_genes.txt",
                            header = TRUE) 

CVH_FST_geneList <- data.frame(geneID = unique(CVH_FST_genes$GeneId),
                               geneName = unique(CVH_FST_genes$GeneName))

head(CVH_FST_geneList)
write.table(CVH_FST_geneList,
            file = "/Users/arteen/Desktop/pop_gen_out/revisions/gene_lists/fst/CVH_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####


# Now for the overlap between FST and CMH


LVH_CMH_FST_Overlap <- data.frame(geneID = unique(pop_gen_UVD_genes$GeneId),
                                  geneName = unique(pop_gen_UVD_genes$GeneName))
head(LVH_CMH_FST_Overlap)
write.table(LVH_CMH_FST_Overlap,
            file = "/Users/arteen/Desktop/pop_gen_out/revisions/gene_lists/cmhFstOlap/LVH_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
##
AVL_CMH_FST_Overlap <- data.frame(geneID = unique(pop_gen_AVL_genes$GeneId),
                                  geneName = unique(pop_gen_AVL_genes$GeneName))
head(AVL_CMH_FST_Overlap)
write.table(AVL_CMH_FST_Overlap,
            file = "/Users/arteen/Desktop/pop_gen_out/revisions/gene_lists/cmhFstOlap/AVL_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####
AVH_CMH_FST_Overlap <- data.frame(geneID = unique(pop_gen_AVH_genes$GeneId),
                                  geneName = unique(pop_gen_AVH_genes$GeneName))
head(AVH_CMH_FST_Overlap)
write.table(AVH_CMH_FST_Overlap,
            file = "/Users/arteen/Desktop/pop_gen_out/revisions/gene_lists/cmhFstOlap//AVH_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)

####
CVL_CMH_FST_Overlap <- data.frame(geneID = unique(pop_gen_CVD_genes$GeneId),
                                  geneName = unique(pop_gen_CVD_genes$GeneName))
head(CVL_CMH_FST_Overlap)
write.table(CVL_CMH_FST_Overlap,
            file = "/Users/arteen/Desktop/pop_gen_out/revisions/gene_lists/cmhFstOlap//CVL_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)
####
CVH_CMH_FST_Overlap <- data.frame(geneID = unique(pop_gen_CVU_genes$GeneId),
                                  geneName = unique(pop_gen_CVU_genes$GeneName))
head(CVH_CMH_FST_Overlap)
write.table(CVH_CMH_FST_Overlap,
            file = "/Users/arteen/Desktop/pop_gen_out/revisions/gene_lists/cmhFstOlap/CVH_geneList.csv",
            quote = FALSE,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)



# OLAP w DGE




