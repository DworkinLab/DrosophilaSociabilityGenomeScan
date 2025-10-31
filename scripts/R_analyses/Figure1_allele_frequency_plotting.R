## Allele Freq plots for paper - post revision

library(data.table)
library(stringr)
library(ggplot2)
library(ggtext)
library(cowplot)

# Some of the subsetting is unnecessary because I already have list of genes

pos_table <- fread(file = "/Users/arteen/Desktop/pop_gen_out/revisions/lowVhigh/alleleF_plotting.txt")
pos_table <- pos_table[,c(1,2,4,5)]
colnames(pos_table) <- c("CHROM", "POS", "REF", "ALT")
pos_table <- as.data.frame(pos_table)
pos_table$match <- paste0(pos_table$CHROM, "_", pos_table$POS)

syncFile <- fread(file = "/Users/arteen/Desktop/pop_gen_out/actual_sync/unmerged_subsetted.sync")
syncFile <- as.data.frame(syncFile)
colnames(syncFile) <- c("CHROM", "POS", "REF",
                        "Ancestor1","Ancestor2","Ancestor3","Ancestor4",
                        "Control1","Control2","Control3","Control4",
                        "Down1","Down2","Down3","Down4", 
                        "Up1","Up2","Up3","Up4")

syncFile$match <- paste0(syncFile$CHROM, "_", syncFile$POS)
syncFile <- syncFile[(syncFile$match %in% pos_table$match), ]
syncFile <- as.data.frame(syncFile) # wants to stay a data.table for some reason





# Give it the named vector
plot_allele_f <- function(name, is_last = FALSE, is_corner = FALSE) {
  #name <- name1
  
  gene_name <- str_split_i(name, ";", i = 2)
  allele <- str_split_i(name, ";", i = 1)
  
  syncRow <- syncFile[which(syncFile$match == allele),]
  # Average and SE of ancestor
  ancestralDF <- data.frame(matrix(nrow = 4, ncol = 5))
  colnames(ancestralDF) <- c("A", "T", "C", "G", "AF")
  row.names(ancestralDF) <- colnames(syncRow)[4:7]
  
  for (j in 4:7) { # for ancestral
    a <- as.integer(str_split_i(syncRow[, j], ":", i = 1))
    t <- as.integer(str_split_i(syncRow[, j], ":", i = 2))
    c <- as.integer(str_split_i(syncRow[, j], ":", i = 3))
    g <- as.integer(str_split_i(syncRow[, j], ":", i = 4))
    
    ancestralDF[c(j - 3), 1] <- a
    ancestralDF[c(j - 3), 2] <- t
    ancestralDF[c(j - 3), 3] <- c
    ancestralDF[c(j - 3), 4] <- g
    
    if (syncRow$REF == "A"){
      ancestralDF[c(j - 3), 5] <- (t + c + g) / (a + t + c + g)
    } else if(syncRow$REF == "T"){
      ancestralDF[c(j - 3), 5] <- (a + c + g) / (a + t + c + g)
    } else if(syncRow$REF == "C"){
      ancestralDF[c(j - 3), 5] <- (a + t + g) / (a + t + c + g)
    } else if(syncRow$REF == "G"){
      ancestralDF[c(j - 3), 5] <- (a + t + c) / (a + t + c + g)
    } else { 
      cat("There is no reference") 
      next
    }
  }
  
  ancestorMean <- mean(ancestralDF$AF)
  ancestorSE <- sd(ancestralDF$AF) / sqrt(nrow(ancestralDF))
  
  # Control
  controlDF <- data.frame(matrix(nrow = 4, ncol = 5))
  colnames(controlDF) <- c("A", "T", "C", "G", "AF")
  row.names(controlDF) <- colnames(syncRow)[8:11]
  
  for (j in 8:11) { # for ancestral
    a <- as.integer(str_split_i(syncRow[, j], ":", i = 1))
    t <- as.integer(str_split_i(syncRow[, j], ":", i = 2))
    c <- as.integer(str_split_i(syncRow[, j], ":", i = 3))
    g <- as.integer(str_split_i(syncRow[, j], ":", i = 4))
    
    controlDF[c(j - 7), 1] <- a
    controlDF[c(j - 7), 2] <- t
    controlDF[c(j - 7), 3] <- c
    controlDF[c(j - 7), 4] <- g
    
    if (syncRow$REF == "A"){
      controlDF[c(j - 7), 5] <- (t + c + g) / (a + t + c + g)
    } else if(syncRow$REF == "T"){
      controlDF[c(j - 7), 5] <- (a + c + g) / (a + t + c + g)
    } else if(syncRow$REF == "C"){
      controlDF[c(j - 7), 5] <- (a + t + g) / (a + t + c + g)
    } else if(syncRow$REF == "G"){
      controlDF[c(j - 7), 5] <- (a + t + c) / (a + t + c + g)
    } else { 
      cat("There is no reference") 
      next
    }
  }
  
  controlMean <- mean(controlDF$AF)
  controlSE <- sd(controlDF$AF) / sqrt(nrow(controlDF))
  
  ## Down
  downDF <- data.frame(matrix(nrow = 4, ncol = 5))
  colnames(downDF) <- c("A", "T", "C", "G", "AF")
  row.names(downDF) <- colnames(syncRow)[12:15]
  
  for (j in 12:15) { # for ancestral
    a <- as.integer(str_split_i(syncRow[, j], ":", i = 1))
    t <- as.integer(str_split_i(syncRow[, j], ":", i = 2))
    c <- as.integer(str_split_i(syncRow[, j], ":", i = 3))
    g <- as.integer(str_split_i(syncRow[, j], ":", i = 4))
    
    downDF[c(j - 11), 1] <- a
    downDF[c(j - 11), 2] <- t
    downDF[c(j - 11), 3] <- c
    downDF[c(j - 11), 4] <- g
    
    if (syncRow$REF == "A"){
      downDF[c(j - 11), 5] <- (t + c + g) / (a + t + c + g)
    } else if(syncRow$REF == "T"){
      downDF[c(j - 11), 5] <- (a + c + g) / (a + t + c + g)
    } else if(syncRow$REF == "C"){
      downDF[c(j - 11), 5] <- (a + t + g) / (a + t + c + g)
    } else if(syncRow$REF == "G"){
      downDF[c(j - 11), 5] <- (a + t + c) / (a + t + c + g)
    } else { 
      cat("There is no reference") 
      next
    }
  }
  
  downMean <- mean(downDF$AF)
  downSE <- sd(downDF$AF) / sqrt(nrow(downDF))
  
  
  ## Up
  upDF <- data.frame(matrix(nrow = 4, ncol = 5))
  colnames(upDF) <- c("A", "T", "C", "G", "AF")
  row.names(upDF) <- colnames(syncRow)[16:19]
  
  for (j in 16:19) { # for ancestral
    a <- as.integer(str_split_i(syncRow[, j], ":", i = 1))
    t <- as.integer(str_split_i(syncRow[, j], ":", i = 2))
    c <- as.integer(str_split_i(syncRow[, j], ":", i = 3))
    g <- as.integer(str_split_i(syncRow[, j], ":", i = 4))
    
    upDF[c(j - 15), 1] <- a
    upDF[c(j - 15), 2] <- t
    upDF[c(j - 15), 3] <- c
    upDF[c(j - 15), 4] <- g
    
    if (syncRow$REF == "A"){
      upDF[c(j - 15), 5] <- (t + c + g) / (a + t + c + g)
    } else if(syncRow$REF == "T"){
      upDF[c(j - 15), 5] <- (a + c + g) / (a + t + c + g)
    } else if(syncRow$REF == "C"){
      upDF[c(j - 15), 5] <- (a + t + g) / (a + t + c + g)
    } else if(syncRow$REF == "G"){
      upDF[c(j - 15), 5] <- (a + t + c) / (a + t + c + g)
    } else { 
      cat("There is no reference") 
      next
    }
  }
  
  upMean <- mean(upDF$AF)
  upSE <- sd(upDF$AF) / sqrt(nrow(upDF))

  # Keep the individual plotting just in case  
  alleleF <- c(ancestorMean,
               controlMean,
               downMean,
               upMean)
               # controlDF[1, 5], controlDF[2, 5], controlDF[3, 5], controlDF[4, 5],
               # downDF[1, 5], downDF[2, 5], downDF[3, 5], downDF[4, 5],
               # upDF[1, 5], upDF[2, 5], upDF[3, 5], upDF[4, 5])
  
  treatment <- c("Ancestor",
                 "Control",
                 "Low",
                 "High")
                 # rep("Control", 4),
                 # rep("Low", 4),
                 # rep("High", 4))
  
  
  
  # treatment <- c("Ancestor",
  #              "Control_1", "Control_2", "Control_3", "Control_4",
  #              "Low_1", "Low_2", "Low_3", "Low_4", 
  #              "High_1", "High_2", "High_3", "High_4")
  
  SE_list <- c(ancestorSE,
               controlSE,
               downSE,
               upSE)
               # rep(NA, 12))
  
  alleleF_df <- data.frame(alleleFreq = alleleF,
                           SE = SE_list,
                           Treatment = treatment)
  
  alleleF_df$Treatment <- factor(alleleF_df$Treatment, levels = c("Ancestor", "Control", "Low", "High"))
  
  max_height <- max(alleleF_df$alleleFreq) + max(alleleF_df$SE) + 0.1
  
  ancestral_df <- alleleF_df[1,]
  alleleF_df <- alleleF_df[-1,]
  
  
  # Wont work currently with mean of lineages
  # alleleF_df$lineage <- c("1", "2", "3", "4", 
  #                         "1", "2", "3", "4", 
  #                         "1", "2", "3", "4")

  #width = 0.45, height = 0,
  pd <-  position_dodge(width = 0.45)
  
  
  if (is_last){
    plot_last <- ggplot() +
      geom_point(data = ancestral_df,
                 aes(x = Treatment, y = alleleFreq, color = Treatment),
                 size = 2.5) +
      geom_errorbar(data = ancestral_df, 
                    aes(x = Treatment, y = alleleFreq, ymin = alleleFreq - SE, ymax = alleleFreq + SE, color = Treatment),
                    width = 0.15, size = 0.7, position = pd, alpha = 1) +
      geom_point(data = alleleF_df,
                 aes(x = Treatment, y = alleleFreq, color = Treatment), # shape = lineage),
                 size = 2.5, alpha = 0.7) +
      geom_errorbar(data = alleleF_df, 
                    aes(x = Treatment, y = alleleFreq, ymin = alleleFreq - SE, ymax = alleleFreq + SE, color = Treatment),
                    width = 0.15, size = 0.7, position = pd, alpha = 1) +
      theme_classic() +
      labs(title = gene_name, subtitle = allele) +
      theme(plot.title = element_text(hjust = 0.5, face = "italic", size = 15),
            plot.subtitle = element_text(size = 10),
            axis.text = element_text(face = "bold"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank())
    return(plot_last)
  }
  
  if (is_corner){
    plot_corner <- ggplot() +
      geom_point(data = ancestral_df,
                  aes(x = Treatment, y = alleleFreq, color = Treatment),
                  size = 2.5) +
      geom_errorbar(data = ancestral_df, 
                    aes(x = Treatment, y = alleleFreq, ymin = alleleFreq - SE, ymax = alleleFreq + SE, color = Treatment),
                    width = 0.15, size = 0.7, position = pd, alpha = 1) +
      geom_point(data = alleleF_df,
                  aes(x = Treatment, y = alleleFreq, color = Treatment), # shape = lineage),
                    size = 2.5, alpha = 0.7) +
      geom_errorbar(data = alleleF_df, 
                    aes(x = Treatment, y = alleleFreq, ymin = alleleFreq - SE, ymax = alleleFreq + SE, color = Treatment),
                    width = 0.15, size = 0.7, position = pd, alpha = 1) +
      theme_classic() +
      labs(title = gene_name, subtitle = allele) +
      theme(plot.title = element_text(hjust = 0.5, face = "italic", size = 15),
            plot.subtitle = element_text(size = 10),
            axis.text = element_text(face = "bold"),
            legend.position = "none") +
      xlab("Treatment") +
      ylab("Allele Frequency") 
    return(plot_corner)
  } 
  
  plot_else <- ggplot() +
    geom_point(data = ancestral_df,
               aes(x = Treatment, y = alleleFreq, color = Treatment),
               size = 2.5) +
    geom_errorbar(data = ancestral_df, 
                  aes(x = Treatment, y = alleleFreq, ymin = alleleFreq - SE, ymax = alleleFreq + SE, color = Treatment),
                  width = 0.15, size = 0.7, position = pd, alpha = 1) +
    geom_point(data = alleleF_df,
               aes(x = Treatment, y = alleleFreq, color = Treatment), # shape = lineage),
               size = 2.5, alpha = 0.7) +
    geom_errorbar(data = alleleF_df, 
                  aes(x = Treatment, y = alleleFreq, ymin = alleleFreq - SE, ymax = alleleFreq + SE, color = Treatment),
                  width = 0.15, size = 0.7, position = pd, alpha = 1) +
    theme_classic() +
    labs(title = gene_name, subtitle = allele) +
    theme(plot.title = element_text(hjust = 0.5, face = "italic", size = 15),
          plot.subtitle = element_text(size = 10),
          axis.text = element_text(face = "bold"),
          legend.position = "none",
          axis.title.x=element_blank(),
          axis.title.y=element_blank())
  return(plot_else)
}


#genes in order

name1 <- "3L_4676265;axo"
name2 <- "3L_3741911;CG32264"
name3 <- "X_14347162;dpr8"
name4 <- "2R_9499163;Wnt2"
name5 <- "3R_31651870;kek6"
name6 <- "2L_18962002;Nak"
name7 <- "3R_7682579;Nlg1"
name8 <- "X_5169837;rg"
name9 <- "X_18800322;S6KL"
name10 <- "3R_29036818;Cpsf100"
name11 <- "3R_29446019;Ptp99A"
name12 <- "3L_13501425;Tgi"

plot1 <- plot_allele_f(name1)
plot2 <- plot_allele_f(name2)
plot3 <- plot_allele_f(name3)
plot4 <- plot_allele_f(name4)
plot5 <- plot_allele_f(name5)
plot6 <- plot_allele_f(name6)
plot7 <- plot_allele_f(name7)
plot8 <- plot_allele_f(name8)
plot9 <- plot_allele_f(name9, is_corner = TRUE)
plot10 <- plot_allele_f(name10)
plot11 <- plot_allele_f(name11)
plot12 <- plot_allele_f(name12, is_last = TRUE)

full_plots <- plot_grid(plot1, plot2, plot3, plot4,
                        plot5, plot6, plot7, plot8,
                        plot9, plot10, plot11, plot12,
                        nrow = 3,
                        ncol = 4)
full_plots
ggsave("/Users/arteen/Desktop/School/pop_gen_revision/figure1_mean_se.png",
       width = 15, height = 7, units = "in", dpi = 600)



