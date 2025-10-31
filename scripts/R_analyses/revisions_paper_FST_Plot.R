library(ggplot2)
library(cowplot)
library(stringr)

fst_in <- read.csv(file = "/Users/arteen/Desktop/pop_gen_out/fst/sync_subsetted/fst.csv") # Mpileup



source("/Users/arteen/Desktop/School/Scripts/SociabilityDNA/R/fst_functions.R")

comparisons <- c("Ancestor.Control", "Ancestor.Low", "Ancestor.High", "Control.Low", "Control.High", "Low.High")


## genes to point to / change shape for
genes_tested_vec <- c("FBgn0027873_3R_29034573_29037062",
                      "FBgn0039862_3R_31622136_31664646",
                      "FBgn0266098_X_5085759_5254864",
                      "FBgn0263930_3L_8827505_8891192",
                      "FBgn0051146_3R_7658247_7694874",
                      "FBgn0004369_3R_29377646_29487131",
                      "FBgn0283473_X_18788464_18806057",
                      "FBgn0052264_3L_3714826_3806459",
                      "FBgn0036373_3L_13486797_13502492",
                      "FBgn0052600_X_14322226_14462351",
                      "FBgn0262975_3R_23185580_23226711",
                      "FBgn0003435_2R_19517864_19631508",
                      "FBgn0004360_2R_9494166_9503188",
                      "FBgn0020251_3L_6495680_6550104",
                      "FBgn0262870_3L_4629604_4687317",
                      "FBgn0038201_3R_14099832_14103122",
                      "FBgn0038200_3R_14097293_14098080",
                      "FBgn0046113_2L_914086_917805",
                      "FBgn0001085_3L_14274343_14368639")


gene_test_df <- do.call(rbind.data.frame, str_split(genes_tested_vec, "_"))
colnames(gene_test_df) <- c("FBgnID", "chrom", "start", "end")

#comparison <- comparisons[6] # to test
plot_fst <- function(comparison) {
  
  init_fst <- fst_in[, c("chrom", "start", "end", comparison)]
  
  sex_chrom_fst <- subset(init_fst, init_fst$chrom == "X")
  auto_chrom_fst <- subset(init_fst, !(init_fst$chrom == "X"))
  
  names(sex_chrom_fst)[4] <-  "mean"
  sex_chrom_fst <- na.omit(sex_chrom_fst)
  sex_chrom_fst$highlight <- ifelse(sex_chrom_fst$mean >= quantile(sex_chrom_fst$mean, 0.95), "high", "normal")
  names(auto_chrom_fst)[4] <-  "mean"
  auto_chrom_fst <- na.omit(auto_chrom_fst)
  auto_chrom_fst$highlight <- ifelse(auto_chrom_fst$mean >= quantile(auto_chrom_fst$mean, 0.95), "high", "normal")
  
  subsetted_fst <- rbind(auto_chrom_fst, sex_chrom_fst)
  
  
  subsetted_fst <- chrNumbering(subsetted_fst)
  
  
  numbering_chr <- middleChr(subsetted_fst)
  
  
  
  chrom_colors <- c('X' = 'grey46', '2L' = 'black', '2R' = 'grey46', '3L' = 'black', '3R' = 'grey46', '4' = 'black')
  subsetted_fst$point_color <- ifelse(subsetted_fst$highlight == "high", "red", chrom_colors[subsetted_fst$chrom])
  
  subsetted_fst$shape <- "none"
  
  full_df <- data.frame()
  ### QUESTION FOR IAN - WHAT IF ITS BETWEEN TWO WINDOWS?
  
  
  # Multiple window for gene
  for (i in 1:nrow(gene_test_df)) {
    #i = 1
    chrom <- gene_test_df[i,]$chrom
    start <- as.integer(gene_test_df[i,]$start)
    end <- as.integer(gene_test_df[i,]$end)
    distance <- end - start
    
    starting_window <- as.integer((floor(start/5000) * 5000 + 1))
    ending_window <- as.integer(ceiling(end/5000) * 5000)
    
    # not a single window
    if (ending_window - starting_window != 4999) {
      num_windows <- ((ending_window - starting_window) + 1) /5000

      for (j in 1:num_windows) {
        constant <- j - 1
        
        end_window_to_grab <- starting_window + (4999*j) + (1*constant)
        if (nrow(subsetted_fst[which(subsetted_fst$end == end_window_to_grab & subsetted_fst$chrom == chrom),]) == 0){
          next
        }
        tmp_df <- subsetted_fst[which(subsetted_fst$end == end_window_to_grab & subsetted_fst$chrom == chrom),]
        full_df <- rbind(full_df, tmp_df)
        subsetted_fst[which(subsetted_fst$end == end_window_to_grab & subsetted_fst$chrom == chrom),]$shape <- "sig"
      }
    }else { #single window so just grab it
      tmp_df <- subsetted_fst[which(subsetted_fst$end == ending_window & subsetted_fst$chrom == chrom),]
      full_df <- rbind(full_df, tmp_df)
      subsetted_fst[which(subsetted_fst$end == end_window_to_grab & subsetted_fst$chrom == chrom),]$shape <- "sig"
    }
  }
  
  # single window per gene
  # for (i in 1:nrow(gene_test_df)) {
  #   #i = 1
  #   chrom_match <- gene_test_df[i,]$chrom
  #   start_match <- as.integer(gene_test_df[i,]$start)
  #   end_match <- as.integer(gene_test_df[i,]$end)
  #   find_window_subset <- subset(subsetted_fst, subsetted_fst$chrom == chrom_match) 
  #   find_window_subset <- subset(find_window_subset, find_window_subset$start <= start_match) 
  #   find_window_subset <- subset(find_window_subset, find_window_subset$end >= start_match) 
  #    
  #   full_df <- rbind(full_df, find_window_subset)
  # }
 
  ggplot(data = subsetted_fst, aes(x = number, y = mean)) + 
    geom_point(aes(color = point_color,
               shape = shape, size = ifelse(shape == "sig", 2, 0.65)), show.legend = F, alpha = 0.3) + 
    theme(panel.background = element_blank()) +
    xlab("Chromosome") +
    ylab(bquote(F[ST]))+
    scale_x_discrete(limits = c(numbering_chr), labels = c("X", "2L", "2R", '3L', '3R', '4')) +
    geom_smooth(aes(group = as.factor(chrom)), colour = "red", size = 0.5) +
    scale_color_identity() +
    scale_shape_manual(values = c("none" = 16, "sig" = 17)) +
    scale_size_identity() +  
    
    
    theme(text = element_text(size = 15),
          axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12))
   
}


Ancestor.Control.plot <- plot_fst(comparison = "Ancestor.Control")
Ancestor.Low.plot <- plot_fst(comparison = "Ancestor.Low")
Ancestor.High.plot <- plot_fst(comparison = "Ancestor.High")
Control.Low.plot <- plot_fst(comparison = "Control.Low")
Control.High.plot <- plot_fst(comparison = "Control.High")
Low.High.plot <- plot_fst(comparison = "Low.High")
ggsave("/Users/arteen/Desktop/School/pop_gen_revision/updated_FST_manhattan_plot_low_vs_high.png", width =7, height = 4, dpi = 600)

