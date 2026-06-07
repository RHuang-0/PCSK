# ==============================================================================
# pcsk.R
# RNA-seq Analysis: PCA, DESeq2 DEA, Heatmap, Violin/Boxplot, Proteomics Integration
#
# Experiment: Time-course (0h, 4h, 8h, 12h) +/- PCI treatment (Ctrl 0/1)
# Species: Mouse 


# ==============================================================================


# ==============================================================================
# SECTION 1: Package Installation & Library Loading
# Run install lines once; comment out after first use
# ==============================================================================

# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# BiocManager::install('PCAtools')
# BiocManager::install("Biobase")
# BiocManager::install('DESeq2')
# install.packages("pheatmap")
# BiocManager::install("ComplexHeatmap")
# install.packages("VennDiagram")

library(readxl)
library(PCAtools)
library(pheatmap)
library(Biobase)
library(dplyr)
library(DESeq2)
library(tibble)
library(ggplot2)
library(ggrepel)
library(tidyr)


# ==============================================================================
# SECTION 2: Data Loading & Preprocessing
# Loads count matrix and metadata; filters to protein-coding genes
# ==============================================================================

meta <- read.csv("Book1.csv", header = TRUE)
haha <- meta[[1]]
rownames(meta) <- haha
meta <- meta[, -1]
rownames(meta) <- haha  


df <- read_excel('Duncan10_Counts (1).xlsx')
df <- df[df$Biotype == "protein_coding", ]
df <- df %>% select(-ENSEMBL)
df <- df %>% select(-Location)
df <- df %>% select(-Strand)
df <- df %>% select(-Biotype)
df <- df %>% filter(!is.na(Gene))
df <- df %>% distinct(Gene, .keep_all = TRUE)
ha <- df$Gene   
rownames(df) <- ha
df <- df %>% select(-Gene)

df <- df[, sort(colnames(df))]
meta <- meta[sort(rownames(meta)), ]


# ==============================================================================
# SECTION 3: PCA
# ==============================================================================

pca_result <- pca(df, metadata = meta)
biplot(pca_result, lab = NULL, colby = "Time", legendPosition = "right")



# ==============================================================================
# SECTION 4: Group Construction (Time x Condition)
# ==============================================================================

group_colors <- c(
  "0h"     = "#E31A1C",  
  "4h CTL" = "#FF7F00",  
  "4h PCI" = "#FDBF6F", 
  "8h CTL" = "#33A02C",  
  "8h PCI" = "#B2DF8A",  
  "12h CTL"= "#6A3D9A", 
  "12h PCI"= "#CAB2D6"   
)

meta$Condition <- ifelse(meta$Ctrl == 0, "CTL", "PCI")
meta$Group <- paste0(meta$Time, "h ", meta$Condition)
meta$Group[meta$Time == 0] <- "0h"

meta$Group <- factor(meta$Group,
                     levels = c("0h", "4h CTL", "4h PCI", "8h CTL", "8h PCI", "12h CTL", "12h PCI"))
pca_result <- pca(df, metadata = meta)

biplot(pca_result,
        lab = NULL,
        colby = "Group",
        colkey = group_colors,
        ellipse = TRUE,
        legendPosition = 'right')

biplot(pca_result,
        lab = NULL,
        colby = "Group",
        colkey = group_colors,
       legendPosition = "right",
        pointSize = 4,
        title = "PCA Analysis: Grouped by Time & Treatment",
        subtitle = "Colors: 0h(Red), 4h(Orange), 8h(Green), 12h(Purple)")


# ==============================================================================
# SECTION 5: ComplexHeatmap — Gene-of-Interest Heatmaps
# Draws Z-score scaled heatmaps for candidate gene sets
# Multiple gene lists below; uncomment the one you want to use
# ==============================================================================

library(ComplexHeatmap)
library(circlize)
library(tidyverse)


# genes_of_interest <- c(
#   "Itgb5", "Lamb3", "Lama1", "Lama4", "Tnc",
#   "Dmp1", "Lamc2", "Thbs1", "Col2a1", "Ibsp",
#   "Itga11", "Col6a3", "Itga5", "Vav3", "Pdgfra",
#   "Jun", "Ccnd2", "Flnc", "Pak3", "Raf1",
#   "Fgf2", "Lama5", "Itgb3", "Itga2", "Itga1", "Col1a1"
# )

# genes_of_interest <- c(
#   "Ret", "Csf1r", "Igsf8", "Grn", "Cemip",
#   "Sema7a", "Mctp1", "Plxnd1", "Cited2", "Sema3b",
#   "Sema3g", "Arid4a", "Rora", "Ptn", "Vsir",
#   "Egfr", "Robo1", "Ppp3ca", "Clec7a", "Plau",
#   "Pdgfd", "Gna12", "Sgk1", "Fermt2", "Vclcavin1",
#   "Scarb1"
# )

# genes_of_interest <- c(
#   "Ramp2", "Calcoco1", "Ramp3", "Ddx3x", "Tfrc",
#   "Cited2", "Irs1", "Hfe", "Nts", "Areg",
#   "Vsir", "Robo1", "C1qtnf1", "Hmgn5", "Clec7a",
#   "Myc", "Ubr5", "Ripk1", "Slc38a2", "Egr1",
#   "Mef2c", "Insr", "Fn1", "F3", "Tgfbr1",
#   "Tgfbr3"
# )

# genes_of_interest <- c("Dc4", "Sema7a", "Itga5", "Tsg6", "Plau", "Plaur", "Tgfbr1", "Fst", "Bmpr2", "Dab2")

# genes_of_interest <- c('Igfbp5', 'Cemip', 'Sfrp2', 'Bmp2', 'Vsir','Clec7', 'Ptk2b', 'Pdgfd',
#                        'F3', 'Plau', 'Rora', 'Cyp1b1', 'Sgk1', 'Robo1', 'Tgfbr1')

# genes_of_interest <- c('Col2a1', 'Lamb3', 'Fgf2', 'Ibsp', 'Dmp1', 'Col6a3', 'Flnc', 'Lama4', 'Pdgfra',
#                         'Lamc2', 'Itga5', 'Vav3', 'Itga11', 'Lama1', 'Itgb5')

# genes_of_interest <- c('Elovl3', 'Adh7', 'Cemip', 'Cyp11a1', 'Itpka', 'Hdc', 'Bmp2', 'Pik3c2g',
#                        'Npr1', 'Cd28', 'Vsir', 'Clec7a', 'Coq8a', 'A2m', 'Gm2a')

# genes_of_interest <- c('Tnfaip6','Synd4')


genes_present <- intersect(genes_of_interest, rownames(df))
if (length(genes_present) < length(genes_of_interest)) {
   warning("Some genes not found in df and will be skipped.")
}

mat_subset <- df[genes_present, ]
mat_numeric <- as.matrix(mat_subset)
mat_scaled <- t(scale(t(mat_numeric)))  # row-wise Z-score
sample_groups <- meta$Group
top_anno <- HeatmapAnnotation(
   Group = sample_groups,
   col = list(Group = group_colors),
   show_annotation_name = TRUE
)

col_fun <- colorRamp2(c(-2, 0, 2), c("#313695", "white", "#a50026"))

Heatmap(mat_scaled,
         name = "Z-score",
         col = col_fun,
         cluster_rows = TRUE,
         cluster_columns = TRUE,
         top_annotation = top_anno,
         show_row_names = TRUE,
         row_names_side = "right",
         row_labels = rownames(mat_scaled),
         row_names_gp = gpar(fontsize = 10),
         show_column_names = TRUE,
         column_names_rot = 90,
         border = TRUE)


target_groups <- c("8h CTL", "8h PCI", "12h CTL", "12h PCI")
keep_cols <- sample_groups %in% target_groups
top_anno_sub <- HeatmapAnnotation(
   Group = as.character(sample_groups[keep_cols]),
   col = list(Group = group_colors[target_groups]),
   show_annotation_name = TRUE
)
Heatmap(mat_scaled[, keep_cols],
         name = "Z-score",
         col = col_fun,
         cluster_rows = TRUE,
         cluster_columns = TRUE,
         top_annotation = top_anno_sub,
         show_row_names = TRUE,
         row_names_side = "right",
         row_labels = rownames(mat_scaled),
         row_names_gp = gpar(fontsize = 10),
         show_column_names = TRUE,
         column_names_rot = 90,
         border = TRUE)


# ==============================================================================
# SECTION 6: Sample Correlation Heatmap (pheatmap)
# Pearson correlation across all samples; annotated by Ctrl and Time
# ==============================================================================

tpm_set <- ExpressionSet(assayData = data.matrix(df), phenoData = AnnotatedDataFrame(meta))
tpm_cor <- cor(exprs(tpm_set))
pData(tpm_set)$Ctrl <- as.factor(pData(tpm_set)$Ctrl)
pData(tpm_set)$Time <- as.factor(pData(tpm_set)$Time)
meta.stage1 <- as.data.frame(cbind(pData(tpm_set)$Ctrl, pData(tpm_set)$Time))
colnames(meta.stage1) <- c("Ctrl", "Time")
rownames(meta.stage1) <- row.names(pData(tpm_set))
heatmap_tpm_adjusted <- pheatmap(tpm_cor, fontsize = 10, height = 20, width = 20,
                                  annotation_row = meta.stage1)
heatmap_tpm_adjusted


# ==============================================================================
# SECTION 7: DESeq2 — Time-Point-Stratified DEA
# Compares PCI (Ctrl=1) vs CTL (Ctrl=0) at each time point
# ==============================================================================

meta$Ctrl <- factor(meta$Ctrl, levels = c(0, 1))


perform_DEA <- function(time_point) {
   meta_subset <- meta %>% filter(Time == time_point)
   df_subset <- df[, rownames(meta_subset)]
   rownames(df_subset) <- ha
   dds <- DESeqDataSetFromMatrix(countData = df_subset,
                                 colData = meta_subset,
                                 design = ~ Ctrl)
   dds <- DESeq(dds)
  res <- results(dds, contrast = c("Ctrl", "1", "0"))
   res_df <- as.data.frame(res) %>%
     rownames_to_column("Gene") %>%
     arrange(padj)
   res_sig <- res_df %>% filter(padj < 0.05)
   write.csv(res_df, paste0("0419DEG_results_time_", time_point, ".csv"), row.names = FALSE)
   write.csv(res_sig, paste0("0419DEG_significant_time_", time_point, ".csv"), row.names = FALSE)
   return(res_sig)
}
deg_time_4  <- perform_DEA(4)
deg_time_8  <- perform_DEA(8)
deg_time_12 <- perform_DEA(12)



# ==============================================================================
# SECTION 8: Volcano Plots
# Visualizes DEA results with log2FC on X and -log10(padj) on Y
# Thresholds: padj < 0.05, |log2FC| > 1
# ==============================================================================

draw_volcano_plot <- function(file_path, time_point) {
   df <- read.csv(file_path)
   df <- df %>% filter(!is.na(padj))
   log2FC_threshold <- 1
   pvalue_threshold <- 0.05
   df <- df %>%
     mutate(Significance = case_when(
       padj < pvalue_threshold & log2FoldChange > log2FC_threshold  ~ "Upregulated",
       padj < pvalue_threshold & log2FoldChange < -log2FC_threshold ~ "Downregulated",
       TRUE ~ "Not Significant"
     ))
   df$Significance <- factor(df$Significance,
                             levels = c("Upregulated", "Downregulated", "Not Significant"))
   p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
     geom_point(alpha = 0.6, size = 2) +
     scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue",
                                   "Not Significant" = "gray")) +
     theme_minimal() +
     geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dashed") +
     geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed") +
     labs(title = paste("Volcano Plot - Time", time_point),
          x = "Log2 Fold Change",
          y = "-Log10 Adjusted P-Value",
          color = "Significance") +
     theme(legend.position = "top")
     # Optional gene labels for top hits:
     # geom_text_repel(data = df %>% filter(padj < 0.01 & abs(log2FoldChange) > 3.5),
     #                 aes(label = Gene), size = 3, max.overlaps = 15)
   return(p)
}
p4  <- draw_volcano_plot("0419DEG_results_time_4.csv", 4)
p8  <- draw_volcano_plot("0419DEG_results_time_8.csv", 8)
p12 <- draw_volcano_plot("0419DEG_results_time_12.csv", 12)
print(p4); print(p8); print(p12)


# ==============================================================================
# SECTION 9: Violin & Boxplot Functions
# ==============================================================================

plot_violin <- function(df, meta, gene_list) {
   df_long <- as.data.frame(t(df)) %>%
     rownames_to_column("Sample") %>%
     pivot_longer(-Sample, names_to = "Gene", values_to = "Expression") %>%
     left_join(meta %>% rownames_to_column("Sample"), by = "Sample")
   df_long <- df_long %>% filter(Gene %in% gene_list)
   df_long <- df_long %>% mutate(Group = paste0("Ctrl_", Ctrl, "_Time_", Time))
   desired_order <- c("Ctrl_0_Time_0", "Ctrl_0_Time_4", "Ctrl_1_Time_4",
                      "Ctrl_0_Time_8", "Ctrl_1_Time_8", "Ctrl_0_Time_12", "Ctrl_1_Time_12")
   df_long$Group <- factor(df_long$Group, levels = desired_order)
   p <- ggplot(df_long, aes(x = Group, y = Expression, fill = as.factor(Ctrl))) +
     geom_violin() +
     geom_jitter(width = 0.2, size = 0.5, alpha = 0.6) +
     facet_wrap(~ Gene, scales = "free_y") +
     labs(x = "Group (Ctrl + Time)", y = "Expression", title = "Violin Plot of Gene Expression") +
     theme_minimal() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
     scale_fill_manual(values = c("0" = "blue", "1" = "red"))
   # ggsave(filename = "violin_plot_ctrl_time.png", plot = p, width = 10, height = 6)
   return(p)
}


library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpubr)

plot_boxplot_batched <- function(df, meta, gene_list, genes_per_plot = 3) {

  library(rstatix)
  df_long <- as.data.frame(t(df)) %>%
    rownames_to_column("Sample") %>%
    pivot_longer(-Sample, names_to = "Gene", values_to = "Expression") %>%
    left_join(meta %>% rownames_to_column("Sample"), by = "Sample") %>%
    filter(Gene %in% gene_list) %>%
    mutate(
      Group = paste0("Ctrl_", Ctrl, "_Time_", Time),
      Time  = as.character(Time),
      Ctrl  = as.factor(Ctrl)
    )

  desired_order <- c(
    "Ctrl_0_Time_0",
    "Ctrl_0_Time_4", "Ctrl_1_Time_4",
    "Ctrl_0_Time_8", "Ctrl_1_Time_8",
    "Ctrl_0_Time_12","Ctrl_1_Time_12"
  )
  pretty_labels <- c(
    "0h",
    "4h CTL", "4h PCI",
    "8h CTL", "8h PCI",
    "12h CTL", "12h PCI"
  )

  df_long <- df_long %>% filter(Group %in% desired_order)
  df_long$Group <- factor(df_long$Group, levels = desired_order)

  gene_chunks <- split(gene_list, ceiling(seq_along(gene_list) / genes_per_plot))
  plot_list <- list()

  comps <- list(
    c("Ctrl_0_Time_4",  "Ctrl_1_Time_4"),
    c("Ctrl_0_Time_8",  "Ctrl_1_Time_8"),
    c("Ctrl_0_Time_12", "Ctrl_1_Time_12")
  )

  for (i in seq_along(gene_chunks)) {

    current_genes <- gene_chunks[[i]]
    df_subset <- df_long %>% filter(Gene %in% current_genes)

    stat_df <- df_subset %>%
      filter(Time %in% c("4","8","12")) %>%
      group_by(Gene, Time) %>%
      wilcox_test(Expression ~ Ctrl, exact = FALSE) %>%
      ungroup() %>%
      mutate(
        group1 = paste0("Ctrl_0_Time_", Time),
        group2 = paste0("Ctrl_1_Time_", Time),
        label  = ifelse(p < 0.1, "*", NA_character_)
      ) %>%
      filter(!is.na(label))

    if (nrow(stat_df) > 0) {
      ymax_df <- df_subset %>%
        group_by(Gene) %>%
        summarise(ymax = max(Expression, na.rm = TRUE))
      stat_df <- stat_df %>%
        left_join(ymax_df, by = "Gene") %>%
        group_by(Gene) %>%
        mutate(y.position = ymax * (1.08 + 0.08 * row_number())) %>%
        ungroup()
    }

    p <- ggplot(df_subset, aes(x = Group, y = Expression, fill = Ctrl)) +
      geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.7, color = "black") +
      geom_jitter(width = 0.2, size = 1, alpha = 0.5, color = "black") +
      facet_wrap(~ Gene, scales = "free_y", ncol = genes_per_plot) +
      scale_x_discrete(labels = pretty_labels) +
      scale_fill_manual(values = c("0" = "white", "1" = "black"),
                        labels = c("0" = "Control", "1" = "PCI")) +
      labs(x = NULL, y = "Normalized Expression", fill = "Condition",
           title = paste("Batch", i)) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
      coord_cartesian(clip = "off") +
      theme_classic(base_size = 14) +
      theme(
        plot.margin = margin(8, 8, 8, 14),
        axis.text.x  = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y  = element_text(color = "black"),
        strip.background = element_rect(fill = "#f0f0f0", color = NA),
        strip.text   = element_text(face = "bold", size = 12),
        legend.position = "top",
        legend.justification = "right"
      )

    if (nrow(stat_df) > 0) {
      p <- p +
        stat_pvalue_manual(
          stat_df,
          label = "label",
          xmin  = "group1",
          xmax  = "group2",
          y.position = "y.position",
          tip.length = 0.01,
          size = 5,
          bracket.size = 0.4
        )
    }

    plot_list[[i]] <- p
  }

  return(plot_list)
}


# ==============================================================================
# SECTION 10: Gene Lists for Boxplot / Violin Plotting
# ==============================================================================

# list1 <- c("Igfbp5", "Sfrp2", "Ptk2b", "Pdgfd", "Plau", "Sgk1")
# list2 <- c("Col2a1", "Lamb3", "Fgf2", "Ibsp", "Dmp1", "Col6a3", "Flnc",
#             "Lama4", "Pdgfra", "Lamc2", "Itga5", "Vav3", "Itga11", "Lama1", "Itgb5")
# list3 <- c("Elov3", "Adh7", "Cyp11a1", "Itpka", "Hdc", "Pik3c2g", "Npr1",
#             "Cd28", "Coq8a", "A2m", "Gm2a")
# list4 <- c("Cemip", "Bmp2", "Vsir", "Clec7a", "F3", "Rora", "Cyp1b1",
#             "Robo1", "Tgfbr1")
# list5 <- c("Has2", "Vcan","Ptx3")
# list  <- c("Col6a3","Itgb5","Lama1","Has2","Ptx3","Vcan",
#             "Cyp11a1","Gm2a","Igfbp5","Ptk2b","Cyp1b1","Cemip")

# gene_list <- c("Furin", "Cemip", "Elmo1","Acvrl1","Gjb3","Adamts9")
# gene_list <- c('Gm21663','Dlx2')
# gene_list <- c('Tmem229a','Serpina12')
# gene_list <- c('Has2','Tnfaip6','Btc','Areg','Ereg','Ptx3','Vcan','Bik',
#                'Ptgs2','Cd44','Hmmr','Itih1','Itih2','Itih3')
# gene_list <- c('Pcsk4','Pcsk3','Pcsk6','Pcsk9','Pcsk1','Pcsk2','Pcsk7','Pcsk8')
# gene_list <- c('Fzd1','Rarg','Ddx3x','Fxd4','Lrrk2','Cav1','Tcf7','Dixdc1',
#                'Lrp6','Nr4a2','Sfrp4','Frat2','Sfrp2','Gpc4','Fermt2')
# gene_list <- c('Jun','Cited1','Itgb5','Hpgd','Smad9','Ptprk','Smad7','Dab2',
#                'Bambi','Id1','Ankrd1','Gcnt2','Sox9','Sox5','Grem1','Xbp1',
#                'Pdzd2','Nog','Pparg','Ldlrad4','Smad6','Id4','Id3','Nbl1','Thbs1')
# gene_list <- c('Bmp2','Bmp1','Bmpr1a','Adamts4','Adamts9','Nts',
#                'Igfbp5','Igf1','Itga5','Itgb5')

plot_boxplot_batched(df, meta, gene_list)
plot_violin(df, meta, gene_list)


# ==============================================================================
# SECTION 11: Sample Correlation Analysis
# Computes Pearson correlation on top-variable genes; outputs heatmap + boxplot
# The boxplot shows which group 12h Control is most similar to
# ==============================================================================

library(tidyverse)
library(pheatmap)

run_correlation_analysis <- function(df, meta, num_variable_genes = 1000) {
   meta_clean <- meta %>%
     rownames_to_column("Sample") %>%
     mutate(Condition = ifelse(Ctrl == 1, "PCI", "Control"),
            Group = paste(Time, "h_", Condition, sep = "")) %>%
     column_to_rownames("Sample")
   gene_vars <- apply(df, 1, var)
   top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:num_variable_genes]
   df_filtered <- df[top_genes, ]
   cor_mat <- cor(df_filtered, method = "pearson")
   ann_colors = list(
     Condition = c(Control = "#377EB8", PCI = "#E41A1C"),
     Time = c("0" = "#f0f0f0", "4" = "#bdbdbd", "8" = "#636363", "12" = "#252525")
   )
   pheatmap(cor_mat,
            annotation_col = meta_clean[, c("Condition", "Time")],
            annotation_row = meta_clean[, c("Condition", "Time")],
            show_rownames = FALSE, show_colnames = FALSE,
            main = "Global Sample-to-Sample Correlation (Top Var Genes)",
            border_color = NA)
}
results <- run_correlation_analysis(df, meta)


# ==============================================================================
# SECTION 12: Cross-Group Comparison — 8h PCI vs 12h CTL
# Loads DEG list from 8h, subsets to those two groups, runs t-test per gene
# ==============================================================================

library(tidyverse)
library(readr)

gene_file <- "0419DEG_significant_time_8.xls"
target_genes_df <- read_excel(gene_file)
gene_list <- target_genes_df$Gene
valid_genes <- intersect(gene_list, rownames(df))
if (length(valid_genes) < length(gene_list)) {
   message(paste("Note:", length(gene_list) - length(valid_genes),
                 "genes were not found in 'df' and will be skipped."))
}
samps_8h_PCI  <- meta %>% filter(Time == 8 & Ctrl == 1) %>% rownames()
samps_12h_Ctrl <- meta %>% filter(Time == 12 & Ctrl == 0) %>% rownames()
sub_df <- df[valid_genes, c(samps_8h_PCI, samps_12h_Ctrl)]
results <- data.frame(Gene = valid_genes, stringsAsFactors = FALSE) %>%
   mutate(
     Mean_8h_PCI   = rowMeans(sub_df[, samps_8h_PCI], na.rm = TRUE),
     Mean_12h_Ctrl = rowMeans(sub_df[, samps_12h_Ctrl], na.rm = TRUE),
     Log2FC_8hPCI_vs_12hCtrl = log2((Mean_8h_PCI + 1) / (Mean_12h_Ctrl + 1)),
     P_Value = apply(sub_df, 1, function(x) {
       g1 <- x[samps_8h_PCI]; g2 <- x[samps_12h_Ctrl]
       if (var(g1) == 0 && var(g2) == 0) return(1)
       res <- try(t.test(g1, g2), silent = TRUE)
       if (inherits(res, "try-error")) return(NA) else return(res$p.value)
     })
   )
results$FDR <- p.adjust(results$P_Value, method = "BH")
results <- results %>%
   mutate(Interpretation = case_when(
     P_Value < 0.05 & Log2FC_8hPCI_vs_12hCtrl >  0.5 ~ "Higher in 8h PCI",
     P_Value < 0.05 & Log2FC_8hPCI_vs_12hCtrl < -0.5 ~ "Higher in 12h Ctrl",
     P_Value >= 0.05 ~ "Similar / No Sig Diff",
     TRUE ~ "Ambiguous"
   ))
write.csv(results, "Comparison_8hPCI_vs_12hCtrl_Results.csv", row.names = FALSE)




# ==============================================================================
# SECTION 13: Proteomics + RNA-seq Integration
# Loads proteomics candidates, integrates with DESeq2 results from each time point,
# then visualizes as a ComplexHeatmap (Log2FC, black = not significant)
# and a scatter plot (RNA 8h vs Proteomics Log2FC)
# ==============================================================================

library(readxl); library(dplyr); library(DESeq2); library(tibble); library(tidyr)
library(ComplexHeatmap); library(circlize); library(Seurat)


prot_data <- read_excel("Candidates_CEK2_150COCs_1to10dilu_directDIA_v2 (1).xlsx", sheet = 2)
target_prot <- prot_data %>%
   dplyr::select("Genes", "AVG Log2 Ratio", "Qvalue") %>%
   dplyr::rename(Gene = "Genes", Prot_Log2FC = "AVG Log2 Ratio", Prot_Qval = "Qvalue") %>%
   dplyr::filter(!is.na(Gene) & Gene != "") %>%
   dplyr::distinct(Gene, .keep_all = TRUE)
target_prot <- as.data.frame(target_prot)

single_cell <- readRDS("10X_dataset_10x_updated_0818.rds")
genes_up   <- target_prot %>% filter(!is.na(Prot_Log2FC) & Prot_Log2FC > 0) %>% pull(Gene)
genes_down <- target_prot %>% filter(!is.na(Prot_Log2FC) & Prot_Log2FC < 0) %>% pull(Gene)
DotPlot(single_cell, features = genes_up) +
   coord_flip() + scale_colour_gradient2(low = "blue", mid = "white", high = "red") +
   labs(title = "Up-regulated Proteins") +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))

perform_DEA <- function(df_subset, meta_subset) {
   dds <- DESeqDataSetFromMatrix(countData = df_subset, colData = meta_subset, design = ~ Ctrl)
   dds <- DESeq(dds)
   res <- results(dds, contrast = c("Ctrl", "1", "0"))
   as.data.frame(res) %>% rownames_to_column("Gene") %>% arrange(padj)
}
time_points <- c("4", "8", "12")
rna_results <- list()
for (tp in time_points) {
   meta_sub <- meta %>% filter(Time == tp)
   df_sub   <- df[, rownames(meta_sub)]
   full_res <- perform_DEA(df_sub, meta_sub)
   res_df <- full_res %>%
     filter(Gene %in% target_prot$Gene) %>%
     select(Gene, log2FoldChange, padj)
   colnames(res_df) <- c("Gene", paste0("RNA_", tp, "h_Log2FC"), paste0("RNA_", tp, "h_Pval"))
   rna_results[[tp]] <- res_df
}

merged_df <- target_prot
for (tp in time_points) merged_df <- left_join(merged_df, rna_results[[tp]], by = "Gene")
merged_df[is.na(merged_df)] <- 1  # set missing p-values to 1 (not significant)
merged_df <- as.data.frame(merged_df)
rownames(merged_df) <- merged_df$Gene
mat_fc   <- as.matrix(merged_df[, c("Prot_Log2FC","RNA_4h_Log2FC","RNA_8h_Log2FC","RNA_12h_Log2FC")])
mat_pval <- as.matrix(merged_df[, c("Prot_Qval","RNA_4h_Pval","RNA_8h_Pval","RNA_12h_Pval")])
colnames(mat_fc) <- colnames(mat_pval) <- c("Proteomics","RNA-4h","RNA-8h","RNA-12h")

col_fun <- colorRamp2(c(-6, 0, 6), c("blue", "white", "red"))
ht <- Heatmap(mat_fc,
               name = "Log2 FC", col = col_fun,
               cluster_columns = FALSE, cluster_rows = TRUE,
               show_row_names = TRUE, row_names_side = "left",
               row_names_gp = gpar(fontsize = 5),
               column_names_gp = gpar(fontsize = 11, fontface = "bold"),
               column_title = "Transcriptome-Proteome Integration",
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if (is.na(mat_pval[i, j]) || mat_pval[i, j] > 0.1)
                   grid.rect(x, y, width, height, gp = gpar(fill = "black", col = "black"))
               })
ht


df_scatter <- merged_df %>%
   filter(!is.na(RNA_8h_Pval) & RNA_8h_Pval <= 0.1) %>%
   filter(!is.na(Prot_Qval) & Prot_Qval <= 0.1)
cor_res <- cor.test(df_scatter$RNA_8h_Log2FC, df_scatter$Prot_Log2FC, method = "pearson")
r_value  <- round(cor_res$estimate, 3)
p_value  <- signif(cor_res$p.value, 3)
cor_label <- paste0("R = ", r_value, "\nP = ", p_value)
p_scatter <- ggplot(df_scatter, aes(x = RNA_8h_Log2FC, y = Prot_Log2FC)) +
   geom_point(color = "black", size = 2, alpha = 0.7) +
   geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
   geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
   geom_smooth(method = "lm", color = "red", fill = "pink", alpha = 0.3) +
   geom_text_repel(aes(label = Gene), size = 3.5, max.overlaps = 20,
                   color = "darkblue", box.padding = 0.5) +
   annotate("text", x = min(df_scatter$RNA_8h_Log2FC, na.rm = TRUE),
            y = max(df_scatter$Prot_Log2FC, na.rm = TRUE),
            label = cor_label, hjust = 0, vjust = 1,
            color = "darkred", fontface = "bold", size = 5) +
   labs(title = "Correlation: Transcriptome (8h) vs Proteome",
        x = "RNA-seq Log2FC (8 hrs)", y = "Proteomics Log2FC") +
   theme_classic(base_size = 14) +
   theme(plot.title = element_text(hjust = 0.5, face = "bold"))
print(p_scatter)
#ggsave("Scatter_RNA8h_vs_Prot.pdf", plot = p_scatter, width = 6, height = 5)


# ==============================================================================
# SECTION 14: GO Biological Process Bar Plots (from Proteomics annotation)
# Explodes GO BP column, counts terms per direction (up/down by Log2FC)
# Plots top 10 terms for each direction
# ==============================================================================

library(tidyr); library(stringr); library(tidytext)
#
target_prot <- prot_data %>%
   dplyr::select("Genes", "AVG Log2 Ratio", "Qvalue", "GO Biological Process") %>%
   dplyr::rename(Gene = "Genes", Prot_Log2FC = "AVG Log2 Ratio",
                 Prot_Qval = "Qvalue", BP = "GO Biological Process") %>%
   dplyr::filter(!is.na(Gene) & Gene != "") %>%
   dplyr::distinct(Gene, .keep_all = TRUE)
target_prot <- as.data.frame(target_prot)


bp_df <- target_prot %>%
   mutate(Direction = ifelse(Prot_Log2FC > 0, "Up-regulated (FC > 0)", "Down-regulated (FC < 0)")) %>%
   separate_rows(BP, sep = ",") %>%
   mutate(BP = str_trim(BP), BP = str_to_title(BP)) %>%
   filter(BP != "")

bp_counts <- bp_df %>% count(Direction, BP, name = "Count")

bp_up   <- bp_counts %>% filter(Direction == "Up-regulated (FC > 0)") %>%
              arrange(desc(Count)) %>% head(10)
bp_down <- bp_counts %>% filter(Direction == "Down-regulated (FC < 0)") %>%
              arrange(desc(Count)) %>% head(10)

plot_go_bar <- function(data, plot_title, bar_color) {
   ggplot(data, aes(x = Count, y = reorder(BP, Count))) +
     geom_bar(stat = "identity", fill = bar_color, color = "black", width = 0.7) +
     theme_bw(base_size = 14) +
     labs(title = plot_title, x = "Frequency (Number of Genes)", y = "") +
     theme(axis.text.y = element_text(color = "black", size = 11),
           plot.title = element_text(face = "bold"))
}
p_up_go   <- plot_go_bar(bp_up,   "Top 10 Pathways in Up-regulated Proteins",   "#de2d26")
p_down_go <- plot_go_bar(bp_down, "Top 10 Pathways in Down-regulated Proteins", "#3182bd")
print(p_up_go)
print(p_down_go)
