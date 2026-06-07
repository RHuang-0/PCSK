# ==============================================================================
# pathway_pcsk.R
# GO/Pathway Enrichment Dot Plot Visualization
#
# Purpose: Generate dot plots from GO enrichment analysis results.
#          Each section corresponds to a different gene set / comparison.
#          Dot size = number of genes; Dot color = pathway category.
# Note: This is one example - all the others are generated through similar code
# ==============================================================================

library(ggplot2)
library(openxlsx)
library(dplyr)
library(stringr)



data <- read.xlsx("12down.xlsx")
data$`Adjusted.p-value` <- as.numeric(data$`Adjusted.p-value`)

data <- data %>%
  mutate(
    plotP = -log10(`Adjusted.p-value`),
    Gene.number = ifelse(is.na(Genes) | Genes == "", 0, str_count(Genes, ";") + 1)
  ) %>%
  arrange(desc(plotP))
data <- data[order(data$plotP, decreasing = TRUE), ]

ggplot(data, aes(x = plotP, y = reorder(Pathway, -plotP))) +
  geom_point(aes(size = Gene.number, color = Category)) +
  scale_size_continuous(range = c(3, 10), limits = c(0, 80)) +
  scale_color_manual(values = c(
    "Protein modification" = '#ba52c0',
    "Extracellular matrix organization" = '#499331',
    "Cell adhesion" = '#44549c',
    "Protein metabolism" = '#6e701e',
    "Signal transduction" = "#e2769a"
  )) +
  labs(
    x = "-log10(Adjusted P-value)",
    y = "Term",
    size = "Gene number",
    color = "Category"
  ) +
  theme_classic() +
  coord_fixed(ratio = 1.5)

