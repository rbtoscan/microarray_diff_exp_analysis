setwd("/path/to/workdir")
# download data
# install packages if needed
#install.packages(c("limma", "Biobase"))
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("GEOquery")

# load packages
library(limma)
library(Biobase)
library(vctrs)
library(GEOquery)
library(readr)
library(data.table)
library(dplyr)
library(pheatmap)
library(ggrepel)
library(ggplot2)

# link to dataset to be used
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32368

# define gse code to be used
gse_code <- "GSE32368"

# fetch data using GEOquery API
gse <- getGEO(gse_code, GSEMatrix = TRUE)

# download files locally for inspection - optional. only for manual inspection
#filePaths = getGEOSuppFiles(gse_code)

# get the expression data matrix
expression <- exprs(gse[[1]])
# get metadata
metadata <- pData(gse[[1]])
# get annotation
annotation <- fData(gse[[1]])

# data overview
#summary(expression)
boxplot(expression,outline=FALSE)

# subset expression, metadata and annotation based on the task interest
metadata_I3C <- subset(metadata, substr(description, nchar(source_name_ch1) - 2, nchar(source_name_ch1)) == "I3C")
metadata_DMSO <- subset(metadata, substr(description, nchar(source_name_ch1) - 3, nchar(source_name_ch1)) == "DMSO")
metadata <- rbind(metadata_I3C,metadata_DMSO)
expression <-  expression[ , c(rownames(metadata))]
rm(metadata_I3C,metadata_DMSO)

# add extra column tagging the relevant classes
metadata <- metadata %>%
  mutate(tag = ifelse(grepl("_DMSO$", source_name_ch1), "DMSO", "I3C"))

# data overview - divide this by color
summary(expression)
colors <- c("DMSO" = "blue", "I3C" = "green")
metadata$tag <- factor(metadata$tag, levels = c("DMSO", "I3C"), labels = c("DMSO", "I3C"))
boxplot(expression, outline = FALSE, col = colors[metadata$tag])
legend("topright", legend = levels(metadata$tag), fill = colors, title = "Sample Type")


# fetch sample identifiers and rename them accordingly
sampleInfo <- select(metadata, tag,source_name_ch1)
sampleInfo <- rename(sampleInfo,treatment=tag,name = source_name_ch1)

# set up the design matrix
design <- model.matrix(~0 + tag, metadata)
colnames(design) <- c('DMSO','I3C')

# perform normalization using limma functions
normExpression <- normalizeBetweenArrays(log2(expression), method = "quantile")
fit <- lmFit(normExpression, design)
fit <- eBayes(fit)

# extract differentially expressed genes
de_genes <- topTable(fit, coef = "I3C", number = Inf)

# filter genes based on adjusted p-value (e.g., p value < 0.05)
de_genes <- de_genes[de_genes$adj.P.Val < 0.05, ]

# view the top differentially expressed genes
head(de_genes)

# ensure row names in de_genes match those in normExpression
common_genes <- intersect(rownames(de_genes), rownames(normExpression))

# extract normalized expression data for common genes
de_expression <- normExpression[common_genes, ]

# check for and handle missing values
de_expression[is.na(de_expression)] <- 0  

# create a heatmap - removed gene names for the purpose of a clearer plot 
heatmap_plot <- pheatmap(de_expression, annotation_col = sampleInfo, cluster_cols = FALSE,show_rownames = FALSE)

# Display the heatmap plot
print(heatmap_plot)



