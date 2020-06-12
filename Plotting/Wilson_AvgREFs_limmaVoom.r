#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
#source("https://bioconductor.org/biocLite.R")
#BiocInstaller::biocLite("forcats")

library(RColorBrewer)
library(matrixStats)
library(ggplot2)
library(edgeR)
library(DESeq)
library(limma)
library(doParallel)
library(variancePartition)
library(org.Hs.eg.db)
library(clusterProfiler)
library(GOSemSim)
library(biomaRt)
library(UpSetR)
library(VennDiagram)
library(ggrepel)
library(dplyr)
library(stringr)
library(forcats)

# Reading in data.
setwd("~/Dropbox (ASU)/NASONIA_take3/DE_limmaVoom")
counts_pseudoNgir <-
  read.table("wilson_counts_VgirRef.txt",
             header = TRUE,
             sep = "\t")
counts_Nvit <-
  read.table("wilson_counts.txt", header = TRUE, sep = "\t")
DF <- data.frame(counts_Nvit, counts_pseudoNgir)
colnames = c(
  "X014444",
  "X014445",
  "X014446",
  "X014447",
  "X014448",
  "X014449",
  "X014450",
  "X014451",
  "X014452",
  "X014453",
  "X014454",
  "X014455"
)
counts <-
  sapply(colnames, function(x)
    rowMeans(DF [, grep(x, names(DF))]))
genes <- read.table("genes.txt", header = TRUE, sep = "\t")
pheno <- read.table("wilson_pheno.csv", header = TRUE, sep = ",")

removals <- c()
samplesToRemove <-
  c(removals) # update depending on comparison being made

#SAMPLE_LENGTH<-as.numeric(length(samplesToRemove)) # to call later
#half_sample_length <- SAMPLE_LENGTH/2 # half the sample length
#removals<-(names(counts) %in% samplesToRemove[1:SAMPLE_LENGTH]) # for matching names create a value of TRUE or FALSES
#counts <-counts[!removals] # create a new counts file that excludes (Ex) the removals IDs
#pheno <- pheno[! pheno$sampleID %in% samplesToRemove[1:SAMPLE_LENGTH],] # update 1:16 depending on size of samples to remove

viralPalette <- brewer.pal(8, "Set2")
VV_Color <- viralPalette[1]
GG_Color <- viralPalette[2]
VG_Color <- viralPalette[3]
GV_Color <- viralPalette[4]

VV_Shape <- c(15)
GG_Shape <- c(16)
VG_Shape <- c(17)
GV_Shape <- c(18)

# Creating the DGEList object.
dge <- DGEList(counts = counts, genes = genes)
dge$samples$strain <- pheno$Strain
table(dge$samples$strain) # Inspecting the N of samples in each group


# Filtering expression data
fpkm <- rpkm(dge, gene.length = dge$genes$Length)

VV_mean_fpkm <- apply(as.data.frame(fpkm)
                      [(dge$samples$strain == "VV")],
                      1, mean, na.rm = TRUE)
GG_mean_fpkm <- apply(as.data.frame(fpkm)
                      [(dge$samples$strain == "GG")],
                      1, mean, na.rm = TRUE)
VG_mean_fpkm <- apply(as.data.frame(fpkm)
                      [(dge$samples$strain == "VG")],
                      1, mean, na.rm = TRUE)
GV_mean_fpkm <- apply(as.data.frame(fpkm)
                      [(dge$samples$strain == "GV")],
                      1, mean, na.rm = TRUE)

keep <- (VV_mean_fpkm > 0.5 | GG_mean_fpkm > 0.5 |
           VG_mean_fpkm > 0.5 | GV_mean_fpkm > 0.5)

dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
keep <- rowSums(dge$counts > 6) >= 2
dge <- dge[keep, , keep.lib.size = FALSE]
dge <- calcNormFactors(dge, method = "TMM")

# N of genes retained after filtering
dim(dge$genes)

# Creating a design model matrix with the variable of interest
design <- model.matrix( ~ 0 + dge$samples$strain)
colnames(design) <-
  gsub("dge\\$samples\\$strain", "", colnames(design))
head(design)

# Running voom with quality weights. Normalizes expression intensities so that
# the log-ratios have similar distributions across a set of samples.
# To quantile normalize, add normalize.method="quantile"
v <- voomWithQualityWeights(dge, design, plot = TRUE)

# MDS plot of the first two dimensions using top50 genes. The dimensions are in the mds object.
png(filename = "figures/wilson_AvgREFs_mds_1and2.png",
    width = 450,
    height = 450)
mds <-
  plotMDS(
    v,
    top = 50,
    ndim = 10,
    dim.plot = c(1, 2),
    plot = TRUE,
    cex = 2,
    pch = ifelse(
      v$targets$strain %in% c("VV"),
      17,
      ifelse(
        v$targets$strain %in% c("GG"),
        15,
        ifelse(v$targets$strain %in% c("VG"), 9,  12)
      )
    ),
    col = ifelse(
      v$targets$strain == "VV",
      VV_Color,
      ifelse(
        v$targets$strain == "GG",
        GG_Color,
        ifelse(v$targets$strain == "VG",
               VG_Color, GV_Color)
      )
    ),
    gene.selection = "common"
  )

# Adding legends
legend(
  "topleft",
  pch = c(17, 15, 9, 12),
  col = c(VV_Color, GG_Color, VG_Color, GV_Color),
  legend = c("VV", "GG", "VG", "GV")
)
dev.off()
dev.off()
# MDS plot of the first two dimensions using top50 genes. The dimensions are in the mds object.
png(filename = "figures/wilson_AvgREFs_mds_2and3.png",
    width = 450,
    height = 450)
mds <-
  plotMDS(
    v,
    top = 50,
    ndim = 10,
    dim.plot = c(2, 3),
    plot = TRUE,
    cex = 2,
    pch = ifelse(
      v$targets$strain %in% c("VV"),
      17,
      ifelse(
        v$targets$strain %in% c("GG"),
        15,
        ifelse(v$targets$strain %in% c("VG"), 9,  12)
      )
    ),
    col = ifelse(
      v$targets$strain == "VV",
      VV_Color,
      ifelse(
        v$targets$strain == "GG",
        GG_Color,
        ifelse(v$targets$strain == "VG",
               VG_Color, GV_Color)
      )
    ),
    gene.selection = "common"
  )

# Adding legends
legend(
  "topleft",
  pch = c(17, 15, 9, 12),
  col = c(VV_Color, GG_Color, VG_Color, GV_Color),
  legend = c("VV", "GG", "VG", "GV")
)
dev.off()
dev.off()
# Select most variable genes based on the biological coefficient of variance
# (mean scaled)

# Voom transformed counts
voomCounts <- v$E
voomCountsMatrix <- data.matrix(voomCounts, rownames.force = NA)

# Setting the N of genes to use
ntop = length(dge$genes$Geneid)

# Sorting by the coefficient of variance
means <- rowMeans(voomCountsMatrix)
Pvars <- rowVars(voomCountsMatrix)
cv2 <- Pvars / means ^ 2
select <-
  order(cv2, decreasing = TRUE)[seq_len(min(ntop, length(cv2)))]
head(select)
highly_variable_exp <- ((voomCountsMatrix)[select,])
dim(highly_variable_exp)

# Running PCA
pca_exp <- prcomp(t(highly_variable_exp), scale = F, center = T)
# scale a logical value indicating whether the variables should be scaled to have unit variance before the analysis takes place.
# a logical value indicating whether the variables should be shifted to be zero centered.
head(pca_exp$x)[, 1:3]

summary(pca_exp)
# Dataframe with the first 10 PCs
dim1_10 <- data.frame(pca_exp$x[, 1:10])
# Adding metadata
pcaWithMetadata <- merge(dim1_10, dge$samples, by = 0, all = TRUE)
pcaWithMetadata$strain <- factor(pcaWithMetadata$strain,
                                 levels = c("VV", "GG", "VG", "GV", NA))

# Plotting
png(filename = "figures/wilson_AvgREFs_PCA.png",
    width = 650,
    height = 650)
ggplot(data = pcaWithMetadata, aes(
  x = PC1,
  y = PC2,
  shape = strain,
  color = strain
)) +
  geom_point(size = 8) +
  theme_bw() +
  xlim(-150, 150) +
  ylim(-150, 150) +
  scale_color_manual(values = c(VV_Color, GG_Color,
                                VG_Color, GV_Color,
                                "azure3")) +
  scale_shape_manual(values = c(VV_Shape, GG_Shape,
                                VG_Shape, GV_Shape)) +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 30),
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22)
  ) +
  #guides(color = guide_legend(order = 2)) +
  theme(legend.title = element_blank()) +
  xlab("PC1 (61.69)%") +
  ylab("PC2 (35.4%)")
dev.off()
dev.off()

# Differential expression analysis with limma
design <- model.matrix( ~ 0 + v$targets$strain)
colnames(design) <-
  gsub("v\\$targets\\$strain", "", colnames(design))
head(design)

# Running voom again with the new design matrix.
v <- voomWithQualityWeights(dge, design, plot = TRUE)
fit <- lmFit(v, design)
# Contrast design for differential expression
# Defining pairwise comparisons
contrasts <- makeContrasts(
  VV_vs_GG = VV - GG,
  VV_vs_VG = VV - VG,
  VV_vs_GV = VV - GV,
  GG_vs_VG = GG - VG,
  GG_vs_GV = GG - GV,
  VG_vs_GV = VG - GV,
  levels = colnames(design)
)
head(contrasts)
# Assigning all comparisons to a vector for later
allComparisons <- colnames(contrasts)

# Running contrast analysis
vfit <- contrasts.fit(fit, contrasts = contrasts)
# Looking at N of DEGs with adj. p <0.01 and log2FC>2
sumTable <-
  summary(decideTests(
    vfit,
    adjust.method = "BH",
    p.value = 0.01,
    lfc = 2
  ))

# Computing differential expression based on the empirical Bayes moderation of
# the standard errors towards a common value. Robust = should the estimation of
# the empirical Bayes prior parameters be robustified against outlier sample
# variances?
veBayesFit <- eBayes(vfit, robust = TRUE)
plotSA(veBayesFit, main = "Final model: Mean-variance trend")

# Getting the DEGs. Log2FC of 1 is equivalent to linear fold change of 2.
# Getting summary statistics for all genes
coef = 1
for (i in allComparisons) {
  vTopTableAll <-
    topTable(
      veBayesFit,
      coef = coef,
      n = Inf,
      p.value = 1,
      lfc = 0
    )
  path <-
    paste("./DEGs/DEGs_wilson_AvgREFs_fpkm05_",
          i,
          "_fdr1_lfc0.txt",
          sep = "")
  write.table(vTopTableAll, path, sep = "\t")
  
  # Adj.p<0.05, log2FC>1
  vTopTable1 <-
    topTable(
      veBayesFit,
      coef = coef,
      n = Inf,
      p.value = 0.05,
      lfc = 1
    )
  path <-
    paste("./DEGs/DEGs_wilson_AvgREFs_fpkm05_",
          i,
          "_fdr05_lfc1.txt",
          sep = "")
  write.table(vTopTable1, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>0
  vTopTable1 <-
    topTable(
      veBayesFit,
      coef = coef,
      n = Inf,
      p.value = 0.01,
      lfc = 0
    )
  path <-
    paste("./DEGs/DEGs_wilson_AvgREFs_fpkm05_",
          i,
          "_fdr001_lfc0.txt",
          sep = "")
  write.table(vTopTable1, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>1
  vTopTable2 <-
    topTable(
      veBayesFit,
      coef = coef,
      n = Inf,
      p.value = 0.01,
      lfc = 1
    )
  path <-
    paste("./DEGs/DEGs_wilson_AvgREFs_fpkm05_",
          i,
          "_fdr001_lfc1.txt",
          sep = "")
  write.table(vTopTable2, path, sep = "\t")
  
  # Adj.p<0.01, log2FC>2
  vTopTable3 <-
    topTable(
      veBayesFit,
      coef = coef,
      n = Inf,
      p.value = 0.01,
      lfc = 2
    )
  path <-
    paste("./DEGs/DEGs_wilson_AvgREFs_fpkm05_",
          i,
          "_fdr001_lfc2.txt",
          sep = "")
  write.table(vTopTable3, path, sep = "\t")
  
  coef = coef + 1
}

# Venn diagrams of DEGs detected in the strain comparisons
VV_GG_DEG <-
  rownames(read.table(
    "./DEGs/DEGs_wilson_AvgREFs_fpkm05_VV_vs_GG_fdr001_lfc2.txt"
  ))
VV_VG_DEG <-
  rownames(read.table(
    "./DEGs/DEGs_wilson_AvgREFs_fpkm05_VV_vs_VG_fdr001_lfc2.txt"
  ))
VV_GV_DEG <-
  rownames(read.table(
    "./DEGs/DEGs_wilson_AvgREFs_fpkm05_VV_vs_GV_fdr001_lfc2.txt"
  ))

venn.plot <- VennDiagram:::venn.diagram(
  x = list(
    "VV_GG_DEG" = VV_GG_DEG,
    "VV_VG_DEG" = VV_VG_DEG,
    "VV_GV_DEG" = VV_GV_DEG
  ),
  filename = "wilson_AvgREFs_VV_comparisons_venn.tiff",
  scaled = TRUE,
  col = "transparent",
  fill = c(VV_Color, VG_Color, GV_Color),
  main.pos = c(0.5, 0.5, 0.5),
  cex = 1.5,
  cat.cex = 1.5,
  main.cex = 2,
  cat.default.pos = "outer",
  cat.pos = c(-15, 15, 180),
  cat.dist = c(0.05, 0.05, 0.05),
  cat.fontfamily = "sans",
  main = "",
  fontfamily = "sans",
  na = "remove",
  inverted = FALSE
)

venn.plot

# GG comparison venn
VV_GG_DEG <-
  rownames(read.table(
    "./DEGs/DEGs_wilson_AvgREFs_fpkm05_VV_vs_GG_fdr001_lfc2.txt"
  ))
GG_VG_DEG <-
  rownames(read.table(
    "./DEGs/DEGs_wilson_AvgREFs_fpkm05_GG_vs_VG_fdr001_lfc2.txt"
  ))
GG_GV_DEG <-
  rownames(read.table(
    "./DEGs/DEGs_wilson_AvgREFs_fpkm05_GG_vs_GV_fdr001_lfc2.txt"
  ))

venn.plot <- VennDiagram:::venn.diagram(
  x = list(
    "VV_GG_DEG" = VV_GG_DEG,
    "GG_VG_DEG" = GG_VG_DEG,
    "GG_GV_DEG" = GG_GV_DEG
  ),
  filename = "wilson_AvgREFs_GG_comparisons_venn.tiff",
  scaled = TRUE,
  col = "transparent",
  fill = c(GG_Color, VG_Color, GV_Color),
  main.pos = c(0.5, 0.5, 0.5),
  cex = 1.5,
  cat.cex = 1.5,
  main.cex = 2,
  cat.default.pos = "outer",
  cat.pos = c(-15, 15, 180),
  cat.dist = c(0.05, 0.05, 0.05),
  cat.fontfamily = "sans",
  main = "",
  fontfamily = "sans",
  na = "remove",
  inverted = FALSE
)

venn.plot

# Volcano plots.
# VV vs GG
VV_vs_GG_DEG <-
  read.table(
    "./DEGs/DEGs_wilson_AvgREFs_fpkm05_VV_vs_GG_fdr1_lfc0.txt",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = F
  )
VV_vs_GG_DF <-
  data.frame(VV_vs_GG_DEG$adj.P.Val,
             VV_vs_GG_DEG$logFC,
             VV_vs_GG_DEG$Chr,
             VV_vs_GG_DEG$Geneid)
colnames(VV_vs_GG_DF) <- c("adj.P.Val", "logFC", "Chr", "Geneid")

VV_vs_GG_DF_Sig <-
  VV_vs_GG_DF[(abs(VV_vs_GG_DF$logFC) >= 2 &
                 VV_vs_GG_DF$adj.P.Val <= 0.01), ]$Geneid

# Finding stain bias genes, assigning color values
nonSig <- subset(VV_vs_GG_DF,!(Geneid %in% VV_vs_GG_DF_Sig))
nonSig <- cbind(nonSig , rep(1, nrow(nonSig)))
colnames(nonSig)[5] <- "Color"

up_VV <-
  subset(
    VV_vs_GG_DF,
    VV_vs_GG_DF$logFC >= 2 &
      VV_vs_GG_DF$adj.P.Val <= 0.01  & (Geneid %in% VV_vs_GG_DF_Sig)
  )
up_VV <- cbind(up_VV , rep(2, nrow(up_VV)))
colnames(up_VV)[5] <- "Color"

up_GG <-
  subset(
    VV_vs_GG_DF,
    VV_vs_GG_DF$logFC <= -2 &
      VV_vs_GG_DF$adj.P.Val <= 0.01  & (Geneid %in% VV_vs_GG_DF_Sig)
  )
up_GG <- cbind(up_GG , rep(3, nrow(up_GG)))
colnames(up_GG)[5] <- "Color"

dfPlot <- rbind(nonSig, up_VV, up_GG)
dfPlot$Color <- as.factor(dfPlot$Color)

# Constructing the plot object.
p <-
  ggplot(data = dfPlot, aes(
    x = logFC,
    y = -log10(adj.P.Val),
    color = Color
  )) +
  geom_point(alpha = 0.5, size = 8) +
  theme_bw() +
  theme(legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 12)) +
  scale_color_manual(values = c("azure3", GG_Color, VV_Color)) +
  labs(
    title = "VV vs GG",
    x = expression(log[2](FC)),
    y = expression(-log[10] ~ "(FDR-adjusted " ~ italic("p") ~ "-value)")
  ) +
  theme(axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))
# Adding lines for significance thresholds
p
png(filename = "figures/wilson_AvgREFs_VV_vs_GG_valcano.png",
    width = 650,
    height = 650)
p + geom_hline(yintercept = 2,
               colour = "#000000",
               linetype = "dashed") +
  geom_vline(xintercept = 2,
             colour = "#000000",
             linetype = "dashed") +
  geom_vline(xintercept = -2,
             colour = "#000000",
             linetype = "dashed")
dev.off()
dev.off()

# VV vs VG
VV_vs_VG_DEG <-
  read.table(
    "./DEGs/DEGs_wilson_AvgREFs_fpkm05_VV_vs_VG_fdr1_lfc0.txt",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = F
  )
VV_vs_VG_DF <-
  data.frame(VV_vs_VG_DEG$adj.P.Val,
             VV_vs_VG_DEG$logFC,
             VV_vs_VG_DEG$Chr,
             VV_vs_VG_DEG$Geneid)
colnames(VV_vs_VG_DF) <- c("adj.P.Val", "logFC", "Chr", "Geneid")

VV_vs_VG_DF_Sig <-
  VV_vs_VG_DF[(abs(VV_vs_VG_DF$logFC) >= 2 &
                 VV_vs_VG_DF$adj.P.Val <= 0.01), ]$Geneid

# Finding stain bias genes, assigning color values
nonSig <- subset(VV_vs_VG_DF,!(Geneid %in% VV_vs_VG_DF_Sig))
nonSig <- cbind(nonSig , rep(1, nrow(nonSig)))
colnames(nonSig)[5] <- "Color"

up_VV <-
  subset(
    VV_vs_VG_DF,
    VV_vs_VG_DF$logFC >= 2 &
      VV_vs_VG_DF$adj.P.Val <= 0.01  & (Geneid %in% VV_vs_VG_DF_Sig)
  )
up_VV <- cbind(up_VV , rep(2, nrow(up_VV)))
colnames(up_VV)[5] <- "Color"

up_VG <-
  subset(
    VV_vs_VG_DF,
    VV_vs_VG_DF$logFC <= -2 &
      VV_vs_VG_DF$adj.P.Val <= 0.01  & (Geneid %in% VV_vs_VG_DF_Sig)
  )
up_VG <- cbind(up_VG , rep(3, nrow(up_VG)))
colnames(up_VG)[5] <- "Color"

dfPlot <- rbind(nonSig, up_VV, up_VG)
dfPlot$Color <- as.factor(dfPlot$Color)

# Constructing the plot object.
p <-
  ggplot(data = dfPlot, aes(
    x = logFC,
    y = -log10(adj.P.Val),
    color = Color
  )) +
  geom_point(alpha = 0.5, size = 8) +
  theme_bw() +
  theme(legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 12)) +
  scale_color_manual(values = c("azure3", VG_Color, VV_Color)) +
  labs(
    title = "VV vs VG",
    x = expression(log[2](FC)),
    y = expression(-log[10] ~ "(FDR-adjusted " ~ italic("p") ~ "-value)")
  ) +
  theme(axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))
# Adding lines for significance thresholds
p
png(filename = "figures/wilson_AvgREFs_VV_vs_VG_valcano.png",
    width = 650,
    height = 650)
p + geom_hline(yintercept = 2,
               colour = "#000000",
               linetype = "dashed") +
  geom_vline(xintercept = 2,
             colour = "#000000",
             linetype = "dashed") +
  geom_vline(xintercept = -2,
             colour = "#000000",
             linetype = "dashed")
dev.off()
dev.off()

# VV vs GV
VV_vs_GV_DEG <-
  read.table(
    "./DEGs/DEGs_wilson_AvgREFs_fpkm05_VV_vs_GV_fdr1_lfc0.txt",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = F
  )
VV_vs_GV_DF <-
  data.frame(VV_vs_GV_DEG$adj.P.Val,
             VV_vs_GV_DEG$logFC,
             VV_vs_GV_DEG$Chr,
             VV_vs_GV_DEG$Geneid)
colnames(VV_vs_GV_DF) <- c("adj.P.Val", "logFC", "Chr", "Geneid")

VV_vs_GV_DF_Sig <-
  VV_vs_GV_DF[(abs(VV_vs_GV_DF$logFC) >= 2 &
                 VV_vs_GV_DF$adj.P.Val <= 0.01), ]$Geneid

# Finding stain bias genes, assigning color values
nonSig <- subset(VV_vs_GV_DF,!(Geneid %in% VV_vs_GV_DF_Sig))
nonSig <- cbind(nonSig , rep(1, nrow(nonSig)))
colnames(nonSig)[5] <- "Color"

up_VV <-
  subset(
    VV_vs_GV_DF,
    VV_vs_GV_DF$logFC >= 2 &
      VV_vs_GV_DF$adj.P.Val <= 0.01  & (Geneid %in% VV_vs_GV_DF_Sig)
  )
up_VV <- cbind(up_VV , rep(2, nrow(up_VV)))
colnames(up_VV)[5] <- "Color"

up_GV <-
  subset(
    VV_vs_GV_DF,
    VV_vs_GV_DF$logFC <= -2 &
      VV_vs_GV_DF$adj.P.Val <= 0.01  & (Geneid %in% VV_vs_GV_DF_Sig)
  )
up_GV <- cbind(up_GV , rep(3, nrow(up_GV)))
colnames(up_GV)[5] <- "Color"

dfPlot <- rbind(nonSig, up_VV, up_GV)
dfPlot$Color <- as.factor(dfPlot$Color)

# Constructing the plot object.
p <-
  ggplot(data = dfPlot, aes(
    x = logFC,
    y = -log10(adj.P.Val),
    color = Color
  )) +
  geom_point(alpha = 0.5, size = 8) +
  theme_bw() +
  theme(legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 12)) +
  scale_color_manual(values = c("azure3", GV_Color, VV_Color)) +
  labs(
    title = "VV vs GV",
    x = expression(log[2](FC)),
    y = expression(-log[10] ~ "(FDR-adjusted " ~ italic("p") ~ "-value)")
  ) +
  theme(axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))
# Adding lines for significance thresholds
png(filename = "figures/wilson_AvgREFs_VV_vs_GV_valcano.png",
    width = 650,
    height = 650)
p + geom_hline(yintercept = 2,
               colour = "#000000",
               linetype = "dashed") +
  geom_vline(xintercept = 2,
             colour = "#000000",
             linetype = "dashed") +
  geom_vline(xintercept = -2,
             colour = "#000000",
             linetype = "dashed")
dev.off()
dev.off()

# GG vs VG
GG_vs_VG_DEG <-
  read.table(
    "./DEGs/DEGs_wilson_AvgREFs_fpkm05_GG_vs_VG_fdr1_lfc0.txt",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = F
  )
GG_vs_VG_DF <-
  data.frame(GG_vs_VG_DEG$adj.P.Val,
             GG_vs_VG_DEG$logFC,
             GG_vs_VG_DEG$Chr,
             GG_vs_VG_DEG$Geneid)
colnames(GG_vs_VG_DF) <- c("adj.P.Val", "logFC", "Chr", "Geneid")

GG_vs_VG_DF_Sig <-
  GG_vs_VG_DF[(abs(GG_vs_VG_DF$logFC) >= 2 &
                 GG_vs_VG_DF$adj.P.Val <= 0.01), ]$Geneid

# Finding stain bias genes, assigning color values
nonSig <- subset(GG_vs_VG_DF,!(Geneid %in% GG_vs_VG_DF_Sig))
nonSig <- cbind(nonSig , rep(1, nrow(nonSig)))
colnames(nonSig)[5] <- "Color"

up_GG <-
  subset(
    GG_vs_VG_DF,
    GG_vs_VG_DF$logFC >= 2 &
      GG_vs_VG_DF$adj.P.Val <= 0.01  & (Geneid %in% GG_vs_VG_DF_Sig)
  )
up_GG <- cbind(up_GG , rep(2, nrow(up_GG)))
colnames(up_GG)[5] <- "Color"

up_VG <-
  subset(
    GG_vs_VG_DF,
    GG_vs_VG_DF$logFC <= -2 &
      GG_vs_VG_DF$adj.P.Val <= 0.01  & (Geneid %in% GG_vs_VG_DF_Sig)
  )
up_VG <- cbind(up_VG , rep(3, nrow(up_VG)))
colnames(up_VG)[5] <- "Color"

dfPlot <- rbind(nonSig, up_GG, up_VG)
dfPlot$Color <- as.factor(dfPlot$Color)

# Constructing the plot object.
p <-
  ggplot(data = dfPlot, aes(
    x = logFC,
    y = -log10(adj.P.Val),
    color = Color
  )) +
  geom_point(alpha = 0.5, size = 8) +
  theme_bw() +
  theme(legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 12)) +
  scale_color_manual(values = c("azure3", VG_Color, GG_Color)) +
  labs(
    title = "GG vs VG",
    x = expression(log[2](FC)),
    y = expression(-log[10] ~ "(FDR-adjusted " ~ italic("p") ~ "-value)")
  ) +
  theme(axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))
# Adding lines for significance thresholds
png(filename = "figures/wilson_AvgREFs_GG_vs_VG_valcano.png",
    width = 650,
    height = 650)
p + geom_hline(yintercept = 2,
               colour = "#000000",
               linetype = "dashed") +
  geom_vline(xintercept = 2,
             colour = "#000000",
             linetype = "dashed") +
  geom_vline(xintercept = -2,
             colour = "#000000",
             linetype = "dashed")
dev.off()
dev.off()

# GG vs GV
GG_vs_GV_DEG <-
  read.table(
    "./DEGs/DEGs_wilson_AvgREFs_fpkm05_GG_vs_GV_fdr1_lfc0.txt",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = F
  )
GG_vs_GV_DF <-
  data.frame(GG_vs_GV_DEG$adj.P.Val,
             GG_vs_GV_DEG$logFC,
             GG_vs_GV_DEG$Chr,
             GG_vs_GV_DEG$Geneid)
colnames(GG_vs_GV_DF) <- c("adj.P.Val", "logFC", "Chr", "Geneid")

GG_vs_GV_DF_Sig <-
  GG_vs_GV_DF[(abs(GG_vs_GV_DF$logFC) >= 2 &
                 GG_vs_GV_DF$adj.P.Val <= 0.01), ]$Geneid

# Finding stain bias genes, assigning color values
nonSig <- subset(GG_vs_GV_DF,!(Geneid %in% GG_vs_GV_DF_Sig))
nonSig <- cbind(nonSig , rep(1, nrow(nonSig)))
colnames(nonSig)[5] <- "Color"

up_GG <-
  subset(
    GG_vs_GV_DF,
    GG_vs_GV_DF$logFC >= 2 &
      GG_vs_GV_DF$adj.P.Val <= 0.01  & (Geneid %in% GG_vs_GV_DF_Sig)
  )
up_GG <- cbind(up_GG , rep(2, nrow(up_GG)))
colnames(up_GG)[5] <- "Color"

up_GV <-
  subset(
    GG_vs_GV_DF,
    GG_vs_GV_DF$logFC <= -2 &
      GG_vs_GV_DF$adj.P.Val <= 0.01  & (Geneid %in% GG_vs_GV_DF_Sig)
  )
up_GV <- cbind(up_GV , rep(3, nrow(up_GV)))
colnames(up_GV)[5] <- "Color"

dfPlot <- rbind(nonSig, up_GG, up_GV)
dfPlot$Color <- as.factor(dfPlot$Color)

# Constructing the plot object.
p <-
  ggplot(data = dfPlot, aes(
    x = logFC,
    y = -log10(adj.P.Val),
    color = Color
  )) +
  geom_point(alpha = 0.5, size = 8) +
  theme_bw() +
  theme(legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 12)) +
  scale_color_manual(values = c("azure3", GV_Color, GG_Color)) +
  labs(
    title = "GG vs GV",
    x = expression(log[2](FC)),
    y = expression(-log[10] ~ "(FDR-adjusted " ~ italic("p") ~ "-value)")
  ) +
  theme(axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))
# Adding lines for significance thresholds
png(filename = "figures/wilson_AvgREFs_GG_vs_GV_valcano.png",
    width = 650,
    height = 650)
p + geom_hline(yintercept = 2,
               colour = "#000000",
               linetype = "dashed") +
  geom_vline(xintercept = 2,
             colour = "#000000",
             linetype = "dashed") +
  geom_vline(xintercept = -2,
             colour = "#000000",
             linetype = "dashed")
dev.off()
dev.off()

# VG vs GV
VG_vs_GV_DEG <-
  read.table(
    "./DEGs/DEGs_wilson_AvgREFs_fpkm05_VG_vs_GV_fdr1_lfc0.txt",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = F
  )
VG_vs_GV_DF <-
  data.frame(VG_vs_GV_DEG$adj.P.Val,
             VG_vs_GV_DEG$logFC,
             VG_vs_GV_DEG$Chr,
             VG_vs_GV_DEG$Geneid)
colnames(VG_vs_GV_DF) <- c("adj.P.Val", "logFC", "Chr", "Geneid")

VG_vs_GV_DF_Sig <-
  VG_vs_GV_DF[(abs(VG_vs_GV_DF$logFC) >= 2 &
                 VG_vs_GV_DF$adj.P.Val <= 0.01), ]$Geneid

# Finding stain bias genes, assigning color values
nonSig <- subset(VG_vs_GV_DF,!(Geneid %in% VG_vs_GV_DF_Sig))
nonSig <- cbind(nonSig , rep(1, nrow(nonSig)))
colnames(nonSig)[5] <- "Color"

up_VG <-
  subset(
    VG_vs_GV_DF,
    VG_vs_GV_DF$logFC >= 2 &
      VG_vs_GV_DF$adj.P.Val <= 0.01  & (Geneid %in% VG_vs_GV_DF_Sig)
  )
up_VG <- cbind(up_VG , rep(2, nrow(up_VG)))
colnames(up_VG)[5] <- "Color"

up_GV <-
  subset(
    VG_vs_GV_DF,
    VG_vs_GV_DF$logFC <= -2 &
      VG_vs_GV_DF$adj.P.Val <= 0.01  & (Geneid %in% VG_vs_GV_DF_Sig)
  )
up_GV <- cbind(up_GV , rep(3, nrow(up_GV)))
colnames(up_GV)[5] <- "Color"

dfPlot <- rbind(nonSig, up_VG, up_GV)
dfPlot$Color <- as.factor(dfPlot$Color)

# Constructing the plot object.
p <-
  ggplot(data = dfPlot, aes(
    x = logFC,
    y = -log10(adj.P.Val),
    color = Color
  )) +
  geom_point(alpha = 0.5, size = 8) +
  theme_bw() +
  theme(legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 12)) +
  scale_color_manual(values = c("azure3", GV_Color, VG_Color)) +
  labs(
    title = "VG vs GV",
    x = expression(log[2](FC)),
    y = expression(-log[10] ~ "(FDR-adjusted " ~ italic("p") ~ "-value)")
  ) +
  theme(axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))
# Adding lines for significance thresholds
png(filename = "figures/wilson_AvgREFs_VG_vs_GV_valcano.png",
    width = 650,
    height = 650)
p + geom_hline(yintercept = 2,
               colour = "#000000",
               linetype = "dashed") +
  geom_vline(xintercept = 2,
             colour = "#000000",
             linetype = "dashed") +
  geom_vline(xintercept = -2,
             colour = "#000000",
             linetype = "dashed")
dev.off()
dev.off()
