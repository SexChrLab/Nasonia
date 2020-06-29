# ==============================================================================
# Author(s) : Heini M Natri, heini.natri@gmail.com
# Date: June 2020
# Description: Plotting regression between Nasonia datasets
# ==============================================================================

# Load libraries
library(RColorBrewer)
library(edgeR)
library(ggplot2)
library(ggExtra)
library(grid)
library(gridExtra)

# Setting the working directory
setwd("/Users/hnatri/Dropbox (ASU)/Nasonia/")

# Defining colors
colorPalette <- brewer.pal(8, "Set2")
VV_Color <- colorPalette[1]
GG_Color <- colorPalette[2]
VG_Color <- colorPalette[3]
GV_Color <- colorPalette[4]
VV_Shape <- c(17)
GG_Shape <- c(15)
VG_Shape <- c(9)
GV_Shape <- c(12)

# Reading data
wilsonCounts <- read.table("wilson_counts.txt", sep = "\t", header=TRUE, check.names=FALSE)
clarkCounts <- read.table("clark_counts.txt", sep = "\t", header=TRUE, check.names=FALSE)
wilsonGirCounts <- read.table("wilson_counts_VgirRef.txt", sep = "\t", header=TRUE, check.names=FALSE)
clarkGirCounts <- read.table("clark_counts_VgirRef.txt", sep = "\t", header=TRUE, check.names=FALSE)

# Taking the mean between reference genomes
clark_df <- data.frame(clarkCounts, clarkGirCounts)
colnames = colnames(clarkCounts)
clarkMeanCounts <- sapply(colnames, function(x) rowMeans(clark_df [, grep(x, names(clark_df))] )  )

wilson_df <- data.frame(wilsonCounts, wilsonGirCounts)
colnames = colnames(wilsonCounts)
wilsonMeanCounts <- sapply(colnames, function(x) rowMeans(wilson_df [, grep(x, names(wilson_df))] )  )

# Combining counts
bothCounts <- merge(wilsonCounts, clarkCounts, by=0)
bothCounts <- bothCounts[, !(colnames(bothCounts) == "Row.names")]

metadata <- read.csv("metadata.csv")
metadata <- metadata[match(colnames(bothCounts),
                           metadata$SampleID),]
metadata$Dataset_Cross <- paste(metadata$Dataset, metadata$Cross, sep="_")

samplesToRemove <- c("SRR2773798") # update depending on comparison being made
metadata <- metadata[!(metadata$SampleID %in% samplesToRemove),]
bothCounts <- bothCounts[, !(colnames(bothCounts) %in% samplesToRemove)]

# Importing gene annotations
genes <- read.table("014451_featurecounts.tsv", header=TRUE, sep="\t")
genes <- data.frame(genes)
genes <- genes[,colnames(genes) %in% c("Geneid", "Chr", "Start", "End", "Strand", "Length")]

# Vectors of samples in each group
WilsonVV <- as.vector(metadata[metadata$Dataset_Cross == "Wilson_VV" , ]$SampleID)
WilsonGG <- as.vector(metadata[metadata$Dataset_Cross == "Wilson_GG" , ]$SampleID)
WilsonVG <- as.vector(metadata[metadata$Dataset_Cross == "Wilson_VG" , ]$SampleID)
WilsonGV <- as.vector(metadata[metadata$Dataset_Cross == "Wilson_GV" , ]$SampleID)
ClarkVV <- as.vector(metadata[metadata$Dataset_Cross == "Clark_VV" , ]$SampleID)
ClarkGG <- as.vector(metadata[metadata$Dataset_Cross == "Clark_GG" , ]$SampleID)
ClarkVG <- as.vector(metadata[metadata$Dataset_Cross == "Clark_VG" , ]$SampleID)
ClarkGG <- as.vector(metadata[metadata$Dataset_Cross == "Clark_GV" , ]$SampleID)

# Calculating mean logCPM for each group of samples in the Wilson and Clark datasets
dge <- DGEList(counts=bothCounts, genes=genes)
colnames(dge) <- colnames(bothCounts)
dge$samples$Dataset_Cross <- metadata$Dataset_Cross

logCPM <- cpm(dge, log = T)
meanlogCPM <- data.frame(matrix(NA, nrow = length(rownames(logCPM)), ncol = 8))
colnames(meanlogCPM) <- c("WilsonVV", "WilsonGG", "WilsonVG", "WilsonGV", "ClarkVV", "ClarkGG", "ClarkVG", "ClarkGV")

meanlogCPM$WilsonVV <- apply(as.data.frame(logCPM)
                             [(dge$samples$Dataset_Cross=="Wilson_VV")],
                             1, mean, na.rm=TRUE)
meanlogCPM$WilsonGG <- apply(as.data.frame(logCPM)
                             [(dge$samples$Dataset_Cross=="Wilson_GG")],
                             1, mean, na.rm=TRUE)
meanlogCPM$WilsonVG <- apply(as.data.frame(logCPM)
                             [(dge$samples$Dataset_Cross=="Wilson_VG")],
                             1, mean, na.rm=TRUE)
meanlogCPM$WilsonGV <- apply(as.data.frame(logCPM)
                             [(dge$samples$Dataset_Cross=="Wilson_GV")],
                             1, mean, na.rm=TRUE)
meanlogCPM$ClarkVV <- apply(as.data.frame(logCPM)
                            [(dge$samples$Dataset_Cross=="Clark_VV")],
                            1, mean, na.rm=TRUE)
meanlogCPM$ClarkGG <- apply(as.data.frame(logCPM)
                            [(dge$samples$Dataset_Cross=="Clark_GG")],
                            1, mean, na.rm=TRUE)
meanlogCPM$ClarkVG <- apply(as.data.frame(logCPM)
                            [(dge$samples$Dataset_Cross=="Clark_VG")],
                            1, mean, na.rm=TRUE)
meanlogCPM$ClarkGV <- apply(as.data.frame(logCPM)
                            [(dge$samples$Dataset_Cross=="Clark_GV")],
                            1, mean, na.rm=TRUE)

# Plotting the Wilson and Clark mean expression values for GG samples
# Calculating R2
p1_r2 <- cor(meanlogCPM$WilsonGG, meanlogCPM$ClarkGG,  method = "pearson", use = "complete.obs")
label <- paste('R-squared', "=", format(round(p1_r2, 2), nsmall = 2), sep="")

p1 <- ggplot(data = meanlogCPM, aes(x = WilsonGG, y = ClarkGG )) +
  geom_point(alpha = 0.5, size = 4, color = GG_Color) +
  theme( legend.position = "none") +
  theme_bw() +
  #xlim(c(-10, 10)) + ylim(c(0, 30)) ++
  labs(x="Wilson GG mean logCPM", y="Clark GG mean logCPM") +
  theme(axis.title.x=element_text(size=12), 
        axis.text.x=element_text(size=10)) +
  theme(axis.title.y=element_text(size=12),
        axis.text.y=element_text(size=10)) +
  scale_x_continuous(breaks=c(-5,0,5,10,15), limits = c(-5, 15)) +
  scale_y_continuous(breaks=c(-5,0,5,10,15), limits = c(-5, 15)) +
  geom_smooth(data = meanlogCPM, aes(x = WilsonGG, y = ClarkGG), method='lm',formula=y~x, color="black") +
  annotate(geom="text", x=1, y=14, label=label, color="black")

p1

# Plotting the Wilson and Clark mean expression values for VV samples
# Calculating R2
p2_r2 <- cor(meanlogCPM$WilsonVV, meanlogCPM$ClarkVV,  method = "pearson", use = "complete.obs")
label <- paste('R-squared', "=", format(round(p2_r2, 2), nsmall = 2), sep="")

p2 <- ggplot(data = meanlogCPM, aes(x = WilsonVV, y = ClarkVV )) +
  geom_point(alpha = 0.5, size = 4, color = VV_Color) +
  theme( legend.position = "none") +
  theme_bw() +
  #xlim(c(-10, 10)) + ylim(c(0, 30)) ++
  labs(x="Wilson VV mean logCPM", y="Clark VV mean logCPM") +
  theme(axis.title.x=element_text(size=12), 
        axis.text.x=element_text(size=10)) +
  theme(axis.title.y=element_text(size=12),
        axis.text.y=element_text(size=10)) +
  scale_x_continuous(breaks=c(-5,0,5,10,15), limits = c(-5, 15)) +
  scale_y_continuous(breaks=c(-5,0,5,10,15), limits = c(-5, 15)) +
  geom_smooth(data = meanlogCPM, aes(x = WilsonVV, y = ClarkVV), method='lm',formula=y~x, color="black") +
  annotate(geom="text", x=1, y=14, label=label, color="black")

p2

# Plotting the Wilson and Clark mean expression values for GV samples
# Calculating R2
p3_r2 <- cor(meanlogCPM$WilsonGV, meanlogCPM$ClarkGV,  method = "pearson", use = "complete.obs")
label <- paste('R-squared', "=", format(round(p3_r2, 2), nsmall = 2), sep="")

p3 <- ggplot(data = meanlogCPM, aes(x = WilsonGV, y = ClarkGV )) +
  geom_point(alpha = 0.5, size = 4, color = GV_Color) +
  theme( legend.position = "none") +
  theme_bw() +
  #xlim(c(-10, 10)) + ylim(c(0, 30)) ++
  labs(x="Wilson GV mean logCPM", y="Clark GV mean logCPM") +
  theme(axis.title.x=element_text(size=12), 
        axis.text.x=element_text(size=10)) +
  theme(axis.title.y=element_text(size=12),
        axis.text.y=element_text(size=10)) +
  scale_x_continuous(breaks=c(-5,0,5,10,15), limits = c(-5, 15)) +
  scale_y_continuous(breaks=c(-5,0,5,10,15), limits = c(-5, 15)) +
  geom_smooth(data = meanlogCPM, aes(x = WilsonGV, y = ClarkGV), method='lm',formula=y~x, color="black") +
  annotate(geom="text", x=1, y=14, label=label, color="black")

p3

# Plotting the Wilson and Clark mean expression values for VG samples
# Calculating R2
p4_r2 <- cor(meanlogCPM$WilsonVG, meanlogCPM$ClarkVG,  method = "pearson", use = "complete.obs")
label <- paste('R-squared', "=", format(round(p4_r2, 2), nsmall = 2), sep="")

p4 <- ggplot(data = meanlogCPM, aes(x = WilsonVG, y = ClarkVG )) +
  geom_point(alpha = 0.5, size = 4, color = VG_Color) +
  theme( legend.position = "none") +
  theme_bw() +
  #xlim(c(-10, 10)) + ylim(c(0, 30)) ++
  labs(x="Wilson VG mean logCPM", y="Clark VG mean logCPM") +
  theme(axis.title.x=element_text(size=12), 
        axis.text.x=element_text(size=10)) +
  theme(axis.title.y=element_text(size=12),
        axis.text.y=element_text(size=10)) +
  scale_x_continuous(breaks=c(-5,0,5,10,15), limits = c(-5, 15)) +
  scale_y_continuous(breaks=c(-5,0,5,10,15), limits = c(-5, 15)) +
  geom_smooth(data = meanlogCPM, aes(x = WilsonVG, y = ClarkVG), method='lm',formula=y~x, color="black") +
  annotate(geom="text", x=1, y=14, label=label, color="black")

p4

# Between species
# Wilson VV and GG
# Calculating R2
p5_r2 <- cor(meanlogCPM$WilsonVV, meanlogCPM$WilsonGG,  method = "pearson", use = "complete.obs")
label <- paste('R-squared', "=", format(round(p5_r2, 2), nsmall = 2), sep="")

p5 <- ggplot(data = meanlogCPM, aes(x = WilsonVV, y = WilsonGG )) +
  geom_point(alpha = 0.5, size = 4, color = "azure3") +
  theme( legend.position = "none") +
  theme_bw() +
  #xlim(c(-10, 10)) + ylim(c(0, 30)) ++
  labs(x="Wilson VV mean logCPM", y="Wilson GG mean logCPM") +
  theme(axis.title.x=element_text(size=12), 
        axis.text.x=element_text(size=10)) +
  theme(axis.title.y=element_text(size=12),
        axis.text.y=element_text(size=10)) +
  scale_x_continuous(breaks=c(-5,0,5,10,15), limits = c(-5, 15)) +
  scale_y_continuous(breaks=c(-5,0,5,10,15), limits = c(-5, 15)) +
  geom_smooth(data = meanlogCPM, aes(x = WilsonVV, y = WilsonGG), method='lm',formula=y~x, color="black") +
  annotate(geom="text", x=1, y=14, label=label, color="black")

p5

# Clark VV and GG
# Calculating R2
p6_r2 <- cor(meanlogCPM$ClarkVV, meanlogCPM$ClarkGG,  method = "pearson", use = "complete.obs")
label <- paste('R-squared', "=", format(round(p6_r2, 2), nsmall = 2), sep="")

p6 <- ggplot(data = meanlogCPM, aes(x = ClarkVV, y = ClarkGG )) +
  geom_point(alpha = 0.5, size = 4, color = "azure3") +
  theme( legend.position = "none") +
  theme_bw() +
  #xlim(c(-10, 10)) + ylim(c(0, 30)) ++
  labs(x="Clark VV mean logCPM", y="Clark GG mean logCPM") +
  theme(axis.title.x=element_text(size=12), 
        axis.text.x=element_text(size=10)) +
  theme(axis.title.y=element_text(size=12),
        axis.text.y=element_text(size=10)) +
  scale_x_continuous(breaks=c(-5,0,5,10,15), limits = c(-5, 15)) +
  scale_y_continuous(breaks=c(-5,0,5,10,15), limits = c(-5, 15)) +
  geom_smooth(data = meanlogCPM, aes(x = ClarkVV, y = ClarkGG), method='lm',formula=y~x, color="black") +
  annotate(geom="text", x=1, y=14, label=label, color="black")

p6

# Saving to a file
ggsave("Clark_Sayres_regression.pdf", 
       grid.arrange(p1, p2, p3, p4, p5, p6, ncol=2, widths=c(3, 3)),
       width = 6, height = 9)


