setwd("~/Dropbox (ASU)/NASONIA_take3/")
library(VennDiagram)
library(RColorBrewer)
viralPalette <- brewer.pal(8, "Set2")
VV_Color <- viralPalette[1]
GG_Color <- viralPalette[2]
VG_Color <- viralPalette[3]
GV_Color <- viralPalette[4]

VV_Shape <- c(17)
GG_Shape <- c(15)
VG_Shape <- c(9)
GV_Shape <- c(12)

# ASE
clark_VG_GV_ASE <-
  read.table("ASE/ASEReadCounter/avgRefGenomes/Clark_ASE_fisherExact.txt",
             header = TRUE)
clark_VG_GV_ASE <-
  subset(clark_VG_GV_ASE, clark_VG_GV_ASE$FDRp <= 0.01)
clark_VG_GV_ASE <- clark_VG_GV_ASE$geneID

wilson_VG_GV_ASE <-
  read.table("ASE/ASEReadCounter/avgRefGenomes/Wilson_ASE_fisherExact.txt",
             header = TRUE)
wilson_VG_GV_ASE <-
  subset(wilson_VG_GV_ASE, wilson_VG_GV_ASE$FDRp <= 0.01)
wilson_VG_GV_ASE <- wilson_VG_GV_ASE$geneID

# DEGs
# GG comparison venn
clark_VG_GV_DEG <-
  read.table("DE_limmaVoom/DEGs/AvgREFs/DEGs_clark_AvgREFs_fpkm05_VG_vs_GV_fdr001_lfc2.txt")
clark_VG_GV_DEG <- clark_VG_GV_DEG$Geneid

wilson_VG_GV_DEG <-
  read.table(
    "DE_limmaVoom/DEGs/AvgREFs/DEGs_wilson_AvgREFs_fpkm05_VG_vs_GV_fdr001_lfc2.txt"
  )
wilson_VG_GV_DEG <- wilson_VG_GV_DEG$Geneid


venn.plot <- VennDiagram:::venn.diagram(
  x = list("ASE" = clark_VG_GV_ASE, "DE" = clark_VG_GV_DEG),
  filename = "VG_to_GV_ASE_DE_clark.tiff",
  scaled = TRUE,
  col = "transparent",
  fill = c("gray", "gray"),
  main.pos = c(0.5, 0.5),
  cex = 1.5,
  cat.cex = .5,
  main.cex = 2,
  cat.pos = 180,
  cat.default.pos = "outer",
  #cat.fontfamily = "sans",
  main = "",
  fontfamily = "sans",
  na = "remove",
  rotation.degree = 180,
  ext.text = FALSE,
  print.mode = c("raw", "percent"),
  margin = .6
)

venn.plot

venn.plot <- VennDiagram:::venn.diagram(
  x = list("ASE" = wilson_VG_GV_ASE, "DE" = wilson_VG_GV_DEG),
  filename = "VG_to_GV_ASE_DE_wilson.tiff",
  scaled = TRUE,
  col = "transparent",
  fill = c("gray", "gray"),
  main.pos = c(0.5, 0.5),
  cex = 1.5,
  cat.cex = .5,
  main.cex = 2,
  cat.pos = 180,
  cat.default.pos = "outer",
  #cat.fontfamily = "sans",
  main = "",
  fontfamily = "sans",
  na = "remove",
  rotation.degree = 180,
  ext.text = FALSE,
  print.mode = c("raw", "percent"),
  margin = .6
)

venn.plot
