
## Script to annotate clusters using known marker genes

## loading libraries
library(tidyverse)
library(Seurat)

BM_combined <- readRDS("~/BM_Combined_Relabeled.RDS")

DefaultAssay(BM_combined) <- "RNA"
FeaturePlot(BM_combined, features = "MME", pt.size = 0.01)

# Undifferentiated - HSC/Progenitor
p1<- FeaturePlot(BM_combined, features = c("EGR1", "MSI2", "CD38", "CD34", "PROM1", "EGFL7"), pt.size = 0.01, combine = F)
p1 <- lapply(X=p1, FUN = function(x) AugmentPlot(x))
CombinePlots(p1)

### Myeloid - GMP/ProMono
p1 <- FeaturePlot(BM_combined, features = c("MPO","ELANE","CTSG","AZU1","LYST","LYZ","CEBPD","MNDA"), pt.size = 0.01, min.cutoff = 'q9', combine = F)
p1 <- lapply(X=p1, FUN = function(x) AugmentPlot(x))
CombinePlots(p1)


### Myeloid - Mono/cDC/pDC
p1 <- lapply(X=FeaturePlot(BM_combined, features = c("FCER1G", "FCN", "CD14", "C5AR1", "CLEC4A", "CLEC10A"),
                           pt.size = 0.01, min.cutoff = 'q1', combine = F), FUN = function(x) AugmentPlot(x))
print(p1)
lapply(X=list(FeaturePlot(BM_combined, features = c("FCER1A", "CLEC4C", "PTPRS", "IRF8", "TCF4"),
                          pt.size = 0.01, min.cutoff = 'q9')), FUN = function(x) AugmentPlot(x))
FeaturePlot(BM_combined, features = c("FCER1A"),
            pt.size = 0.01, min.cutoff = 'q9')

## Erythroid - Late/Early
DefaultAssay(BM_combined) <- "integrated"
p1 <- lapply(X=FeaturePlot(BM_combined, features = c("CSF1", "KIT", "HBB", "HBD", "GYPA"),
                           pt.size = 0.01, min.cutoff = 'q1', combine = F), FUN = function(x) AugmentPlot(x))
CombinePlots(p1)

## Lymphoid - ProB/B/Plasma
lapply(X=list(FeaturePlot(BM_combined, features = c("CD24", "EBF1", "MME", "VPREB1", "PAX5", "CD19"),
                          pt.size = 0.01, min.cutoff = 'q9')), FUN = function(x) AugmentPlot(x))
lapply(X=list(FeaturePlot(BM_combined, features = c("CD79A", "MS4A1", "BANK1", "MZB1", "IGLL5", "JCHAIN"),
                          pt.size = 0.01, min.cutoff = 'q9')), FUN = function(x) AugmentPlot(x))
FeaturePlot(BM_combined, features = "CCL24", split.by="Condition", min.cutoff = "q9")
## Lymphoid - T, CTL, K
lapply(X=list(FeaturePlot(BM_combined, features = c("CD3D", "CD3G", "IL32", "IL7R", "TCF7", "CCL5"),
                          pt.size = 0.01, min.cutoff = 'q9')), FUN = function(x) AugmentPlot(x))
lapply(X=list(FeaturePlot(BM_combined, features = c("GZMK", "CD8A", "KLRB1", "KLRD1",  "GZMB", "NCAM1"),
                          pt.size = 0.01, min.cutoff = 'q9')), FUN = function(x) AugmentPlot(x))
FeaturePlot(BM_combined, features = c("CD14", "CD68", "TFRC", "GYPA"),
            pt.size = 0.01, min.cutoff = 'q9')
FeaturePlot(BM_combined, features = c("GP1BA", "ITGA2B", "ITGB3", "ITGAM"),
            pt.size = 0.01, min.cutoff = 'q9')
DimPlot(BM_combined, label=F, group.by = "orig.ident") #+ guides(color=FALSE)
### Heatmap

Genes <- c("MEIS1", "EGR1", "MSI2", "CD34", "PROM1", "EGFL7", "CD38",
           "MPO","ELANE","CTSG","AZU1","LYST","LYZ","CEBPD","MNDA",
           "FCER1G", "FCN1", "CD14", "C5AR1", "CLEC4A", "CLEC10A",
           "FCER1A", "CLEC4C", "PTPRS", "IRF8", "TCF4",
           "CSF1", "KIT", "HBB", "HBD", "GYPA",
           "CD24", "EBF1", "MME", "VPREB1", "PAX5", "CD19",
           "CD79A", "MS4A1", "BANK1", "MZB1", "IGLL5", "JCHAIN",
           "CD3D", "CD3G", "IL32", "IL7R", "TCF7", "CCL5",
           "GZMK", "CD8A", "KLRB1", "KLRD1", "GZMB", "NCAM1")


DefaultAssay(BM_combined) <- "integrated"
DimPlot(BM_combined)
p1 <- DoHeatmap(BM_combined, cells = sample(colnames(BM_combined), 10000), features = Genes, label = T, size = 6, draw.lines = F) + NoLegend()

## Annotation of column
AllGenes <- read.table("~/Genes.txt", 
                       sep = "\t", 
                       stringsAsFactors = F,
                       header = F)
colnames(AllGenes) <- c("Genes", "CellType")
AllGenes$Genes <- gsub("(.*)[[:space:]]\\(.*", "\\1", AllGenes$Genes)
Groups <- unique(AllGenes$CellType)

MegEryth1 <- AllGenes$Genes[AllGenes$CellType == Groups[1]]
MegEryth2 <- AllGenes$Genes[AllGenes$CellType == Groups[2]]
MEP <- AllGenes$Genes[AllGenes$CellType == Groups[3]]
Mega <- AllGenes$Genes[AllGenes$CellType == Groups[4]]
ImmatureMyeloid <- AllGenes$Genes[AllGenes$CellType == Groups[5]]
PreBCell <- AllGenes$Genes[AllGenes$CellType == Groups[6]]
ImmatureMyeloid2 <- AllGenes$Genes[AllGenes$CellType == Groups[7]]
BM_combined <- AddModuleScore(BM_combined, features = list(MegEryth1), name = "MegEryth_1_")
MegEryth_plot <- FeaturePlot(BM_combined, "MegEryth_1_1", min.cutoff = 'q20') + ggtitle("MegEryth1")
BM_combined <- AddModuleScore(BM_combined, features = list(MegEryth2), name = "MegEryth_2_")
MegEryth1_plot <- FeaturePlot(BM_combined, "MegEryth_2_1", min.cutoff = 'q20') + ggtitle("MegEryth2")
BM_combined <- AddModuleScore(BM_combined, features = list(MEP), name = "MEP")
MEP_plot <- FeaturePlot(BM_combined, "MEP1", min.cutoff = 'q20') + ggtitle("MEP")
BM_combined <- AddModuleScore(BM_combined, features = list(Mega), name = "Mega")
Mega_plot <- FeaturePlot(BM_combined, "Mega1", min.cutoff = 'q20') + ggtitle("Mega")
Mega_plot
BM_combined <- AddModuleScore(BM_combined, features = list(ImmatureMyeloid), name = "ImmatureMyeloid")
ImmatureMyeloid_plot <- FeaturePlot(BM_combined, "ImmatureMyeloid1", min.cutoff = 'q20') + ggtitle("ImmatureMyeloid")
BM_combined <- AddModuleScore(BM_combined, features = list(ImmatureMyeloid2), name = "ImmatureMyeloid_2")
ImmatureMyeloid2_plot <- FeaturePlot(BM_combined, "ImmatureMyeloid_21", min.cutoff = 'q20') + ggtitle("ImmatureMyeloid1")
BM_combined <- AddModuleScore(BM_combined, features = list(PreBCell), name = "PreBCell")
PreBCell_plot <- FeaturePlot(BM_combined, "PreBCell1", min.cutoff = 'q20') + ggtitle("PreBCell")