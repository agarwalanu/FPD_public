
# Script to integrate the FPD and HD samples using the Seurat objects
# created during preprocessing

## laoding libraries
library(Seurat)
library(rgl)
library(tidyverse)

## list with seurat objects
SeuratObject.List <- list(HD19002_1, HD19002_Enriched,
                          FPD1_1, FPD1_2, FPD1_3, FPD1_4,
                          FPD1_Enriched_1, FPD1_Enriched_2,
                          FPD2_Enriched_1, FPD2_Enriched_2,
                          FPD3_Enriched_1, FPD3_Enriched_2,
                          FPD4_Enriched_1, FPD4_Enriched_2,
                          FPD6_Enriched_1,  
                          FPD7_1, FPD7_Enriched_1,  FPD7_CD11BCD14_1, FPD7_CD11BCD14minus_1, 
                          FPD8_1, FPD8_Enriched_1, 
                          FPD9_1, FPD9_Enriched_1, 
                          FPD10_1, FPD10_Enriched_1, 
                          HD4_1, HD4_Enriched_1,
                          HD6_1, HD5, HD_4_2, FPD12)

referenceDatasets <- list(Healthy_BM1, Healthy_BM2,
                          HD4_1, HD4_Enriched_1,
                          FPD8_1, FPD8_Enriched_1)

featuresToUse<- SelectIntegrationFeatures(SeuratObject.List)

BM.anchors <- FindIntegrationAnchors(object.list = SeuratObject.List, 
                                     anchor.features = featuresToUse,
                                     reference = c(1,2,20,21), 
                                     dims = 1:50)

system.time(BM_combined <- IntegrateData(anchorset = BM.anchors, dims = 1:50))
all.genes <- rownames(BM_combined)
BM_combined <- ScaleData(object = BM_combined)

## Perform PCA
BM_combined <- RunPCA(object = BM_combined, npcs = 75)
ElbowPlot(BM_combined, ndims = 75) + geom_vline(xintercept=50, color = "red")

## Perform Clustering
BM_combined <- FindNeighbors(object = BM_combined, dims = 1:50, nn.eps = 0.5)
BM_combined <- FindClusters(object = BM_combined, resolution = 0.7, n.start = 10)

## UMAP
BM_combined <- RunUMAP(object = BM_combined, dims = 1:50, n.components = 2, n.neighbors = 40, min.dist = 0.3, verbose = T)

rm(BM.anchors)

BM_combined_minDist0.3.Nneighbors40_3D <- RunUMAP(object = BM_combined_minDist0.3.Nneighbors60, dims = 1:50, n.components = 3, n.neighbors = 40, min.dist = 0.4, verbose = T)

BM_combined_Normal <- readRDS("~/BM_combined_May12.RDS")
set.seed(2)
BM_combined_Normal <- subset(BM_combined_Normal, cells = sample(colnames(BM_combined_Normal), size = 25000, replace =F))
DefaultAssay(BM_combined_Normal) <- "integrated"
rm(BM.anchors.tranfer)
BM_combined_Normal <- FindVariableFeatures(BM_combined_Normal)
AnchorFeatures <- intersect(VariableFeatures(BM_combined_Normal), VariableFeatures(BM_combined))

BM_combined_Normal <- RenameCells(BM_combined_Normal, add.cell.id = "First")
BM.anchors.tranfer <- FindTransferAnchors(reference = BM_combined_Normal, 
                                          query = BM_combined, 
                                          features = AnchorFeatures,
                                          dims = 1:30)
predictions <- TransferData(anchorset = BM.anchors.tranfer, 
                            refdata = BM_combined_Normal$CelltypeProgenitor_2, 
                            dims = 1:30)
BM_combined <- AddMetaData(BM_combined, metadata = predictions$predicted.id, "predicted")
DimPlot(BM_combined, reduction = "umap", group.by = "predicted", label =T)

## Get rid of cd11b 
BM_combined <- subset(BM_combined, cells = names(grep("CD11B", BM_combined$Enrichment, value = T)), invert=T)
DimPlot(BM_combined, group.by = "Celltype")

## Add predicted IDs as Ident
Idents(BM_combined) <- "predicted"
BM_combined$Celltype <- Idents(BM_combined)
Idents(BM_combined)

## Add METADATA columns
SampleNames <- names(table(BM_combined$orig.ident))
SampleNames

#Condition
Condition <- c(rep("FPD", 24), rep("Healthy", 7))
names(Condition) <- SampleNames
BM_combined$Condition <- plyr::revalue(BM_combined$orig.ident, Condition)
FeaturePlot(BM_combined, "CD74", split.by = "Condition")
VlnPlot(BM_combined, "MIF", group.by =  "Celltype", split.by = "Condition", pt.size=0)

#Donor
Donor <- gsub("(FPD[0-9]+)_.*", "\\1", SampleNames)
Donor <- gsub("(HD[0-9]+)_.*", "\\1", Donor)
Donor <- gsub("Healthy_", "", Donor)
Donor <- gsub("HD_4_2", "HD4", Donor)
names(Donor) <- SampleNames
BM_combined$Donor <- plyr::revalue(BM_combined$orig.ident, Donor)

# Origin
Origin <- c(rep("CEDAR", 31))
names(Origin) <- SampleNames
BM_combined$Origin <- plyr::revalue(BM_combined$orig.ident, Origin)

## Celltype Condition
BM_combined$celltype.Condition <- paste(BM_combined$Celltype, BM_combined$Condition, sep = "_")

## Condition Origin
BM_combined$DataSource <- paste(BM_combined$Condition, BM_combined$Origin, sep="_")

## Enrichment
Enrichment <- rep("Unenriched", length(SampleNames))
Enrichment[grep("Enriched", SampleNames)] <- "Enriched CD34+"
Enrichment[grep("CD11BCD14_1", SampleNames)] <- "Enriched CD11B+ CD14+"
Enrichment[grep("CD11BCD14minus", SampleNames)] <- "Enriched CD11B+"
names(Enrichment) <- SampleNames
BM_combined$Enrichment <- plyr::revalue(BM_combined$orig.ident, Enrichment)

## plotting
DimPlot(Factor(BM_combined, factor = "Condition", cellsPerIdent = 40000), pt.size = 0.001, label=T, split.by = "Condition") + NoLegend()
DimPlot(BM_combined, pt.size = 0.001, split.by = "Donor", ncol = 4) + NoLegend()
DimPlot(BM_combined, pt.size = 0.001, label=T) + NoLegend() 


