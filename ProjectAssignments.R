setwd("~/Box Sync/Alex_Data/Anupriya_cellR3/Analysis/Integrated_Realign/")
###################################
####         HD19002           ####
###################################
library(Seurat)
library(tidyverse)
####### HD19002 unenriched rep 1
HD19002_1.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/SCC190416AA/SCC190416AA_HD19002_BM/filtered_feature_bc_matrix/")
HD19002_1 <- CreateSeuratObject(counts = HD19002_1.data, project = "HD19002_1", min.cells = 3, min.features = 200)
HD19002_1[["percent.mt"]] <- PercentageFeatureSet(object = HD19002_1, pattern = "^MT-")
plot1 <- FeatureScatter(object = HD19002_1, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=10, linetype="dashed")
plot2 <- FeatureScatter(object = HD19002_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_hline(yintercept = c(200, 5000), linetype="dashed")
CombinePlots(plots = list(plot1, plot2))
HD19002_1 <- subset(x = HD19002_1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
HD19002_1 <- NormalizeData(object = HD19002_1)
HD19002_1 <- FindVariableFeatures(object = HD19002_1, selection.method = "vst", nfeatures = 2000)

####### HD19002 unenriched rep 2
HD19002_Enriched.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/SCC190416AA/SCC190416AA_HD19002_BM_CD34_MSC/filtered_feature_bc_matrix/")
HD19002_Enriched <- CreateSeuratObject(counts = HD19002_Enriched.data, project = "HD19002_Enriched", min.cells = 3, min.features = 200)
HD19002_Enriched[["percent.mt"]] <- PercentageFeatureSet(object = HD19002_Enriched, pattern = "^MT-")
plot1 <- FeatureScatter(object = HD19002_Enriched, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=10, linetype="dashed")
plot2 <- FeatureScatter(object = HD19002_Enriched, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_hline(yintercept = c(200, 4700), linetype="dashed")
CombinePlots(plots = list(plot1, plot2))
HD19002_Enriched <- subset(x = HD19002_Enriched, subset = nFeature_RNA > 200 & nFeature_RNA < 4700 & percent.mt < 11)
HD19002_Enriched <- NormalizeData(object = HD19002_Enriched)
HD19002_Enriched <- FindVariableFeatures(object = HD19002_Enriched, selection.method = "vst", nfeatures = 2000)


#####################################################
####        CEL181127AA (FPD1 Unenriched)        ####
#####################################################

####### FPD1 unenriched rep 1
FPD1_1.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/CEL181127AA/CEL181127AA_TB-BM_rep_1/filtered_feature_bc_matrix/")
FPD1_1 <- CreateSeuratObject(counts = FPD1_1.data, project = "FPD1_1", min.cells = 3, min.features = 200)
FPD1_1[["percent.mt"]] <- PercentageFeatureSet(object = FPD1_1, pattern = "^MT-")
plot1 <- FeatureScatter(object = FPD1_1, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=11, linetype="dashed")
plot2 <- FeatureScatter(object = FPD1_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_hline(yintercept = c(200, 5000), linetype="dashed")
CombinePlots(plots = list(plot1, plot2))
FPD1_1 <- subset(x = FPD1_1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 11)
FPD1_1 <- NormalizeData(object = FPD1_1)
FPD1_1 <- FindVariableFeatures(object = FPD1_1, selection.method = "vst", nfeatures = 2000)

####### FPD1 unenriched rep 2
FPD1_2.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/CEL181127AA/CEL181127AA_TB-BM_rep_2/filtered_feature_bc_matrix/")
FPD1_2 <- CreateSeuratObject(counts = FPD1_2.data, project = "FPD1_2", min.cells = 3, min.features = 200)
FPD1_2[["percent.mt"]] <- PercentageFeatureSet(object = FPD1_2, pattern = "^MT-")
plot1 <- FeatureScatter(object = FPD1_2, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=12, linetype="dashed")
plot2 <- FeatureScatter(object = FPD1_2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_hline(yintercept = c(200, 5000), linetype="dashed")
CombinePlots(plots = list(plot1, plot2))
FPD1_2 <- subset(x = FPD1_2, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 12)
FPD1_2 <- NormalizeData(object = FPD1_2)
FPD1_2 <- FindVariableFeatures(object = FPD1_2, selection.method = "vst", nfeatures = 2000)

####### FPD1 unenriched rep 3
FPD1_3.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/CEL181127AA/CEL181127AA_TB-BM_rep_3/filtered_feature_bc_matrix/")
FPD1_3 <- CreateSeuratObject(counts = FPD1_3.data, project = "FPD1_3", min.cells = 3, min.features = 200)
FPD1_3[["percent.mt"]] <- PercentageFeatureSet(object = FPD1_3, pattern = "^MT-")
plot1 <- FeatureScatter(object = FPD1_3, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=15, linetype="dashed")
plot2 <- FeatureScatter(object = FPD1_3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_hline(yintercept = c(200, 5000), linetype="dashed")
CombinePlots(plots = list(plot1, plot2))
FPD1_3 <- subset(x = FPD1_3, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
FPD1_3 <- NormalizeData(object = FPD1_3)
FPD1_3 <- FindVariableFeatures(object = FPD1_3, selection.method = "vst", nfeatures = 2000)

####### FPD1 unenriched rep 4
FPD1_4.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/CEL181127AA/CEL181127AA_TB-BM_rep_4/filtered_feature_bc_matrix/")
FPD1_4 <- CreateSeuratObject(counts = FPD1_4.data, project = "FPD1_4", min.cells = 3, min.features = 200)
FPD1_4[["percent.mt"]] <- PercentageFeatureSet(object = FPD1_4, pattern = "^MT-")
plot1 <- FeatureScatter(object = FPD1_4, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=12, linetype="dashed")
plot2 <- FeatureScatter(object = FPD1_4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_hline(yintercept = c(200, 5000), linetype="dashed")
CombinePlots(plots = list(plot1, plot2))
FPD1_4 <- subset(x = FPD1_4, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 12)
FPD1_4 <- NormalizeData(object = FPD1_4)
FPD1_4 <- FindVariableFeatures(object = FPD1_4, selection.method = "vst", nfeatures = 2000)









#####################################################
#####        SCL190221AA (FPD1 Enriched)        #####
#####################################################

####### FPD1 Enriched rep 1
FPD1_Enriched_1.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/SCL181127AA/SCL181127AA_AO_BM_rep_1/filtered_feature_bc_matrix/")
FPD1_Enriched_1 <- CreateSeuratObject(counts = FPD1_Enriched_1.data, project = "FPD1_Enriched_1", min.cells = 3, min.features = 200)
FPD1_Enriched_1[["percent.mt"]] <- PercentageFeatureSet(object = FPD1_Enriched_1, pattern = "^MT-")
plot1 <- FeatureScatter(object = FPD1_Enriched_1, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=15, linetype="dashed")
plot2 <- FeatureScatter(object = FPD1_Enriched_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_hline(yintercept = c(200, 5000), linetype="dashed")
CombinePlots(plots = list(plot1, plot2))
FPD1_Enriched_1 <- subset(x = FPD1_Enriched_1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
FPD1_Enriched_1 <- NormalizeData(object = FPD1_Enriched_1)
FPD1_Enriched_1 <- FindVariableFeatures(object = FPD1_Enriched_1, selection.method = "vst", nfeatures = 2000)

####### FPD1 Enriched rep 2
FPD1_Enriched_2.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/SCL181127AA/SCL181127AA_AO_BM_rep_2/filtered_feature_bc_matrix/")
FPD1_Enriched_2 <- CreateSeuratObject(counts = FPD1_Enriched_2.data, project = "FPD1_Enriched_2", min.cells = 3, min.features = 200)
FPD1_Enriched_2[["percent.mt"]] <- PercentageFeatureSet(object = FPD1_Enriched_2, pattern = "^MT-")
plot1 <- FeatureScatter(object = FPD1_Enriched_2, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=15, linetype="dashed")
plot2 <- FeatureScatter(object = FPD1_Enriched_2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_hline(yintercept = c(200, 5000), linetype="dashed")
CombinePlots(plots = list(plot1, plot2))
FPD1_Enriched_2 <- subset(x = FPD1_Enriched_2, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
FPD1_Enriched_2 <- NormalizeData(object = FPD1_Enriched_2)
FPD1_Enriched_2 <- FindVariableFeatures(object = FPD1_Enriched_2, selection.method = "vst", nfeatures = 2000)








########################################################
#####        SCC190314AA (FPD2, FPD3, FPD4)        #####
########################################################

####### FPD2 Enriched rep 1
FPD2_Enriched_1.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/SCC190314AA/SCC190314AA_FPD_2_1/filtered_feature_bc_matrix/")
FPD2_Enriched_1 <- CreateSeuratObject(counts = FPD2_Enriched_1.data, project = "FPD2_Enriched_1", min.cells = 3, min.features = 200)
FPD2_Enriched_1[["percent.mt"]] <- PercentageFeatureSet(object = FPD2_Enriched_1, pattern = "^MT-")
plot1 <- FeatureScatter(object = FPD2_Enriched_1, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=10, linetype="dashed")
plot2 <- FeatureScatter(object = FPD2_Enriched_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_hline(yintercept = c(200, 4600), linetype="dashed")
CombinePlots(plots = list(plot1, plot2))
FPD2_Enriched_1 <- subset(x = FPD2_Enriched_1, subset = nFeature_RNA > 200 & nFeature_RNA < 4600 & percent.mt < 10)
FPD2_Enriched_1 <- NormalizeData(object = FPD2_Enriched_1)
FPD2_Enriched_1 <- FindVariableFeatures(object = FPD2_Enriched_1, selection.method = "vst", nfeatures = 2000)

####### FPD2 Enriched rep 2
FPD2_Enriched_2.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/SCC190314AA/SCC190314AA_FPD_2_2/filtered_feature_bc_matrix/")
FPD2_Enriched_2 <- CreateSeuratObject(counts = FPD2_Enriched_2.data, project = "FPD2_Enriched_2", min.cells = 3, min.features = 200)
FPD2_Enriched_2[["percent.mt"]] <- PercentageFeatureSet(object = FPD2_Enriched_2, pattern = "^MT-")
plot1 <- FeatureScatter(object = FPD2_Enriched_2, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=9, linetype="dashed")
plot2 <- FeatureScatter(object = FPD2_Enriched_2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_hline(yintercept = c(200, 4300), linetype="dashed")
CombinePlots(plots = list(plot1, plot2))
FPD2_Enriched_2 <- subset(x = FPD2_Enriched_2, subset = nFeature_RNA > 200 & nFeature_RNA < 4300 & percent.mt < 9)
FPD2_Enriched_2 <- NormalizeData(object = FPD2_Enriched_2)
FPD2_Enriched_2 <- FindVariableFeatures(object = FPD2_Enriched_2, selection.method = "vst", nfeatures = 2000)

####### FPD3 Enriched rep 1
FPD3_Enriched_1.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/SCC190314AA/SCC190314AA_FPD_3_1/filtered_feature_bc_matrix/")
FPD3_Enriched_1 <- CreateSeuratObject(counts = FPD3_Enriched_1.data, project = "FPD3_Enriched_1", min.cells = 3, min.features = 200)
FPD3_Enriched_1[["percent.mt"]] <- PercentageFeatureSet(object = FPD3_Enriched_1, pattern = "^MT-")
plot1 <- FeatureScatter(object = FPD3_Enriched_1, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=10, linetype="dashed")
plot2 <- FeatureScatter(object = FPD3_Enriched_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_hline(yintercept = c(200, 4600), linetype="dashed")
CombinePlots(plots = list(plot1, plot2))
FPD3_Enriched_1 <- subset(x = FPD3_Enriched_1, subset = nFeature_RNA > 200 & nFeature_RNA < 4600 & percent.mt < 10)
FPD3_Enriched_1 <- NormalizeData(object = FPD3_Enriched_1)
FPD3_Enriched_1 <- FindVariableFeatures(object = FPD3_Enriched_1, selection.method = "vst", nfeatures = 2000)

####### FPD3 Enriched rep 2
FPD3_Enriched_2.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/SCC190314AA/SCC190314AA_FPD_3_2/filtered_feature_bc_matrix/")
FPD3_Enriched_2 <- CreateSeuratObject(counts = FPD3_Enriched_2.data, project = "FPD3_Enriched_2", min.cells = 3, min.features = 200)
FPD3_Enriched_2[["percent.mt"]] <- PercentageFeatureSet(object = FPD3_Enriched_2, pattern = "^MT-")
plot1 <- FeatureScatter(object = FPD3_Enriched_2, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=9, linetype="dashed")
plot2 <- FeatureScatter(object = FPD3_Enriched_2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_hline(yintercept = c(200, 4300), linetype="dashed")
CombinePlots(plots = list(plot1, plot2))
FPD3_Enriched_2 <- subset(x = FPD3_Enriched_2, subset = nFeature_RNA > 200 & nFeature_RNA < 4300 & percent.mt < 9)
FPD3_Enriched_2 <- NormalizeData(object = FPD3_Enriched_2)
FPD3_Enriched_2 <- FindVariableFeatures(object = FPD3_Enriched_2, selection.method = "vst", nfeatures = 2000)

####### FPD4 Enriched rep 1
FPD4_Enriched_1.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/SCC190314AA/SCC190314AA_FPD_4_1/filtered_feature_bc_matrix/")
FPD4_Enriched_1 <- CreateSeuratObject(counts = FPD4_Enriched_1.data, project = "FPD4_Enriched_1", min.cells = 3, min.features = 200)
FPD4_Enriched_1[["percent.mt"]] <- PercentageFeatureSet(object = FPD4_Enriched_1, pattern = "^MT-")
plot1 <- FeatureScatter(object = FPD4_Enriched_1, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=10, linetype="dashed")
plot2 <- FeatureScatter(object = FPD4_Enriched_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_hline(yintercept = c(200, 4600), linetype="dashed")
CombinePlots(plots = list(plot1, plot2))
FPD4_Enriched_1 <- subset(x = FPD4_Enriched_1, subset = nFeature_RNA > 200 & nFeature_RNA < 4600 & percent.mt < 10)
FPD4_Enriched_1 <- NormalizeData(object = FPD4_Enriched_1)
FPD4_Enriched_1 <- FindVariableFeatures(object = FPD4_Enriched_1, selection.method = "vst", nfeatures = 2000)

####### FPD4 Enriched rep 2
FPD4_Enriched_2.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/SCC190314AA/SCC190314AA_FPD_4_2/filtered_feature_bc_matrix/")
FPD4_Enriched_2 <- CreateSeuratObject(counts = FPD4_Enriched_2.data, project = "FPD4_Enriched_2", min.cells = 3, min.features = 200)
FPD4_Enriched_2[["percent.mt"]] <- PercentageFeatureSet(object = FPD4_Enriched_2, pattern = "^MT-")
plot1 <- FeatureScatter(object = FPD4_Enriched_2, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=10, linetype="dashed")
plot2 <- FeatureScatter(object = FPD4_Enriched_2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_hline(yintercept = c(200, 4600), linetype="dashed")
CombinePlots(plots = list(plot1, plot2))
FPD4_Enriched_2 <- subset(x = FPD4_Enriched_2, subset = nFeature_RNA > 200 & nFeature_RNA < 4600 & percent.mt < 10)
FPD4_Enriched_2 <- NormalizeData(object = FPD4_Enriched_2)
FPD4_Enriched_2 <- FindVariableFeatures(object = FPD4_Enriched_2, selection.method = "vst", nfeatures = 2000)





####### FPD6 enriched rep 1
FPD6_Enriched_1.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/CEL190731AA/CEL190731AA_FPD6_BM_MNC_CD34_MSC/filtered_feature_bc_matrix/")
FPD6_Enriched_1 <- CreateSeuratObject(counts = FPD6_Enriched_1.data, project = "FPD6_Enriched_1", min.cells = 3, min.features = 200)
FPD6_Enriched_1[["percent.mt"]] <- PercentageFeatureSet(object = FPD6_Enriched_1, pattern = "^MT-")
plot1 <- FeatureScatter(object = FPD6_Enriched_1, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=c(1,14), linetype="dashed")
plot2 <- FeatureScatter(object = FPD6_Enriched_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_hline(yintercept = c(200, 6200), linetype="dashed")
CombinePlots(plots = list(plot1, plot2))
FPD6_Enriched_1 <- subset(x = FPD6_Enriched_1, subset = nFeature_RNA > 200 & nFeature_RNA < 6200 & percent.mt < 14 & percent.mt >1)
FPD6_Enriched_1 <- NormalizeData(object = FPD6_Enriched_1)
FPD6_Enriched_1 <- FindVariableFeatures(object = FPD6_Enriched_1, selection.method = "vst", nfeatures = 2000)

####### FPD7 unenriched rep 1
FPD7_1.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/CEL190823AA/CEL190823AA_FPD7_MNC/filtered_feature_bc_matrix/")
FPD7_1 <- CreateSeuratObject(counts = FPD7_1.data, project = "FPD7_1", min.cells = 3, min.features = 200)
FPD7_1[["percent.mt"]] <- PercentageFeatureSet(object = FPD7_1, pattern = "^MT-")
plot1 <- FeatureScatter(object = FPD7_1, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=c(20), linetype="dashed")
plot2 <- FeatureScatter(object = FPD7_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_hline(yintercept = c(200, 6200), linetype="dashed")
CombinePlots(plots = list(plot1, plot2))
FPD7_1 <- subset(x = FPD7_1, subset = nFeature_RNA > 200 & nFeature_RNA < 6200 & percent.mt < 20)
FPD7_1 <- NormalizeData(object = FPD7_1)
FPD7_1 <- FindVariableFeatures(object = FPD7_1, selection.method = "vst", nfeatures = 2000)

####### FPD7 enriched rep 1
FPD7_Enriched_1.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/CEL190823AA/CEL190823AA_FPD7_MNC_plus_CD34_plus_MSC/filtered_feature_bc_matrix/")
FPD7_Enriched_1 <- CreateSeuratObject(counts = FPD7_Enriched_1.data, project = "FPD7_Enriched_1", min.cells = 3, min.features = 200)
FPD7_Enriched_1[["percent.mt"]] <- PercentageFeatureSet(object = FPD7_Enriched_1, pattern = "^MT-")
plot1 <- FeatureScatter(object = FPD7_Enriched_1, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=c(20), linetype="dashed")
plot2 <- FeatureScatter(object = FPD7_Enriched_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_hline(yintercept = c(200, 6200), linetype="dashed")
CombinePlots(plots = list(plot1, plot2))
FPD7_Enriched_1 <- subset(x = FPD7_Enriched_1, subset = nFeature_RNA > 200 & nFeature_RNA < 6200 & percent.mt < 20)
FPD7_Enriched_1 <- NormalizeData(object = FPD7_Enriched_1)
FPD7_Enriched_1 <- FindVariableFeatures(object = FPD7_Enriched_1, selection.method = "vst", nfeatures = 2000)

####### FPD7 unenriched rep 1
FPD7_CD11BCD14_1.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/CEL190823AA/CEL190823AA_FPD7_Cd11b_plus_CD14_plus/filtered_feature_bc_matrix/")
FPD7_CD11BCD14_1 <- CreateSeuratObject(counts = FPD7_CD11BCD14_1.data, project = "FPD7_CD11BCD14_1", min.cells = 3, min.features = 200)
FPD7_CD11BCD14_1[["percent.mt"]] <- PercentageFeatureSet(object = FPD7_CD11BCD14_1, pattern = "^MT-")
plot1 <- FeatureScatter(object = FPD7_CD11BCD14_1, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=c(20), linetype="dashed")
plot2 <- FeatureScatter(object = FPD7_CD11BCD14_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_hline(yintercept = c(200, 3600), linetype="dashed")
CombinePlots(plots = list(plot1, plot2))
FPD7_CD11BCD14_1 <- subset(x = FPD7_CD11BCD14_1, subset = nFeature_RNA > 200 & nFeature_RNA < 6200 & percent.mt < 20)
FPD7_CD11BCD14_1 <- NormalizeData(object = FPD7_CD11BCD14_1)
FPD7_CD11BCD14_1 <- FindVariableFeatures(object = FPD7_CD11BCD14_1, selection.method = "vst", nfeatures = 2000)

####### FPD7 enriched rep 1
FPD7_CD11BCD14minus_1.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/CEL190823AA/CEL190823AA_FPD7_Cd11b_plus_CD14_minus/filtered_feature_bc_matrix/")
FPD7_CD11BCD14minus_1 <- CreateSeuratObject(counts = FPD7_CD11BCD14minus_1.data, project = "FPD7_CD11BCD14minus_1", min.cells = 3, min.features = 200)
FPD7_CD11BCD14minus_1[["percent.mt"]] <- PercentageFeatureSet(object = FPD7_CD11BCD14minus_1, pattern = "^MT-")
plot1 <- FeatureScatter(object = FPD7_CD11BCD14minus_1, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=c(20), linetype="dashed") + NoLegend()
plot2 <- FeatureScatter(object = FPD7_CD11BCD14minus_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_hline(yintercept = c(200, 4000), linetype="dashed") + NoLegend()
CombinePlots(plots = list(plot1, plot2))
FPD7_CD11BCD14minus_1 <- subset(x = FPD7_CD11BCD14minus_1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
FPD7_CD11BCD14minus_1 <- NormalizeData(object = FPD7_CD11BCD14minus_1)
FPD7_CD11BCD14minus_1 <- FindVariableFeatures(object = FPD7_CD11BCD14minus_1, selection.method = "vst", nfeatures = 2000)

# ####### FPD7 unenriched rep 1
# HD3_1.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/CEL190823AA/CEL190823AA_Batch3_v2_HD3_1/filtered_feature_bc_matrix/")
# HD3_1 <- CreateSeuratObject(counts = HD3_1.data, project = "HD3_1", min.cells = 3, min.features = 200)
# HD3_1[["percent.mt"]] <- PercentageFeatureSet(object = HD3_1, pattern = "^MT-")
# plot1 <- FeatureScatter(object = HD3_1, feature1 = "nCount_RNA", feature2 = "percent.mt") +
#   geom_hline(yintercept=c(9), linetype="dashed")
# plot2 <- FeatureScatter(object = HD3_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
#   geom_hline(yintercept = c(200, 3000), linetype="dashed")
# CombinePlots(plots = list(plot1, plot2))
# HD3_1 <- subset(x = HD3_1, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 9)
# HD3_1 <- NormalizeData(object = HD3_1)
# HD3_1 <- FindVariableFeatures(object = HD3_1, selection.method = "vst", nfeatures = 2000)
# 
# ####### FPD7 unenriched rep 1
# HD3_2.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/CEL190823AA/CEL190823AA_Batch3_v2_HD3_2/filtered_feature_bc_matrix/")
# HD3_2 <- CreateSeuratObject(counts = HD3_2.data, project = "HD3_2", min.cells = 3, min.features = 200)
# HD3_2[["percent.mt"]] <- PercentageFeatureSet(object = HD3_2, pattern = "^MT-")
# plot1 <- FeatureScatter(object = HD3_2, feature1 = "nCount_RNA", feature2 = "percent.mt") +
#   geom_hline(yintercept=c(1,9), linetype="dashed")
# plot2 <- FeatureScatter(object = HD3_2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
#   geom_hline(yintercept = c(200, 1700), linetype="dashed")
# CombinePlots(plots = list(plot1, plot2))
# HD3_2 <- subset(x = HD3_2, subset = nFeature_RNA > 200 & nFeature_RNA < 1700 & percent.mt < 9 & percent.mt > 1)
# HD3_2 <- NormalizeData(object = HD3_2)
# HD3_2 <- FindVariableFeatures(object = HD3_2, selection.method = "vst", nfeatures = 2000)
# 
# ####### FPD7 unenriched rep 1
# HD3_3.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/CEL190823AA/CEL190823AA_Batch3_v2_HD3_3/filtered_feature_bc_matrix/")
# HD3_3 <- CreateSeuratObject(counts = HD3_3.data, project = "HD3_3", min.cells = 3, min.features = 200)
# HD3_3[["percent.mt"]] <- PercentageFeatureSet(object = HD3_3, pattern = "^MT-")
# plot1 <- FeatureScatter(object = HD3_3, feature1 = "nCount_RNA", feature2 = "percent.mt") +
#   geom_hline(yintercept=c(1,9), linetype="dashed")
# plot2 <- FeatureScatter(object = HD3_3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
#   geom_hline(yintercept = c(200, 800), linetype="dashed")
# CombinePlots(plots = list(plot1, plot2))
# HD3_3 <- subset(x = HD3_3, subset = nFeature_RNA > 200 & nFeature_RNA < 1700 & percent.mt < 9 & percent.mt > 1)
# HD3_3 <- NormalizeData(object = HD3_3)
# HD3_3 <- FindVariableFeatures(object = HD3_3, selection.method = "vst", nfeatures = 2000)
# 
# ####### FPD7 unenriched rep 1
# HD3_4.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/CEL190823AA/CEL190823AA_Batch3_v2_HD3_4/filtered_feature_bc_matrix/")
# HD3_4 <- CreateSeuratObject(counts = HD3_4.data, project = "HD3_4", min.cells = 3, min.features = 200)
# HD3_4[["percent.mt"]] <- PercentageFeatureSet(object = HD3_4, pattern = "^MT-")
# plot1 <- FeatureScatter(object = HD3_4, feature1 = "nCount_RNA", feature2 = "percent.mt") +
#   geom_hline(yintercept=c(1,7), linetype="dashed")
# plot2 <- FeatureScatter(object = HD3_4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
#   geom_hline(yintercept = c(200, 1250), linetype="dashed")
# CombinePlots(plots = list(plot1, plot2))
# HD3_4 <- subset(x = HD3_4, subset = nFeature_RNA > 200 & nFeature_RNA < 1200 & percent.mt < 9 & percent.mt > 1)
# HD3_4 <- NormalizeData(object = HD3_4)
# HD3_4 <- FindVariableFeatures(object = HD3_4, selection.method = "vst", nfeatures = 2000)

####### FPD7 unenriched rep 1
FPD8_1.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/CEL190905AA/CEL190905AA_FPD8_MNC/filtered_feature_bc_matrix/")
FPD8_1 <- CreateSeuratObject(counts = FPD8_1.data, project = "FPD8_1", min.cells = 3, min.features = 200)
FPD8_1[["percent.mt"]] <- PercentageFeatureSet(object = FPD8_1, pattern = "^MT-")
plot1 <- FeatureScatter(object = FPD8_1, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=c(1,15), linetype="dashed")
plot2 <- FeatureScatter(object = FPD8_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_hline(yintercept = c(200, 6000), linetype="dashed")
CombinePlots(plots = list(plot1, plot2))
FPD8_1 <- subset(x = FPD8_1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15 & percent.mt > 1)
FPD8_1 <- NormalizeData(object = FPD8_1)
FPD8_1 <- FindVariableFeatures(object = FPD8_1, selection.method = "vst", nfeatures = 2000)

####### FPD7 unenriched rep 1
FPD8_Enriched_1.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/CEL190905AA/CEL190905AA_FPD8_MNC_plus_Cd34_plus_MSC/filtered_feature_bc_matrix/")
FPD8_Enriched_1 <- CreateSeuratObject(counts = FPD8_Enriched_1.data, project = "FPD8_Enriched_1", min.cells = 3, min.features = 200)
FPD8_Enriched_1[["percent.mt"]] <- PercentageFeatureSet(object = FPD8_Enriched_1, pattern = "^MT-")
plot1 <- FeatureScatter(object = FPD8_Enriched_1, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=c(15), linetype="dashed")
plot2 <- FeatureScatter(object = FPD8_Enriched_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_hline(yintercept = c(200, 5000), linetype="dashed")
CombinePlots(plots = list(plot1, plot2))
FPD8_Enriched_1 <- subset(x = FPD8_Enriched_1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
FPD8_Enriched_1 <- NormalizeData(object = FPD8_Enriched_1)
FPD8_Enriched_1 <- FindVariableFeatures(object = FPD8_Enriched_1, selection.method = "vst", nfeatures = 2000)

####### FPD7 unenriched rep 1
HD4_1.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/CEL190905AA/CEL190905AA_HD19-004_MNC/filtered_feature_bc_matrix/")
HD4_1 <- CreateSeuratObject(counts = HD4_1.data, project = "HD4_1", min.cells = 3, min.features = 200)
HD4_1[["percent.mt"]] <- PercentageFeatureSet(object = HD4_1, pattern = "^MT-")
plot1 <- FeatureScatter(object = HD4_1, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=c(12), linetype="dashed")
plot2 <- FeatureScatter(object = HD4_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_hline(yintercept = c(200, 6000), linetype="dashed")
CombinePlots(plots = list(plot1, plot2))
HD4_1 <- subset(x = HD4_1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 12)
HD4_1 <- NormalizeData(object = HD4_1)
HD4_1 <- FindVariableFeatures(object = HD4_1, selection.method = "vst", nfeatures = 2000)

####### FPD7 unenriched rep 1
HD4_Enriched_1.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/CEL190905AA/CEL190905AA_HD19-004_MNC_plus_Cd34_plus_MSC/filtered_feature_bc_matrix/")
HD4_Enriched_1 <- CreateSeuratObject(counts = HD4_Enriched_1.data, project = "HD4_Enriched_1", min.cells = 3, min.features = 200)
HD4_Enriched_1[["percent.mt"]] <- PercentageFeatureSet(object = HD4_Enriched_1, pattern = "^MT-")
plot1 <- FeatureScatter(object = HD4_Enriched_1, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=c(15), linetype="dashed")
plot2 <- FeatureScatter(object = HD4_Enriched_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_hline(yintercept = c(200, 5400), linetype="dashed")
CombinePlots(plots = list(plot1, plot2))
HD4_Enriched_1 <- subset(x = HD4_Enriched_1, subset = nFeature_RNA > 200 & nFeature_RNA < 5400 & percent.mt < 15)
HD4_Enriched_1 <- NormalizeData(object = HD4_Enriched_1)
HD4_Enriched_1 <- FindVariableFeatures(object = HD4_Enriched_1, selection.method = "vst", nfeatures = 2000)

####### FPD7 unenriched rep 1
FPD9_1.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/CEL190911AA/CEL190911AA_FPD9/filtered_feature_bc_matrix/")
FPD9_1 <- CreateSeuratObject(counts = FPD9_1.data, project = "FPD9_1", min.cells = 3, min.features = 200)
FPD9_1[["percent.mt"]] <- PercentageFeatureSet(object = FPD9_1, pattern = "^MT-")
plot1 <- FeatureScatter(object = FPD9_1, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=c(12), linetype="dashed")
plot2 <- FeatureScatter(object = FPD9_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_hline(yintercept = c(200, 4700), linetype="dashed")
CombinePlots(plots = list(plot1, plot2))
FPD9_1 <- subset(x = FPD9_1, subset = nFeature_RNA > 200 & nFeature_RNA < 4700 & percent.mt < 12)
FPD9_1 <- NormalizeData(object = FPD9_1)
FPD9_1 <- FindVariableFeatures(object = FPD9_1, selection.method = "vst", nfeatures = 2000)

####### FPD7 unenriched rep 1
FPD9_Enriched_1.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/CEL190911AA/CEL190911AA_FPD9_CD34plus_MSC/filtered_feature_bc_matrix/")
FPD9_Enriched_1 <- CreateSeuratObject(counts = FPD9_Enriched_1.data, project = "FPD9_Enriched_1", min.cells = 3, min.features = 200)
FPD9_Enriched_1[["percent.mt"]] <- PercentageFeatureSet(object = FPD9_Enriched_1, pattern = "^MT-")
plot1 <- FeatureScatter(object = FPD9_Enriched_1, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=c(15), linetype="dashed")
plot2 <- FeatureScatter(object = FPD9_Enriched_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_hline(yintercept = c(200, 4500), linetype="dashed")
CombinePlots(plots = list(plot1, plot2))
FPD9_Enriched_1 <- subset(x = FPD9_Enriched_1, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 15)
FPD9_Enriched_1 <- NormalizeData(object = FPD9_Enriched_1)
FPD9_Enriched_1 <- FindVariableFeatures(object = FPD9_Enriched_1, selection.method = "vst", nfeatures = 2000)

####### FPD7 unenriched rep 1
FPD10_1.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/CEL190911AA/CEL190911AA_FPD10/filtered_feature_bc_matrix/")
FPD10_1 <- CreateSeuratObject(counts = FPD10_1.data, project = "FPD10_1", min.cells = 3, min.features = 200)
FPD10_1[["percent.mt"]] <- PercentageFeatureSet(object = FPD10_1, pattern = "^MT-")
plot1 <- FeatureScatter(object = FPD10_1, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=c(12), linetype="dashed")
plot2 <- FeatureScatter(object = FPD10_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_hline(yintercept = c(200, 6000), linetype="dashed")
CombinePlots(plots = list(plot1, plot2))
FPD10_1 <- subset(x = FPD10_1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 12)
FPD10_1 <- NormalizeData(object = FPD10_1)
FPD10_1 <- FindVariableFeatures(object = FPD10_1, selection.method = "vst", nfeatures = 2000)

####### FPD7 unenriched rep 1
FPD10_Enriched_1.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/CEL190911AA/CEL190911AA_FPD10_CD34plus_MSC/filtered_feature_bc_matrix/")
FPD10_Enriched_1 <- CreateSeuratObject(counts = FPD10_Enriched_1.data, project = "FPD10_Enriched_1", min.cells = 3, min.features = 200)
FPD10_Enriched_1[["percent.mt"]] <- PercentageFeatureSet(object = FPD10_Enriched_1, pattern = "^MT-")
plot1 <- FeatureScatter(object = FPD10_Enriched_1, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=c(1,15), linetype="dashed")
plot2 <- FeatureScatter(object = FPD10_Enriched_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_hline(yintercept = c(200, 5800), linetype="dashed")
CombinePlots(plots = list(plot1, plot2))
FPD10_Enriched_1 <- subset(x = FPD10_Enriched_1, subset = nFeature_RNA > 200 & nFeature_RNA < 5800 & percent.mt < 15 & percent.mt > 1)
FPD10_Enriched_1 <- NormalizeData(object = FPD10_Enriched_1)
FPD10_Enriched_1 <- FindVariableFeatures(object = FPD10_Enriched_1, selection.method = "vst", nfeatures = 2000)


###### FPD7 unenriched rep 1
HD_4_2.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/CEL190823AA/CEL190823AA_Batch_3_HD_4_v2/filtered_feature_bc_matrix/")
HD_4_2 <- CreateSeuratObject(counts = HD_4_2.data, project = "HD_4_2", min.cells = 3, min.features = 200)
HD_4_2[["percent.mt"]] <- PercentageFeatureSet(object = HD_4_2, pattern = "^MT-")
plot1 <- FeatureScatter(object = HD_4_2, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=c(10), linetype="dashed")
plot2 <- FeatureScatter(object = HD_4_2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_hline(yintercept = c(200, 3300), linetype="dashed")
CombinePlots(plots = list(plot1, plot2))
HD_4_2 <- subset(x = HD_4_2, subset = nFeature_RNA > 200 & nFeature_RNA < 3300 & percent.mt < 10)
HD_4_2 <- NormalizeData(object = HD_4_2)
HD_4_2 <- FindVariableFeatures(object = HD_4_2, selection.method = "vst", nfeatures = 2000)

####### FPD7 unenriched rep 1
HD5.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/CEL190916AA/CEL190916AA_HD19005/filtered_feature_bc_matrix/")
HD5 <- CreateSeuratObject(counts = HD5.data, project = "HD5", min.cells = 3, min.features = 200)
HD5[["percent.mt"]] <- PercentageFeatureSet(object = HD5, pattern = "^MT-")
plot1 <- FeatureScatter(object = HD5, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=c(20), linetype="dashed")
plot2 <- FeatureScatter(object = HD5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_hline(yintercept = c(200, 3800), linetype="dashed")
CombinePlots(plots = list(plot1, plot2))
HD5 <- subset(x = HD5, subset = nFeature_RNA > 200 & nFeature_RNA < 3800 & percent.mt < 20)
HD5 <- NormalizeData(object = HD5)
HD5 <- FindVariableFeatures(object = HD5, selection.method = "vst", nfeatures = 2000)

####### FPD7 unenriched rep 1
FPD12.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/CEL191022AA/CEL191022AA_FPD12_MNC/filtered_feature_bc_matrix/")
FPD12 <- CreateSeuratObject(counts = FPD12.data, project = "FPD12", min.cells = 3, min.features = 200)
FPD12[["percent.mt"]] <- PercentageFeatureSet(object = FPD12, pattern = "^MT-")
plot1 <- FeatureScatter(object = FPD12, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=c(16), linetype="dashed")
plot2 <- FeatureScatter(object = FPD12, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_hline(yintercept = c(200, 4000), linetype="dashed")
CombinePlots(plots = list(plot1, plot2))
FPD12 <- subset(x = FPD12, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 16)
FPD12 <- NormalizeData(object = FPD12)
FPD12 <- FindVariableFeatures(object = FPD12, selection.method = "vst", nfeatures = 2000)

####### FPD7 unenriched rep 1
HD6_1.data <- Read10X(data.dir = "~/Box Sync/Alex_Data/Anupriya_cellR3/CEL190911AA/CEL190911AA_FPD10_CD34plus_MSC/filtered_feature_bc_matrix/")
HD6_1 <- CreateSeuratObject(counts = HD6_1.data, project = "HD6_1", min.cells = 3, min.features = 200)
HD6_1[["percent.mt"]] <- PercentageFeatureSet(object = HD6_1, pattern = "^MT-")
plot1 <- FeatureScatter(object = HD6_1, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=c(1,10), linetype="dashed")
plot2 <- FeatureScatter(object = HD6_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_hline(yintercept = c(200, 5000), linetype="dashed")
CombinePlots(plots = list(plot1, plot2))
HD6_1 <- subset(x = HD6_1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10 & percent.mt > 1)
HD6_1 <- NormalizeData(object = HD6_1)
HD6_1 <- FindVariableFeatures(object = HD6_1, selection.method = "vst", nfeatures = 2000)



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
