
# Script to perform differential expression considering the cluster and
# the condition of the samples (FPD or HD)

# load libraries
library(SeuratData)
library(tidyverse)

equalSampleSizeDownsample <- function(Object, downsampleFactor, downsampleCellsNumber, seed = 8) {
  set.seed(seed)
  
  ## Cells To Subset0
  cellBarcodes <- c()
  
  ## DE cells 
  downsampleFactors <- Object[[downsampleFactor, drop=TRUE]]
  downsampleFactorsUniq <- unique(downsampleFactors)
  
  for (i in 1:length(downsampleFactorsUniq)) {
    downsampleCellsByFactor <- names(grep(paste0("^", downsampleFactorsUniq[i], "$"), downsampleFactors, value = T))
    
    ## Intersect the barcodes to actually sample
    barcodesToDownsample <- downsampleCellsByFactor
    
    ## Skip samples that have no cells of that type 
    if (length(barcodesToDownsample) > 0) {
      if(length(barcodesToDownsample) > downsampleCellsNumber) {
        sampledBarcodes<- sample(barcodesToDownsample, size = downsampleCellsNumber)
        cellBarcodes <- c(cellBarcodes, sampledBarcodes)
      } else {
        cellBarcodes <- c(cellBarcodes, barcodesToDownsample)
      }
    }
  }
  Object <- Seurat::SubsetData(Object, cells = cellBarcodes)
  return(Object)
}
BM_combined_Unenriched <- equalSampleSizeDownsample(BM_combined_Unenriched, downsampleFactor = "Condition", downsampleCellsNumber = 31450)

table(BM_combined_Unenriched$Condition)
equalSampleSizeDownsample <- function(Object, downsampleFactor, DEFoctorsOfInterest, DEGroupBy, downsampleCellsNumber, seed = 8) {
  set.seed(seed)
  
  ## Cells To Subset
  cellBarcodes <- c()
  
  ## DE cells of interest
  DEFactors  <- Object[[DEGroupBy, drop=TRUE]]
  DECellsOfInterest <- names(grep(paste(paste0("^", DEFoctorsOfInterest, "$"), collapse = "|"), DEFactors, value = T))
  #print(DECellsOfInterest)
  ## DE cells 
  downsampleFactors <- Object[[downsampleFactor, drop=TRUE]]
  downsampleFactorsUniq <- unique(downsampleFactors)
  
  for (i in 1:length(downsampleFactorsUniq)) {
    downsampleCellsByFactor <- names(grep(paste0("^", downsampleFactorsUniq[i], "$"), downsampleFactors, value = T))
    
    ## Intersect the barcodes to actually sample
    barcodesToDownsample <- intersect(downsampleCellsByFactor, DECellsOfInterest)
    
    ## Skip samples that have no cells of that type 
    if (length(barcodesToDownsample) > 0) {
      if(length(barcodesToDownsample) > downsampleCellsNumber) {
        sampledBarcodes<- sample(barcodesToDownsample, size = downsampleCellsNumber)
        cellBarcodes <- c(cellBarcodes, sampledBarcodes)
      } else {
        cellBarcodes <- c(cellBarcodes, barcodesToDownsample)
      }
    }
  }
  Object <- Seurat::SubsetData(Object, cells = cellBarcodes)
  return(Object)
}


DifferentialExpression <- function(Object, ident.1, ident.2, DE.group.by, outFolder = "DE", maxCellsIdent = Inf, plotCellProportions = FALSE) {
  if (!dir.exists(outFolder)) {
    system(sprintf("mkdir %s", outFolder))
  }
  if (maxCellsIdent != Inf) {
    message(paste("Subsetting Idents to Equal Sizes By",  DE.group.by))
    Object <- downsampleIdentByFactor(BM_combined, factor = DE.group.by, cellsPerIdent = maxCellsIdent)
  }
  
  
  ## DE by group
  message(sprintf("Calculating DE of %s vs %s", ident.1, ident.2))
  Idents(Object) <- DE.group.by
  tmpDETable <- FindMarkers(Object, 
                            ident.1 = ident.1, 
                            ident.2 = ident.2)#,
  # logfc.threshold = -Inf,
  # min.pct = -Inf,
  # min.diff.pct = -Inf )
  # tmpDETable <- tmpDETable[(tmpDETable[,5] < 0.01),]
  outFile <- paste("DE", ident.1, ident.2, "Markers.txt", sep=".")
  assign(outFile, tmpDETable, envir = .GlobalEnv)
  write.table(x = tmpDETable, 
              file = sprintf("%s/%s", outFolder, outFile), 
              quote = F,
              col.names = NA, 
              sep = "\t")
  
  ## Get avg exp for scatterPlot ident.1
  cellsOfInterest <- subset(Object, idents = c(ident.1,ident.2))
  avgExp.Cells <- log1p(AverageExpression(cellsOfInterest, verbose = FALSE)$RNA)
  avgExp.Cells$Gene <- rownames(avgExp.Cells)
  colnames(avgExp.Cells)[1:2] <- levels(cellsOfInterest)
  avgExp.Cells <- avgExp.Cells[,c(ident.1, ident.2, "Gene")]
  
  ## Plot Avg Expression
  xaxis <- sym(ident.1)
  yaxis <- sym(ident.2)
  p1 <- ggplot(avgExp.Cells, aes(!!xaxis, !!yaxis)) +
    geom_point() +
    ggtitle(sprintf("DE of %s vs %s", ident.1, ident.2)) +
    xlab(ident.1) +
    ylab(ident.2) +
    theme_classic()
  if (nrow(tmpDETable) > 0) {
    genes.to.label = rownames(tmpDETable)[1:30]
    p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
  }
  
  ## Save PDF and PNG
  pdf(file = sprintf("%s/%s.pdf", outFolder, outFile),
      width = 5,
      height = 5)
  plot(p1)
  dev.off()
  png(file = sprintf("%s/%s.png", outFolder, outFile),
      width = 5,
      height = 5,
      units = "in",
      res =100)
  plot(p1)
  dev.off()
  if (plotCellProportions != FALSE) {
    barGraphCellType(SeuratObject = Object,
                     Group.By = DE.group.by, 
                     PlotIdents = c(ident.1, ident.2), 
                     file = sprintf("%s/%s_cellProportions.pdf", outFolder, outFile))
  }
  return(p1)
}

barGraphCellType <- function(SeuratObject, Group.By, PlotIdents = NULL, file = NULL, Title = "Percent Celltype by Sample") {
  CellTypeDF <- as.data.frame(SeuratObject@meta.data %>% 
                                group_by_(Group.By, "Celltype") %>% 
                                tally())
  if (!is.null(PlotIdents)) {
    keep <-  CellTypeDF[,Group.By] %in% PlotIdents 
    CellTypeDF <- CellTypeDF[keep,]
  }
  library(RColorBrewer)
  set.seed(2019)
  n <- 15
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  p1 <- ggplot(as.data.frame(CellTypeDF), aes_string(x=Group.By, y="n", fill="Celltype")) + 
    geom_bar(position = "fill",stat = "identity") + 
    scale_fill_manual(values=colorRampPalette(brewer.pal(8, name = "Set1"))(18)) +
    scale_y_continuous(labels = scales::percent_format(), expand = c(0,0)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), 
          axis.ticks.x = element_blank(),
          axis.line.x=element_blank()) +
    xlab("Sample") +
    ylab("Percentage") +
    ggtitle(Title)
  if (!is.null(file)) {
    pdf(file, width = 11, height = 8.5)
    print(p1)
    dev.off()
  } else {
    print(p1) 
  }
}
DefaultAssay(BM_combined) <- "RNA"

DifferentialExpression(Object = BM_combined, ident.1 = "Stroma", ident.2 = "MKP", DE.group.by = "CellType", outFolder = "DETest", plotCellProportions = T)

## Differential Expression Across Condition
DifferentialExpression(Object = BM_combined, ident.1 = "Healthy", ident.2 = "FPD", DE.group.by = "Condition", outFolder = "DECondition", plotCellProportions = T)

possibleIdents <- unique(BM_combined$CellType)

for (i in 1:length(possibleIdents)) {
  ident1 <- paste(possibleIdents[i], "FPD", sep="_")
  ident2 <- paste(possibleIdents[i], "Healthy", sep = "_")
  message(sprintf("Calculating DE of %s vs %s", ident1, ident2))
  DifferentialExpression(Object = BM_combined,
                         ident.1 = ident1,
                         ident.2 = ident2,
                         DE.group.by = "celltype.Condition",
                         outFolder = "DE_ClusterCondition")
}

for (i in 1:length(possibleIdents)) {
  ident1 <- paste(possibleIdents[i], "FPD", sep="_")
  ident2 <- paste(possibleIdents[i], "Healthy", sep = "_")
  message(sprintf("Calculating DE of %s vs %s", ident1, ident2))
  DifferentialExpression(Object = BM_combined,
                         ident.1 = ident1,
                         ident.2 = ident2,
                         DE.group.by = "celltype.Condition",
                         outFolder = "DE_ClusterCondition_NewProgenitors")
}

equalSampleSizeDownsample <- function(Object, downsampleFactor, DEFoctorsOfInterest, DEGroupBy, downsampleCellsNumber, seed = 8) {
  set.seed(seed)
  
  ## Cells To Subset
  cellBarcodes <- c()
  
  ## DE cells of interest
  DEFactors  <- Object[[DEGroupBy, drop=TRUE]]
  DECellsOfInterest <- names(grep(paste(paste0("^", DEFoctorsOfInterest, "$"), collapse = "|"), DEFactors, value = T))
  #print(DECellsOfInterest)
  ## DE cells 
  downsampleFactors <- Object[[downsampleFactor, drop=TRUE]]
  downsampleFactorsUniq <- unique(downsampleFactors)
  
  for (i in 1:length(downsampleFactorsUniq)) {
    downsampleCellsByFactor <- names(grep(paste0("^", downsampleFactorsUniq[i], "$"), downsampleFactors, value = T))
    
    ## Intersect the barcodes to actually sample
    barcodesToDownsample <- intersect(downsampleCellsByFactor, DECellsOfInterest)
    
    ## Skip samples that have no cells of that type 
    if (length(barcodesToDownsample) > 0) {
      if(length(barcodesToDownsample) > downsampleCellsNumber) {
        sampledBarcodes<- sample(barcodesToDownsample, size = downsampleCellsNumber)
        cellBarcodes <- c(cellBarcodes, sampledBarcodes)
      } else {
        cellBarcodes <- c(cellBarcodes, barcodesToDownsample)
      }
    }
  }
  Object <- Seurat::SubsetData(Object, cells = cellBarcodes)
  return(Object)
}

plotByFactor <- function(Object, downsampleFactor, DEGroupBy, outFile = NULL ) {
  plot.df <- Object@meta.data %>% group_by_(downsampleFactor, DEGroupBy) %>% tally()
  downsampleFactor<- sym(downsampleFactor)
  DEGroupBy<- sym(DEGroupBy)
  if (length(outFile) > 0) {
    pdf(outFile, width = 6, height = 6)
    p1 <- ggplot(data.frame(plot.df), aes(x = fct_reorder(!!downsampleFactor, !!DEGroupBy), y=n,  fill = !!DEGroupBy)) + 
      geom_bar(stat = "identity", position = "dodge") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      xlab(downsampleFactor) + 
      ggtitle(paste("Number of", downsampleFactor, "Per Comparison"))
    print(p1)
    dev.off()
    print("Worked!")
  }
  ggplot(data.frame(plot.df), aes(x = fct_reorder(!!downsampleFactor, !!DEGroupBy), y=n,  fill = !!DEGroupBy)) + 
    geom_bar(stat = "identity", position = "dodge") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab(downsampleFactor) + 
    ggtitle(paste("Number of", downsampleFactor, "Per Comparison"))
}

DifferentialExpression <- function(Object, 
                                   ident.1, 
                                   ident.2, 
                                   DE.group.by, 
                                   outFolder = "DE", 
                                   maxCellsIdent = Inf, 
                                   plotCellProportions = FALSE,
                                   downsampleCellsNumber = NULL,
                                   plotCellDonorDistrubition = FALSE) {
  if (!dir.exists(outFolder)) {
    system(sprintf("mkdir %s", outFolder))
  }
  if (maxCellsIdent != Inf) {
    message(paste("Subsetting Idents to Equal Sizes By",  DE.group.by))
    Object <- downsampleIdentByFactor(BM_combined, factor = DE.group.by, cellsPerIdent = maxCellsIdent)
  }
  if (length(downsampleCellsNumber) > 0) {
    cat("subsetting Data...")
    Object <- equalSampleSizeDownsample(Object = Object, 
                                        downsampleFactor = "Donor", 
                                        DEGroupBy = DE.group.by, 
                                        DEFoctorsOfInterest = c(ident.1, ident.2),
                                        downsampleCellsNumber = downsampleCellsNumber)
  }
  
  ## DE by group
  message(sprintf("Calculating DE of %s vs %s", ident.1, ident.2))
  Idents(Object) <- DE.group.by
  tmpDETable <- FindMarkers(Object, 
                            ident.1 = ident.1, 
                            ident.2 = ident.2,
                            logfc.threshold = -Inf,
                            min.pct = -Inf,
                            min.diff.pct = -Inf)
  # tmpDETable <- tmpDETable[(tmpDETable[,5] < 0.01),]
  outFile <- paste("DE", ident.1, ident.2, "Markers.txt", sep=".")
  assign(outFile, tmpDETable, envir = .GlobalEnv)
  write.table(x = tmpDETable, 
              file = sprintf("%s/%s", outFolder, outFile), 
              quote = F,
              col.names = NA, 
              sep = "\t")
  
  ## Get avg exp for scatterPlot ident.1
  cellsOfInterest <- subset(Object, idents = c(ident.1,ident.2))
  avgExp.Cells <- log1p(AverageExpression(cellsOfInterest, verbose = FALSE)$RNA)
  avgExp.Cells$Gene <- rownames(avgExp.Cells)
  colnames(avgExp.Cells)[1:2] <- levels(cellsOfInterest)
  avgExp.Cells <- avgExp.Cells[,c(ident.1, ident.2, "Gene")]
  
  ## Plot Avg Expression
  xaxis <- sym(ident.1)
  yaxis <- sym(ident.2)
  p1 <- ggplot(avgExp.Cells, aes(!!xaxis, !!yaxis)) +
    geom_point() +
    ggtitle(sprintf("DE of %s vs %s", ident.1, ident.2)) +
    xlab(ident.1) +
    ylab(ident.2) +
    theme_classic()
  if (nrow(tmpDETable) > 0) {
    genes.to.label = rownames(tmpDETable)[1:30]
    p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
  }
  
  ## Save PDF and PNG
  pdf(file = sprintf("%s/%s.pdf", outFolder, outFile),
      width = 5,
      height = 5)
  plot(p1)
  dev.off()
  png(file = sprintf("%s/%s.png", outFolder, outFile),
      width = 5,
      height = 5,
      units = "in",
      res =100)
  plot(p1)
  dev.off()
  if (plotCellProportions != FALSE) {
    barGraphCellType(SeuratObject = Object,
                     Group.By = DE.group.by, 
                     PlotIdents = c(ident.1, ident.2), 
                     file = sprintf("%s/%s_cellProportions.pdf", outFolder, outFile))
  }
  if (plotCellDonorDistrubition == TRUE) {
    plotByFactor(Object, 
                 downsampleFactor = "Donor", 
                 DEGroupBy = DE.group.by,
                 outFile = sprintf("%s/%s_DonorProportions.pdf", outFolder, outFile))
  }
  return(p1)
}
possibleIdents <- unique(BM_combined$Celltype)
possibleIdents <- unique(BM_combined$CelltypeProgenitor_2)
BM_combined$celltype.Condition <- paste(BM_combined$CelltypeProgenitor_2, BM_combined$Condition, sep = "_")
BM_combined$celltype.Condition

for (i in 1:length(possibleIdents)) {
  ident1 <- paste(possibleIdents[i], "FPD", sep="_")
  ident2 <- paste(possibleIdents[i], "Healthy", sep = "_")
  message(sprintf("Calculating DE of %s vs %s", ident1, ident2))
  DifferentialExpression(Object = BM_combined,
                         ident.1 = ident1,
                         ident.2 = ident2,
                         DE.group.by = "celltype.Condition",
                         outFolder = "DE_ClusterConditionEqualDonor150_CellType_Progenitors_AllGenes/",
                         downsampleCellsNumber = 150,
                         plotCellDonorDistrubition = TRUE)
}
VlnPlot(equalSampleSizeDownsample(Object = BM_combined, 
                                  downsampleFactor = "Donor", 
                                  DEGroupBy = "celltype.Condition", 
                                  DEFoctorsOfInterest = c("Mono_Healthy", "Mono_FPD"),
                                  downsampleCellsNumber = 150), "DEFA3", group.by = "Donor", ncol = 4) + geom_violin(scale = "count")
for (i in 1:length(possibleIdents)) {
  ident1 <- paste(possibleIdents[i], "FPD", sep="_")
  ident2 <- paste(possibleIdents[i], "Healthy", sep = "_")
  message(sprintf("Calculating DE of %s vs %s", ident1, ident2))
  DifferentialExpression(Object = BM_combined,
                         ident.1 = ident1,
                         ident.2 = ident2,
                         DE.group.by = "celltype.Condition",
                         outFolder = "DE_ClusterConditionEqualDonor150/",
                         downsampleCellsNumber = 150,
                         plotCellDonorDistrubition = TRUE)
}

## calling the function
DifferentialExpression(Object = BM_combined, ident.1 = "HSC_FPD", ident.2 = "HSC_Healthy", DE.group.by = "celltype.Condition", outFolder = "DE_ClusterConditionEqualDonor/",  plotCellDonorDistrubition = TRUE)

DifferentialExpression(Object = BM_combined,
                       ident.1 = "MKP_FPD",
                       ident.2 = "MKP_Healthy",
                       DE.group.by = "celltype.Condition",
                       outFolder = "DE_ClusterConditionEqualDonor150/",
                       downsampleCellsNumber = 30,
                       plotCellDonorDistrubition = TRUE)
DifferentialExpression(Object = BM_combined,
                       ident.1 = "CLP_FPD",
                       ident.2 = "CLP_Healthy",
                       DE.group.by = "celltype.Condition",
                       outFolder = "DE_ClusterConditionEqualDonor150/",
                       downsampleCellsNumber = 50,
                       plotCellDonorDistrubition = TRUE)
DifferentialExpression(Object = BM_combined,
                       ident.1 = "Progenitor_FPD",
                       ident.2 = "Progenitor_Healthy",
                       DE.group.by = "celltype.Condition",
                       outFolder = "DE_ClusterConditionEqualDonor150/",
                       downsampleCellsNumber = 150,
                       plotCellDonorDistrubition = TRUE)

DifferentialExpression(Object = BM_combined,
                       ident.1 = "Early Eryth_FPD",
                       ident.2 = "Early Eryth_Healthy",
                       DE.group.by = "celltype.Condition",
                       outFolder = "DE_ClusterConditionEqualDonor150/",
                       downsampleCellsNumber = 150,
                       plotCellDonorDistrubition = TRUE)
DifferentialExpression(Object = BM_combined,
                       ident.1 = "HSC_FPD",
                       ident.2 = "HSC_Healthy",
                       DE.group.by = "celltype.Condition",
                       outFolder = "DE_ClusterConditionEqualDonor150_ALL/",
                       downsampleCellsNumber = 150,
                       plotCellDonorDistrubition = TRUE)
DifferentialExpression(Object = BM_combined,
                       ident.1 = "Progenitor_FPD",
                       ident.2 = "Progenitor_Healthy",
                       DE.group.by = "celltype.Condition",
                       outFolder = "DE_ClusterConditionEqualDonor150_ALL/",
                       downsampleCellsNumber = 150,
                       plotCellDonorDistrubition = TRUE)
DifferentialExpression(Object = BM_combined,
                       ident.1 = "FPD",
                       ident.2 = "Healthy",
                       DE.group.by = "Condition",
                       outFolder = "DE_ClusterConditionEqualDonor150_ALL/",
                       downsampleCellsNumber = 3000,
                       plotCellDonorDistrubition = TRUE)
DifferentialExpression(Object = BM_combined,
                       ident.1 = "FPD",
                       ident.2 = "Healthy",
                       DE.group.by = "Condition",
                       outFolder = "DE_ClusterConditionEqualDonor150_ALL_Max_Cells/",
                       maxCellsIdent = 150,
                       downsampleCellsNumber = 3000,
                       plotCellDonorDistrubition = TRUE)