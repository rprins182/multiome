# script to perform standard workflow steps to analyze single cell RNA-Seq data
# data: 20k Mixture of NSCLC DTCs from 7 donors, 3' v3.1
# data source: https://www.10xgenomics.com/resources/datasets/10-k-human-pbm-cs-multiome-v-1-0-chromium-controller-1-standard-2-0-0         

# setwd("~/Desktop/demo/single_cell_RNASeq/scripts")
#BiocManager::install("GenomeInfoDb")
#BiocManager::install('BSgenome')
#BiocManager::install('rtracklayer')
#BiocManager::install('BSgenome.Mmusculus.UCSC.mm39')
#BiocManager::install('BSgenome.Mmusculus.Ensembl.mm39')

#BiocManager::install('PWMEnrich.Mmusculus.background')
#BiocManager::install('PWMEnrich')
#BiocManager::install("glmGamPoi")

# install sctransform from Github
#install.packages("sctransform")
#BiocManager::install('MotifDb')
#BiocManager::install('Biostrings')
#BiocManager::install('JASPAR2020')
#install.packages("ggpubr")
#install.packages(fcoex)
#install.packages('BiocManager')
#BiocManager::install('multtest')
#install.packages('metap')
#install.packages('Signac')
#BiocManager::install("EnsDb.Mmusculul.v79")
#BiocManager::install("Rsamtools")
# if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install("biovizBase")

#BiocManager::install("EnsDb.Mmusculus.v79")'
#BiocManager::install("SingleR")
#install.packages("remotes")
#remotes::install_github("LTLA/celldex")
#devtools::install_github('exaexa/scattermore')

# load libraries
#install.packages("matrixStats")
library(Seurat)
library(tidyverse)
library(cowplot)
library(ggplot2)
library(patchwork)
library(dplyr)
library(fcoex)
library(ggpubr)
library(SingleR)
library(glmGamPoi)
library(sctransform)
library("pheatmap")

#multiome
library(Signac)
library(EnsDb.Mmusculus.v79)
library("LOLA")
library(Biostrings)
library(JASPAR2020)
library(MotifDb)
library(BSgenome.Mmusculus.UCSC.mm10)

library(rtracklayer)
library(BSgenome)
library(GenomeInfoDb)
library(PWMEnrich)
library(PWMEnrich.Mmusculus.background)
library('glmGamPoi')
library("DESeq2")

library("ggrepel")

###############################################################

#sample_file_path = '/Users/rubenprins/Library/Mobile Documents/com~apple~CloudDocs/USC stuff/IL23 Transcription Factor work/10xGenomics Multiome sequencing/HPV-2-GEX_out/2filtered_feature_bc_matrix.h5' 
#sample_atac_fragments <- "/Users/rubenprins/Library/Mobile Documents/com~apple~CloudDocs/USC stuff/IL23 Transcription Factor work/10xGenomics Multiome sequencing/HPV-2-GEX_out/atac_fragments.tsv.gz"
project_label = "HPV"
search_pattern = "^mt-"
genome = "mm39"
control_file_path = '/Users/rubenprins/Library/Mobile Documents/com~apple~CloudDocs/USC stuff/IL23 Transcription Factor work/10xGenomics Multiome sequencing/GFP-1-GEX_out/2filtered_feature_bc_matrix.h5'
control_atac_fragments <- 
# Setting reference annotations
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- genome

##########################################################

read_10x <- function (file_path){
  data_file <- Read10X_h5 (finemane = file_path,
                            use.names = TRUE,
                            unique.features = TRUE)
  return data_file
}

columns_selection <- function (s4_matrices, column_name){
    selected_column <- s4_matrices$column_name
    return selected_column
}

atac_table_settings <- function (s4_matrix){
  table_chr_coord <- StringToGRanges(rownames(s4_matrix), sep = c(":", "-"))
  table_chr_named_coord <- seqnames(table_chr_coord) %in% standardChromosomes(table_chr_coord)
  s4_matrix_atac_named_coord <- s4_matrix_atac[as.vector(table_chr_named_coord), ]
  return s4_matrix_atac_named_coord
}

chrom_assay <- function (s4_matrix_atac_named_coord, annotations, genome, atac_fragments){
  chrom_assay <- CreateChromatinAssay(
    counts = s4_matrix_atac_named_coord,
    sep = c(":", "-"),
    genome = genome,
    fragments = atac_fragments,
    min.cells = 10,
    annotation = annotations
  )
  return chrom_assay
}

add_to_seurat <- function(seurat_object, column_to_add, data_to_add = "", pattern = FALSE) {
    if (pattern) {
        data_to_add <- PercentageFeatureSet(seurat_obj, pattern = pattern)  # Assuming this function call is correct
    }
    seurat_object[[column_to_add]] <- data_to_add
    return(seurat_object)
}

##########################################################


sample_s4_matrices <- read_10x(sample_file_path)
control_s4_matrices <- read_10(control_file_path)

sample_s4_matrix_rna <- columns_selection(sample_s4_matrices, `Gene Expression`)
sample_s4_matrix_atac <- columns_selection(sample_s4_matrices, `Peaks`)

control_s4_matrix_rna <- columns_selection(control_s4_matrices, `Gene Expression`)
control_s4_matrix_atac <- columns_selection(control_s4_matrices, `Peaks`)

sample_seurat_obj <- CreateSeuratObject(counts = sample_s4_matrix_rna, project = sample_project_label, min.cells = 5)
control_seurat_obj <- CreateSeuratObject(counts = control_s4_matrix_rna, project = control_project_label, min.cells = 5)

# Setting the table
sample_s4_matrix_atac_named_coord <- atac_table_settings(sample_s4_matrix_atac)
control_s4_matrix_atac_named_coord <- atac_table_settings(control_s4_matrix_atac)

sample_chrom_assay <- (sample_s4_matrix_atac_named_coord, annotation, genome, sample_atac_fragments)
control_chrom_assay <- (control_s4_matrix_atac_named_coord, annotation, genome, control_atac_fragments)

# FILTERING
# Add new column with chrom_assay object
sample_seurat_obj <- add_to_seurat(sample_seurat_obj, "ATAC", sample_chrom_assay)
control_seurat_obj <- add_to_seurat(control_seurat_obj, "ATAC", control_chrom_assay)

sample_seurat_obj <- add_to_seurat(sample_seurat_obj, "percent.mt", pattern=search_pattern)
control_seurat_obj <- add_to_seurat(control_seurat_obj, "percent.mt", pattern=search_pattern)


#Violin plot visualization for the percentage
sample_qc_vlnplot <- VlnPlot(sample_seaurat_obj, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0)
control_qc_vlnplot <- VlnPlot(control_seaurat_obj, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0)






#Sets boundaries for the cell events to be included        
seurat_obj_rna <- subset(seurat_obj_rna, subset = nCount_ATAC >5e3 & nFeature_RNA > 1000 & nFeature_RNA < 3500 & percent.mt < 30 & nCount_RNA >2000)






ctrl.datarna <-  ctrl.data12$`Gene Expression`
ctrl.dataatac <-  ctrl.data12$`Peaks`

ctrl12 <- CreateSeuratObject(counts = ctrl.datarna, project = "GFP", min.cells = 5)

# ATAC analysis add gene annotation information
grange.counts <- StringToGRanges(rownames(ctrl.dataatac), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
ctrl.dataatac <- ctrl.dataatac[as.vector(grange.use), ]
#annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
#seqlevelsStyle(annotations) <- "UCSC"
#genome(annotations) <- "mm39"

frag.filegfp <- "/Users/rubenprins/Library/Mobile Documents/com~apple~CloudDocs/USC stuff/IL23 Transcription Factor work/10xGenomics Multiome sequencing/GFP-1-GEX_out/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = ctrl.dataatac,
  sep = c(":", "-"),
  genome = 'mm39',
  fragments = frag.filegfp,
  min.cells = 10,
  annotation = annotations
)

ctrl12[["ATAC"]] <- chrom_assay

#QC
ctrl12[["percent.mt"]] <- PercentageFeatureSet(ctrl12, pattern = "^mt-")
VlnPlot(ctrl12, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0)
ctrl12 <- subset(ctrl12, subset = nCount_ATAC >5e3 & nFeature_RNA > 1000 & nFeature_RNA < 3500 & percent.mt < 30 & nCount_RNA >2000)

# We exclude the first dimension as this is typically correlated with sequencing depth
#ctrl12 <- RunTFIDF(ctrl12)
#ctrl12 <- FindTopFeatures(ctrl12, min.cutoff = "q0")
#ctrl12 <- RunSVD(ctrl12)
#ctrl12 <- RunUMAP(ctrl12, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

#RNA analysis
# RNA analysis
#DefaultAssay(ctrl12) <- "RNA"
#ctrl12 <- SCTransform(ctrl12, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
#
#umap
#ctrl12 <- FindMultiModalNeighbors(ctrl12, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
#ctrl12 <- RunUMAP(ctrl12, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
#ctrl12 <- FindClusters(ctrl12, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

#p1 <- DimPlot(ctrl12, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
#p2 <- DimPlot(ctrl12, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
#p3 <- DimPlot(ctrl12, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
#p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

#p4 <- FeaturePlot(object = ctrl12, 
#                  features = "Il23a",
#                  cols = c("grey", "green", "red","purple"), 
#                  reduction = "wnn.umap",
#                  order = TRUE,
#                  pt.size = .5)
#p23 <- VlnPlot(object = stim21, features = "Il23a", pt.size = 1) + scale_y_continuous(limits = c(0.001,3))# + stat_compare_means(comparisons = comparisons1, label = "p.signif")

#plot_grid(p3,p4,p23)

#Il23expressionfreqGFP <- data.frame(FetchData(object = ctrl12, vars =  'Il23a'))# + stat_compare_means(comparisons = comparisons1, label = "p.signif")
#write.csv(Il23expressionfreqGFP, "/Users/rubenprins/Downloads/Il23expressionfreqgfp.csv")


#CoveragePlot(ctrl12, region = 'Il23a', features = 'Il23a', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE, ymax = 20)





# Set up stimulated object
#stim1 <- CreateSeuratObject(counts = stim.data1, project = "HPV", min.cells = 5)
#stim3 <- CreateSeuratObject(counts = stim.data3, project = "HPV", min.cells = 5)
#stim21 <- CreateSeuratObject(counts = stim.datarna, project = "HPV", min.cells = 5)
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
#stim2[["percent.mt"]] <- PercentageFeatureSet(stim2, pattern = "^mt-")
#stim21[["percent.mt"]] <- PercentageFeatureSet(stim21, pattern = "^mt-")
# Visualize QC metrics as a violin plot
#VlnPlot(stim2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#VlnPlot(stim21, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#stim2 <- merge(stim2, y = stim21, project = "HPV")

#stim1$stim <- "HPV"
#stim1 <- subset(stim1, subset = nFeature_RNA > 500)
#stim1 <- subset(stim1, subset = Il23a >0)#|Tmsb4x <5
#stim1 <- NormalizeData(stim1, verbose = FALSE)
#stim1 <- FindVariableFeatures(stim1, selection.method = "vst", nfeatures = 2000)

stim21$stim <- "HPV2"
#stim2 <- subset(stim2, subset = nFeature_RNA > 1000 & nFeature_RNA < 3500 & percent.mt < 30 & nCount_RNA >2000)

#stim2 <- subset(stim2, subset = Il23a >0)
stim21 <- NormalizeData(stim21, assay = 'RNA', verbose = FALSE)
stim21 <- NormalizeData(stim21, assay = 'ATAC', verbose = FALSE)
stim21 <- FindVariableFeatures(stim21, assay = "RNA", selection.method = "vst", nfeatures = 2000)
stim21 <- FindVariableFeatures(stim21, assay = "ATAC", selection.method = "vst", nfeatures = 2000)
#stim <- merge(stim1, y = stim2, project = "PBMC12K")# add.cell.ids = c("HPV1", "HPV2"),
#ctrl
#stim

stim21sct <- SCTransform(stim21,method = "glmGamPoi", vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE, reduction.name = 'umap.rna') %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)



#ctrl.data <- Read10X_h5(filename = '//Users/rubenprins/Downloads/GFP-1-GEX_out/filtered_feature_bc_matrix.h5', use.names = TRUE, unique.features = TRUE)
#ctrl.data11 <- Read10X_h5(filename = '//Users/rubenprins/Downloads/GFP-1-GEX_out/2filtered_feature_bc_matrix.h5', use.names = TRUE, unique.features = TRUE)
#ctrl.data11 <-  ctrl.data11$`Gene Expression`
#ctrl.data2 <- Read10X_h5(filename = '//Users/rubenprins/Downloads/GFP-2-GEX_out/filtered_feature_bc_matrix.h5', use.names = TRUE, unique.features = TRUE)


#ctrl1 <- CreateSeuratObject(counts = ctrl.data11, project = "GFP", min.cells = 5)
#ctrl11 <- CreateSeuratObject(counts = ctrl.data, project = "GFP", min.cells = 5)
#ctrl2 <- CreateSeuratObject(counts = ctrl.data2, project = "GFP", min.cells = 5)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
#ctrl1[["percent.mt"]] <- PercentageFeatureSet(ctrl1, pattern = "^mt-")
# Visualize QC metrics as a violin plot
#VlnPlot(ctrl1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


#plot1 <- FeatureScatter(ctrl1, feature1 = "nCount_RNA", feature2 = "percent.mt")
#to check limit for ncount_RNA
#plot1 + scale_x_continuous(breaks=seq(0,10000,500), waiver(), limits = c(0,10000))
#plot2 <- FeatureScatter(ctrl1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot_grid(plot1,plot2)

#ctrl1 <- merge(ctrl1, y = ctrl11, project = "HPV")


ctrl12$stim <- "GFP1"
#ctrl1 <- subset(ctrl1, subset = nFeature_RNA > 1000 & nFeature_RNA < 3500 & percent.mt < 30 & nCount_RNA >2000)
#ctrl1 <- subset(ctrl1, subset = Il23a >0)#| Epb41l3 >0|Tmsb4x <5
ctrl12 <- NormalizeData(ctrl12, assay = "RNA", verbose = FALSE)
ctrl12 <- NormalizeData(ctrl12, assay = "ATAC", verbose = FALSE)
ctrl12 <- FindVariableFeatures(ctrl12, assay = "RNA", selection.method = "vst", nfeatures = 2000)
ctrl12 <- FindVariableFeatures(ctrl12, assay = "ATAC", selection.method = "vst", nfeatures = 2000)

ctrl12sct <- SCTransform(ctrl12,method = "glmGamPoi", vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE, reduction.name = 'umap.rna') %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)

apc.list <- list(ctrl = ctrl12sct, stim = stim21sct)
features <- SelectIntegrationFeatures(object.list = apc.list, nfeatures = 3000)
apc.list <- PrepSCTIntegration(object.list = apc.list, anchor.features = features)
immune.anchors <- FindIntegrationAnchors(object.list = apc.list, normalization.method = "SCT",
                                         anchor.features = features)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")

immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30, verbose = FALSE)
immune.combined.sct <- FindNeighbors(immune.combined.sct, reduction = "pca", dims = 1:30)
immune.combined.sct <- FindClusters(immune.combined.sct, resolution = 0.3)

DimPlot(immune.combined.sct, reduction = "umap", split.by = "stim")

DefaultAssay(immune.combined.sct) <- "ATAC"
immune.combined.sct <- RunTFIDF(immune.combined.sct)
immune.combined.sct <- FindTopFeatures(immune.combined.sct, min.cutoff = 'q0')
immune.combined.sct <- RunSVD(immune.combined.sct)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

immune.combined.sct <- FindMultiModalNeighbors(immune.combined.sct, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
immune.combined.sct <- RunUMAP(immune.combined.sct, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
immune.combined.sct <- FindClusters(immune.combined.sct, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1 <- DimPlot(immune.combined.sct, reduction = "umap", label = TRUE, label.size = 5, label.box = TRUE, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(immune.combined.sct, reduction = "umap.atac",label = TRUE, label.size = 5,label.box = TRUE, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(immune.combined.sct, reduction = "wnn.umap", label = TRUE, label.size = 5,label.box = TRUE, repel = TRUE) + ggtitle("WNN")

p1 + p2 + p3 & NoLegend() & theme(axis.text=element_text(size =20), plot.title = element_text(hjust = 0.5, size = 40,face="bold"),
                                  axis.title.x = element_text(size = 30),
                                  axis.title.y = element_text(size = 30))


#ctrl2$stim <- "GFP"
#ctrl2 <- subset(ctrl2, subset = nFeature_RNA > 500)
#ctrl2 <- subset(ctrl2, subset = Il23a >0) #|Tmsb4x <5
#ctrl2 <- NormalizeData(ctrl2, verbose = FALSE)
#ctrl2 <- FindVariableFeatures(ctrl2, selection.method = "vst", nfeatures = 2000)

#ctrl <- merge(ctrl1, y = ctrl2, project = "PBMC12K")


immune.anchors <- FindIntegrationAnchors(object.list = list(ctrl12, stim21), dims = 1:10)
#to_integrate <- Reduce(intersect, lapply(immune.anchors@object.list, rownames))
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)







# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim", pt.size = .1)
p3 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
#p4 <- FeaturePlot(object = immune.combined, 
#                  features = c("Gzmb", "Siglech", "Ly6d", "Irf8", "Ccr7", "Cxcl9", "Itgax"),
#                  cols = c("grey", "green", "red","purple"), 
#                  reduction = "umap",
#                  order = TRUE,
#                  pt.size = 2)

plot_grid(p1,p2,p4)
plot_grid(p1, p2)

ref <- celldex::MouseRNAseqData()
results <- SingleR(test= as.SingleCellExperiment(immune.combined), ref=ref, labels= ref$label.main)
immune.combined$singlr_labels<-results$labels
immune.combined[[]]
# If you want to find genes making up each cell type
metadata(results)$de.genes
#visualize those genes in a heatmap
beta.markers <- unique(unlist(all.markers$Macrophages))
plotScoreHeatmap(results, order_columns_by="labels", features=beta.markers)



pR <- DimPlot(immune.combined, reduction = "umap", group.by = "singlr_labels", label = TRUE, pt.size = .1)


DefaultAssay(immune.combined2) <- 'RNA'
pall23 <- FeaturePlot(object = immune.combined, 
                  features = c("Ccl22","Il23r"), 
                  cols = c("black", "yellow"), 
                  reduction = "umap",
                  order = TRUE,
                  pt.size = 1)
pall23
plot_grid(p1, p2, pall23,p4)

DefaultAssay(immune.combined.sct) <- "SCT"
p4 <- FeaturePlot(object = immune.combined.sct, 
                  features = "Ccr2",
                  cols = c("black", "yellow"), 
                  reduction = "wnn.umap",
                  order = TRUE,
                  pt.size = 1, ncol=1) #,"Dhcr7"
p4&DarkTheme()

FeaturePlot(immune.combined.sct, features = "Il23a", split.by = "stim", max.cutoff = 3,
            cols = c("black", "yellow"), order = TRUE, reduction = "wnn.umap", pt.size = 1.5) &DarkTheme()


immune.combined$stim[immune.combined$stim == "HPV2"] <- "HPV"
immune.combined$stim[immune.combined$stim == "GFP1"] <- "GFP"
immune.combined2$stim[immune.combined2$stim == "HPV2"] <- "HPV"
immune.combined2$stim[immune.combined2$stim == "GFP1"] <- "GFP"


trackcolors <- c("green3","salmon")
#Il23 chr10-128286139-128298084
cov<-CoveragePlot(immune.combined2, region = 'Il23a', assay = 'ATAC', 
             group.by = "stim", extend.downstream = 1500, peaks = FALSE) & scale_fill_manual(values = c("green3", "salmon",
                                                                                         "#353436", "#f6f805")) & theme(text = element_text(size = 20)) 
  

markers.to.plot <- c("Il23a", "Gas6","Mertk", "Mapk14")
DotPlot(immune.combined, assay = "RNA", split.by = "stim", features = rev(markers.to.plot), cols = c("blue", "orange","red"), dot.scale = 8) + RotatedAxis()
 

Il23cells <- AverageExpression(immune.combined1, assays = "RNA")
Il23cells
write.csv(Il23cells, "/Users/rubenprins/Downloads/new_file.csv")

#mergedil23 <- merge(stim23, y=ctrl23)
#mergedil23a <- subset(mergedil23, Il23a>0)
#mergedil23r <- subset(mergedil23, Il23r>0)
comparisons1 <- list(c("GFP", "HPV"))
#p38 p35 
#p23 <- VlnPlot(object = mergedil23a, assay = "RNA",features = "Il23r", pt.size = .1) + scale_y_continuous(limits = c(0.00,8)) + stat_compare_means(comparisons = comparisons1, label = "p.signif", aes(label = paste0("p = ", after_stat(p.format))))
p23 <- VlnPlot(object = immune.combined2, assay = "RNA",features = "Il23r", pt.size = .1)# + scale_y_continuous(limits = c(0.00,8)) + stat_compare_means(comparisons = comparisons1, label = "p.signif", aes(label = paste0("p = ", after_stat(p.format))))

p23
plot_grid(p23,p2)

DefaultAssay(immune.combined.sct) <- "SCT"
FeaturePlot(immune.combined.sct, features = "Il23a", split.by = "stim", max.cutoff = 3,
            cols = c("black", "yellow"), order = TRUE, reduction = "wnn.umap")

new.cluster.ids <-c("0","1","1")
names(new.cluster.ids) <- levels(immune.combined2)
immune.combined3 <- RenameIdents(immune.combined2, new.cluster.ids)

vlnplotdata <- data.frame(FetchData(object = mergedil23, vars = c("stim","Il23a")))
write.csv(vlnplotdata, "/Users/rubenprins/Downloads/Il23adatahpvvsgfpmultiomecombsct.csv")


DefaultAssay(immune.combined.sct) <- "SCT"
f1 = FeaturePlot(immune.combined.sct,reduction = "wnn.umap", features = c("Itgax", "Itgam"),pt.size = .5, blend = TRUE,order = TRUE, cols = c("black","yellow", "blue"), blend.threshold = 0, max.cutoff = 1)
f1 & DarkTheme()
FeaturePlot(immune.combined2, features = c("Il23a", "Klf6"),pt.size = 1, blend = TRUE,order = TRUE, cols = c("white","yellow","purple"), blend.threshold = 0, ncol = 5)

c("Il23a","Rbpj","Mef2c","Junb","Maf","Mef2a","Tcf12","Cebpb","Atf3","Jun", "Fos","Creb5","Mitf", "Mafb","Jund", "Klf2","Klf7","Klf4","Atf4","Nfe2l2","Fosb","Fli1", "Egr1","Arnt","Sp3","Batf3","Srebf2","Bhlhe40","Nfil3","Zbtb18","Nfkb1","Stat5b","Stat3","Elf1","Creb1","Nr4a2","Relb","Klf6","Zeb1")

levels(immune.combined)

# Identify gene markers
all_markers <-FindAllMarkers(immune.combined,
                             min.pct =  0.65,
                             ident)

View(all_markers)

# Return top 10 markers for cluster specified 'x'
gen_marker_table <- function(x){
  all_markers[all_markers$cluster == x, ] %>%
    head(n=10)
}

# Create a data frame of results for clusters 0-6
top10_markers <- map_dfr(0:15, gen_marker_table)

View(top10_markers)

top10_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

top10_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(immune.combined.sct, features = top10$gene) 

DefaultAssay(immune.combined.sct) <- "SCT"


immune.combined <- RunTSNE(immune.combined,dims.use = 1:20,reduction.use = "pca", dim_embed = 2)



p3<-FeaturePlot(object = immune.combined, 
            features = c("Il23r", "Ccr7"), 
            cols = c("grey", "blue"), 
            reduction = "umap",
            pt.size = 2,
            order = TRUE)
p4<-FeaturePlot(object = immune.combined, 
                features = c("Il17r"),
                cols = c("grey", "green","red", "purple"), 
                reduction = "umap",
                pt.size = 2,
                order = TRUE)
p4
p5<-FeaturePlot(object = immune.combined, 
                features = c("Cd34", "Flt3", "Irf8","Kit", "Klf4"), 
                cols = c("grey", "blue"), 
                reduction = "umap",
                pt.size = 1,
                order = TRUE)

p6<-FeaturePlot(object = immune.combined, 
                features = "Pycard", 
                cols = c("grey", "blue"), 
                reduction = "umap",
                pt.size = 2,
                order = TRUE)
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")

plot_grid(p3, p2, p4, p5, p6)


FeaturePlot(object = immune.combined, 
            features = c(top10_markers[top10_markers$cluster == 12, "gene"]), 
            cols = c("grey", "blue"), 
            reduction = "umap",
            order = TRUE)


DefaultAssay(immune.combined2) <- "RNA"
nk.markers <- FindConservedMarkers(immune.combined1, ident.1 = 1, grouping.var = "stim", verbose = FALSE, exact = FALSE)
head(nk.markers)

p11 <- FeaturePlot(immune.combined1, features = c("Bst2","Ccr7","Ly6d"), min.cutoff = "q9", order = TRUE)
p12 <- DimPlot(immune.combined1, reduction = "umap", label = TRUE)
p11 + p12

data(immune.combined1)

exprs <- data.frame(GetAssayData(immune.combined1))
target <- Idents(immune.combined1)

fc <- new_fcoex(data.frame(exprs),target)
fc <- discretize(fc)
fc <- find_cbf_modules(fc,n_genes = 70, verbose = FALSE, is_parallel = FALSE)
fc <- get_nets(fc)
mod_names(fc)
network_plots <- show_net(fc)
network_plots


immune.combined2 <- subset(immune.combined1, subset = c(2,4,0,8,5))
immune.combined2 <- ScaleData(immune.combined2, verbose = FALSE)
immune.combined2 <- RunPCA(immune.combined2, npcs = 30, verbose = FALSE)

immune.combined2 <- RunUMAP(immune.combined2, reduction = "pca", dims = 1:20)
immune.combined2 <- FindNeighbors(immune.combined2, reduction = "pca", dims = 1:20)
immune.combined2 <- FindClusters(immune.combined2, resolution = 0.5)

p1 <- DimPlot(immune.combined2, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined2, reduction = "umap", label = TRUE)

p3 <- FeaturePlot(object = immune.combined2, 
                  features = c("Il23a"), 
                  cols = c("black", "yellow"), 
                  reduction = "umap",
                  order = TRUE)
p4 <- FeaturePlot(object = immune.combined2, 
                  features = c( "Il23r"), 
                  cols = c("black", "yellow"), 
                  reduction = "umap",
                  order = TRUE)
p2+p3+p4 & DarkTheme()



#immune.combined <- subset(immune.combined, idents = c("0","1","6"))

library(data.table)
data_to_write_out <- as.data.frame(as.matrix(immune.combined1@assays$RNA@data))
fwrite(x = data_to_write_out, row.names = TRUE, file = "/Users/rubenprins/Downloads/outfile.csv")

immune.combined10 <- FindVariableFeatures(immune.combined1, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(immune.combined10), 20)
plot1 <- VariableFeaturePlot(immune.combined10)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

DimHeatmap(immune.combined1, dims = 1, cells = 500, balanced = TRUE)


#Split data sets due to difference in sequencing rather than sample. Il23 subsetted data shows GFP2 and HPV1 (cluster 1) cluster based onPycard and other immune cell markers
# but 10xgenomics recommends GFP-1 and HPV-2 are better since they have high UMI count.
immune.combined1 <- subset(immune.combined2, idents = 0)#1, 0, 5, 4
immune.combined1 <- ScaleData(immune.combined1, verbose = FALSE)
immune.combined1 <- RunPCA(immune.combined1, npcs = 30, verbose = FALSE)

immune.combined1 <- RunUMAP(immune.combined1, reduction = "pca", dims = 1:20)
immune.combined1 <- FindNeighbors(immune.combined1, reduction = "pca", dims = 1:20)
immune.combined1 <- FindClusters(immune.combined1, resolution = 0.5)

p1 <- DimPlot(immune.combined1, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined1, reduction = "umap", label = TRUE)
plot_grid(p1, p2)




#Subset only Il23a cells
ctrl23 <- subset(ctrl12, subset = Il23a >0| Il23r >0)#|Tmsb4x <5
stim23 <- subset(stim21, subset = Il23a >0| Il23r >0)#|Tmsb4x <5
ctrl23 <- FindVariableFeatures(ctrl23, selection.method = "vst", nfeatures = 2000)
stim23 <- FindVariableFeatures(stim23, selection.method = "vst", nfeatures = 2000)

immune.anchors2 <- FindIntegrationAnchors(object.list = list(ctrl23, stim23), dims = 1:10, k.filter = 71)
#to_integrate <- Reduce(intersect, lapply(immune.anchors@object.list, rownames))
immune.combined2 <- IntegrateData(anchorset = immune.anchors2, dims = 1:20, k.weight = 70)
DefaultAssay(immune.combined2) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined2 <- ScaleData(immune.combined2, verbose = FALSE)
immune.combined2 <- RunPCA(immune.combined2, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
immune.combined2 <- RunUMAP(immune.combined2, reduction = "pca", dims = 1:20)
immune.combined2 <- FindNeighbors(immune.combined2, reduction = "pca", dims = 1:20)
immune.combined2 <- FindClusters(immune.combined2, resolution = 0.5)

# Visualization
p1 <- DimPlot(immune.combined2, reduction = "umap", group.by = "stim", pt.size = .4)
p22 <- DimPlot(immune.combined2, reduction = "umap", label = FALSE,pt.size = .4, cols = c("steelblue","plum"))
#upregulated
p4 <- FeaturePlot(object = immune.combined.sct, 
                  features = c("Ccl2"),
                  cols = c("grey", "green", "red","purple"), 
                  reduction = "umap",
                  order = TRUE,
                  pt.size = .4) #,"Dhcr7"
p5 <- FeaturePlot(object = immune.combined2, 
                  features = c( "Il23r"),
                  cols = c("grey", "green", "red","purple"), 
                  reduction = "umap",
                  order = TRUE,
                  pt.size = .4) #,"Dhcr7"
plot_grid(p1,p22)
plot_grid(p1, p2,p4,p5)


markers.to.plot <- "Il23a"
DotPlot(immune.combined2, assay = "RNA", features = c("Il23a","Dhcr7","Itgax",""), cols = c("blue", "orange","red"), dot.scale = 8, split.by = "stim") + RotatedAxis()
all_markers1 <-FindAllMarkers(immune.combined2,
                              assay = "ATAC",
                              min.pct =  0.25, 
                              min.diff.pct = 0.25,
                              test.use=T)

View(all_markers1)

all_markers2 <- all_markers1[all_markers1$cluster == 0, ]
view(all_markers2)
region_of_interest <- all_markers2$gene
regions1 <- region_of_interest[1:4]

cluster_list <- 0:2 # replace n with the total number of clusters in your Seurat object
markers <- FindAllMarkers(immune.combined2, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, test.use = "LR", cluster.list = cluster_list)

markers <- FindMarkers(immune.combined2, idents = "0", test.use = "LR", only.pos = TRUE, logfc.threshold = 0.25, group.by = "stim")
markers

immune.combined.sct <- PrepSCTFindMarkers(immune.combined.sct)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
t.cells <- subset(immune.combined.sct, idents = c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19"))
Idents(t.cells) <- "stim"
avg.t.cells <- log1p(AverageExpression(t.cells, verbose = FALSE)$RNA)
avg.t.cells <- as.data.frame(avg.t.cells)
p101 <- ggplot(avg.t.cells, aes(GFP1, HPV2)) + geom_point() + ggtitle("Cluster all")
p101


p44 <- VlnPlot(object = immune.combined.sct, 
               assay = "RNA",
               idents = "0",
               features = "Ccl2",
               pt.size = 0.1, group.by = "stim") + geom_boxplot()

# Show the plot
p44
library("ggpubr")
p44_with_comparison <- p44 +
  geom_signif(comparisons = list(c("GFP1", "HPV2")),
              map_signif_level = F,
              y_position = 2.5,test = t.test)  # Adjust the y_position as needed
VlnPlot(immune.combined, features = "Mafb", pt.size = .5, assay= "RNA", group.by = "stim") + scale_y_continuous(limits = c(0.1,5)) + geom_boxplot()

#Show the plot with the comparison bar
p44_with_comparison

immune.combined.1 <-immune.combined2
immune.combined.1$celltype.stim <- paste(Idents(immune.combined.1), immune.combined.1$stim, sep = "_")
immune.combined.1$celltype <- Idents(immune.combined.1)
Idents(immune.combined.1) <- "celltype.stim"
hpv.response <- FindMarkers(immune.combined.1, ident.1 = "4_HPV2", ident.2 = "0_HPV2", test.use = "DESeq2",
                             verbose = FALSE, slot = "counts", min.pct = 0.15)
hpv.response$gene <- rownames(hpv.response)
view(hpv.response)

ggplot(hpv.response, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 1e-55, gene,
                                                                   
                                                                                                                                                "")), colour = "red", size = 3)

immune.combined.1$condition <- sample(c("HPV2", "GFP1"), size = ncol(immune.combined.1), replace = TRUE)
orig.levels <- levels(immune.combined.1)
Idents(immune.combined.1) <- gsub(pattern = " ", replacement = "_", x = Idents(immune.combined.1))
orig.levels <- gsub(pattern = " ", replacement = "_", x = orig.levels)
levels(immune.combined.1) <- orig.levels
cluster.averages <- AverageExpression(immune.combined.1, return.seurat = TRUE)
DefaultAssay(cluster.averages) <- "RNA"
CellScatter(immune.combined.1, cell1 = "0_GFP", cell2 = "0_HPV", )

> FeatureScatter(immune.combined2, feature1 = "Il23a", feature2 = "Mrc1", group.by = "stim",)


pplot <- FeaturePlot(immune.combined2, features = c("Il23a", "Cd38"),blend = T, pt.size = 1, max.cutoff = 2, order = T, cols = c("gray", "red", "green"), split.by = "stim")
pplot & theme(panel.background = element_rect(fill = "gray"))
#IL23a location chr10-128296139-128298084 
CoveragePlot(immune.combined2, region = "chr15-97854425-97908310", features = 'Vdr', assay = 'ATAC', expression.assay = 'RNA', peaks = TRUE, split.by = "stim")
CoveragePlot(immune.combined, region = "chr8-72319033-72321656", idents = "4", assay = 'ATAC', expression.assay = 'RNA', peaks = TRUE, split.by = "stim")


all_markersrna <-FindAllMarkers(immune.combined2,
                              assay = "RNA",
                              min.pct =  0.25, 
                              min.diff.pct = 0.25)

View(all_markersrna)
DotPlot(immune.combined2, assay = "RNA", features = c("Il23a","Nostrin","Polr3c","Arg1","P2rx4","Camk4","Tbc1d8","Mapk14","Dlst","Il17ra","Mapk3","Mapk1","Man1a","Atf6","Rheb","Adgre1","Cd84","Itgax","Dhcr7","Ifng","Il12a","Il12b","Il1b", "Il23r","Ccr7","Rela", "Relb","Rel","Nfkb1","Nfkb2"), cols = c("blue", "orange","red"), dot.scale = 8) + RotatedAxis()
DotPlot(immune.combined2, assay = "RNA", features = c("Il23a","Adgre1","Cd84","Itgax","Dhcr7","Ifng","Il12a","Il12b","Il1b", "Il23r","Ccr7","Rela", "Relb","Rel","Nfkb1","Nfkb2"), cols = c("blue", "orange","red"), dot.scale = 8) + RotatedAxis()
DotPlot(immune.combined2, assay = "RNA", features = c("Il23a","Man1a","Gm16998","Il23r"), cols = c("blue", "orange","red"), dot.scale = 8, split.by = "stim") + RotatedAxis()


FeaturePlot(immune.combined2, features = c("Il23a", "Cd84"),pt.size = 2, blend = TRUE,order = TRUE, cols = c("black","yellow", "purple"), blend.threshold = 0)


# Return top 10 markers for cluster specified 'x'
gen_marker_table <- function(x){
  all_markers1[all_markers1$cluster == 0, ] 
}

# Create a data frame of results for clusters 0-6
top10_markers <- map_dfr(0:1, gen_marker_table)

View(top10_markers)

write.csv(top10_markers, "/Users/rubenprins/Downloads/Il23subsetpeakssignificant.csv")

ponly23 <- VlnPlot(object = immune.combined3, assay = "RNA",features = "Il23a", pt.size = .1,cols = c("steelblue","plum1")) + NoLegend() + xlab('Cluster')+ scale_y_continuous(limits = c(0.00,3))+ stat_compare_means(aes(label = "p.signif"))
 #& scale_fill_manual(values = c("green3", "salmon", "#353436", "#f6f805")) & theme(text = element_text(size = 20)) 
ponly24 <- VlnPlot(object = immune.combined3, assay = "RNA",features = "Il23r", pt.size = .1,cols = c("steelblue","plum1"))+ NoLegend() + xlab('Cluster') + scale_y_continuous(limits = c(0.00,3.5))+ stat_compare_means(aes(label = paste0("p = ")))
ponly25 <- VlnPlot(object = immune.combined2, assay = "RNA",features = "Il23r", pt.size = .1) #& scale_fill_manual(values = c("green3", "salmon","#353436", "#f6f805")) & theme(text = element_text(size = 20)) 

ponly25 <- VlnPlot(object = immune.combined2, assay = "RNA",features = c("Igf1","Mir99ahg"), split.by = "stim", pt.size = .1) & scale_fill_manual(values = c("green3", "salmon",
                                                                                                                                                             "#353436", "#f6f805")) & theme(text = element_text(size = 20)) 


plot_grid(p22,ponly23,ponly24, ncol = 3)

immune.combined3 <- subset(immune.combined2, idents = "0","1","2")

DefaultAssay(immune.combined2) <- "ATAC"
nk.markers <- FindConservedMarkers(immune.combined2, ident.1 = 0, grouping.var = "stim", verbose = FALSE, exact = FALSE)
view(nk.markers)
write.csv(nk.markers, "/Users/rubenprins/Downloads/multiome_findconservedcluster0.csv")
geniusgene<-multiome_findconservedcluster0[,1]
view(geniusgene)
DotPlot(immune.combined2, assay = "RNA", features = c(geniusgene[0:40]), cols = c("blue", "orange","red"), dot.scale = 8, split.by = "stim") + RotatedAxis()
DotPlot(immune.combined2, assay = "RNA", features = c("Cebpg","Nfatc2","Maf","Tcf12","Mafb","Il23a","Socs2","Alx1","Pou2f2","Stat2","Stat3","Cebpa","Cebpb","Iqgap2" ,"Mapk14","Stab1","Hoxb3","Arid5b","Jdp2","Pparg", "Ppargc1b","Hmbox1","Atf4"), cols = c("blue", "orange","red"), dot.scale = 8) + RotatedAxis()
#

Dmbx1, HEY1, HES5, BHLHE23, HEY2, MYBL1, HES1, CEBPG, HES2, MAF
#Chr2_1322 Pou2f2/Mafb
#
#
#
#
#
#
#
#Chr10_8769 Nfatc2/Tcf12/Mafb
#Chr11: Tcf12
#
#Chr13: Cebpb

#how to find TFs of interest
chr11 <- '/Users/rubenprins/Downloads/ATACchr11.fa'
chr11 <- readDNAStringSet(chr11)
#has mafb
chr13 <- '/Users/rubenprins/Downloads/ATACchr13.fa'
chr13 <- readDNAStringSet(chr13)
#has Mafb
chr7 <- '/Users/rubenprins/Downloads/ATACchr7.fa'
chr7 <- readDNAStringSet(chr7)

chr16 <- '/Users/rubenprins/Downloads/ATACchr16.fa'
chr16 <- readDNAStringSet(chr16)

chr2 <- '/Users/rubenprins/Downloads/ATACchr2.fa'
chr2 <- readDNAStringSet(chr2)

chr2_1322 <- '/Users/rubenprins/Downloads/ATACchr2_1322.fa'
chr2 <- readDNAStringSet(chr2_1322)


chr18 <- '/Users/rubenprins/Downloads/ATACchr18.fa'
chr18 <- readDNAStringSet(chr18)

chr10_8769 <- '/Users/rubenprins/Downloads/ATACchr10_8769.fa'
chr10 <- readDNAStringSet(chr10_8769)

chr10_9115s <- '/Users/rubenprins/Downloads/ATACchr10_9115s.fa'
chr10 <- readDNAStringSet(chr10_9115s)


chr10 <- '/Users/rubenprins/Downloads/ATACchr10il23.fa'
chr10 <- readDNAStringSet(chr10)

chr10il23prom <- '/Users/rubenprins/Downloads/ATACchr10_1282s-7.fa'
chr10il23prom <- readDNAStringSet(chr10il23prom)


chr10_1282 <- '/Users/rubenprins/Downloads/ATACchr10_1282s.fa'
chr10 <- readDNAStringSet(chr10_1282)


chr1 <- '/Users/rubenprins/Downloads/ATACchr1.fa'
chr10 <- readDNAStringSet(chr1)

chr10_12766 <- '/Users/rubenprins/Downloads/ATACchr10_12766.fa'
chr10 <- readDNAStringSet(chr10_12766)

chr10_9153s <- '/Users/rubenprins/Downloads/ATACchr10_9153s.fa'
chr10s <- readDNAStringSet(chr10_9153s)

data(PWMLogn.mm9.MotifDb.Mmus)
res = motifEnrichment(chr10, PWMLogn.mm9.MotifDb.Mmus)
report = sequenceReport(res, 1)
plot(report[report$p.value < 0.10], fontsize=7, id.fontsize=6)


# extract PWM for TFs we are interested in
ids = c("UP00093","M6046_1.02","M6038_1.02", "M6045_1.02")
ids = "M0208_1.02"
sel.pwms = PWMLogn.mm9.MotifDb.Mmus$pwms[ids]
names(sel.pwms) = c("Klf7","Mafb", "Jdp2", "Mafb200")
scores = PWMEnrich::motifScores(chr10il23prom, sel.pwms, raw.scores=TRUE, verbose=TRUE)
dim(scores[[1]])
head(scores[[1]])
plotMotifScores(scores, cols = c("orange","yellow","green","blue","red","purple"))

BiocManager::install("PWMEnrich.Hsapiens.background")
library("PWMEnrich.Hsapiens.background")

motifs.klf2 = readMotifs(system.file(package = "PWMEnrich", dir='extdata', file="MA1515.1.transfac"), remove.acc = TRUE)
genomic.acgt = getBackgroundFrequencies(("hg19"))
pwms.klf2 = toPWM(motifs.klf2, prior = genomic.acgt)
bg.klf2 = makeBackground(pwms.klf2, organism = "hg19",type = "logn", quick = TRUE)
  
motif_names <- names(PWMLogn.hg19.MotifDb.Hsap$pwms)
mdb.klf2 <- query (MotifDb, 'Jun')
view(names(mdb.klf2))
Klf2 <- mdb.klf2
ids = c("Klf2")
sel.pwms = pwms.klf2
names(sel.pwms) = c("Klf2")
scores = PWMEnrich::motifScores(chr10il23prom, sel.pwms,raw.scores = TRUE)
dim(scores[[1]])
head(scores[[1]])
plotMotifScores(scores, cols = c("orange","yellow","green","blue","red","purple"))
view(scores)


#see which tf binding site is significant

mafb.ecdf = PWMEnrich::motifEcdf(sel.pwms, "mm9", quick=TRUE)[[1]]
threshold.1e3 = log2(quantile(mafb.ecdf,1 - 1e-3))
threshold.1e3
plotMotifScores(scores, cols = "green", sel.motifs = "M0208_1.02", cutoff = threshold.1e3, seq.len.spacing = 20)
pvals = 1- mafb.ecdf(scores[[1]][,"Klf2"])
head(pvals)
write.csv(pvals, file = "/Users/rubenprins/Downloads/pvals_scores1.csv", row.names = TRUE)
view(plotMotifScores(scores, cols = "green", sel.motifs = "Klf2", cutoff = threshold.1e3, seq.len.spacing = 20))






library(BSgenome.Mmusculus.UCSC.mm10)
library("reticulate")



# call peaks using MACS2
peaks <- Signac::CallPeaks(
  immune.combined2,
  group.by = "stim",
  assay = "ATAC",
  macs2.path = "/usr/local/Caskroom/miniconda/base/envs/py27/bin/macs2"
)



# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)


DefaultAssay(immune.combined2) <- "ATAC"

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(immune.combined2),
  features = peaks,
  cells = colnames(immune.combined2)
)


# create a new assay using the MACS2 peak set and add it to the Seurat object
immune.combined2[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = Fragments(immune.combined2),
  annotation = annotations
)



DefaultAssay(immune.combined2) <- "RNA"
immune.combined2 <- SCTransform(immune.combined2)
immune.combined2 <- RunPCA(immune.combined2)

DefaultAssay(immune.combined2) <- "peaks"


# first compute the GC content for each peak
immune.combined2 <- RegionStats(immune.combined2, genome = BSgenome.Mmusculus.UCSC.mm10)

# link peaks to genes
immune.combined2 <- LinkPeaks(
  object = immune.combined2,
  peak.assay = "peaks",
  expression.assay = "SCT",
  genes.use = c("Mafb","Il23a","Cebpa","Cebpb","Mitf","Arnt","Mef2c","Maf","Tfec","Tcf12"),
  pvalue_cutoff = .05,
  distance = 1e+09,
  min.cell = 50
)

p100 <- CoveragePlot(
  object = immune.combined2,
  region = "Il23a",
  features = "Il23a",
  expression.assay = "SCT",
  extend.upstream = 10000000,
  extend.downstream = 1000000,
  split.by = "stim",
  ranges = peaks
)

p102 <- CoveragePlot(
  object = immune.combined2,
  region = "Mafb",
  features = "Mafb",
  expression.assay = "SCT",
  extend.upstream = 100,
  extend.downstream = 100,
  split.by = "stim"
)


patchwork::wrap_plots(p100, p102, ncol = 1)

p100
regionsil23a <- as.data.frame(Links(immune.combined2))
highz<- subset(regionsil23a, zscore >0)
view(highz)
lowz<- subset(regionsil23a, zscore <0)

#Alternative matif analysis
library(Signac)
library(Seurat)
library(JASPAR2022)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
#library(BSgenome.Mmusculus.UCSC.mm39)
library(patchwork)
set.seed(1234)
library(motifmatchr)
library(ggseqlogo)
library(chromVAR)
#BiocManager::install('motifmatchr')
#BiocManager::install('TFBSTools')
#BiocManager::install('BSgenome.Mmusculus.UCSC.mm10')
#BiocManager::install('chromVAR')


# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2022,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

DefaultAssay(immune.combined2) <- "peaks"
# add motif information
immune.combined2 <- AddMotifs(
  object = immune.combined2,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

da_peaks <- FindMarkers(
  object = immune.combined2,
  assay = "peaks",
  ident.1 = '0',
  ident.2 = '1',
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_ATAC'
)

# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])

# find peaks open in Pvalb or Sst cells
open.peaks <- AccessiblePeaks(immune.combined2, idents = c("0", "1"))

# match the overall GC content in the peak set
meta.feature <- GetAssayData(immune.combined2, assay = "peaks", slot = "meta.features")
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[highz$peak,],
  n = 199
)

RegionStats(immune.combined2$ATAC,
            genome = BSgenome.Mmusculus.UCSC.mm10,
            assay = 'ATAC',
            pfm = pfm
)




top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])

# Assuming you have already defined the 'top.da.peak' and 'peaks.matched' variables appropriately

# Set the default assay to 'peaks'
DefaultAssay(immune.combined2) <- "peaks"

# Perform motif enrichment analysis
enriched.motifs <- FindMotifs(
  object = immune.combined2,
  features = highz$peak,
  background = peaks.matched
)

MotifPlot(
  object = immune.combined2,
  motifs = head(rownames(enriched.motifs))
)

#immune.combined2 <- RunChromVAR(
  object = immune.combined2,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

DefaultAssay(immune.combined2) <- 'chromvar'

checkregions <- enriched.motifs[1:5,1]
# look at the activity of Mef2c
p2 <- FeaturePlot(
  object = immune.combined2,
  features = "MA1116.1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.5
)

plowzinterestdot <- DotPlot(immune.combined2, assay = "RNA", features = c("Relb","Klf6","Zeb1","Jun", "Fos", "Mafb","Jund","Il23a", "Klf2", "Egr1","Klf7","Sp3","Klf4"), cols = c("blue", "orange","red"), dot.scale = 8) + RotatedAxis()
plot_grid(p22,plowzinterestdot)
plowzinterest <- FeaturePlot(
  object = immune.combined2,
  features = c("MA1117.1","MA1517.1","MA0103.3", "MA1126.1" ,"MA0117.2","MA0492.1"), #Relb, Klf6,Zeb1, JUN::FOS, Mafb
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
plowzinterest

phighzinterestdot <- DotPlot(immune.combined2, assay = "RNA", features = c("Il23a","Sp3", "Klf2","Mef2c","Rbpj","Mitf","Mafb","Fos","Jun","Junb","Jund","Maf","Mef2a","Tcf12","Cebpb","Nfkb1","Stat5b","Stat3","Zeb1","Elf1","Creb1","Nr4a2"), cols = c("blue", "orange","red"), dot.scale = 8) + RotatedAxis()
plot_grid(p22, phighzinterestdot)
phighzinterest <- FeaturePlot(
  object = immune.combined2,
  features = c("MA0746.2","MA1515.1","MA0497.1", "MA1621.1" ,"MA0620.3","MA0117.2", "MA1951.1","MA1126.1","MA1140.2","MA0492.1","MA1520.1","MA0052.4","MA0521.2","MA1648.1","MA0466.3","MA0105.4","MA1625.1","MA0144.2","MA0103.3","MA0473.3","MA0018.4","MA0160.2"), #Sp3, Klf2, Mef2c,Rbpjl,Mitf,Mafb, Fos, FOS:Jun,Junb,Jund, Maf, Mef2a,Tcf12, TCF12,Cebpb,Nfkb1,stat5b,STAT3,Zeb1,Elf1, Creb1, Nr4a2 
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 2,
  order = TRUE
)
phighzinterest& DarkTheme()

ptopinterestdot <- DotPlot(immune.combined2, assay = "RNA", features = c("Il23a","Cebpb","Maf","Mitf","Mafb","Atf4","Mef2a","Mef2c","Sp3","Arnt","Tcf12","Junb","Klf2","Nfe2l2","Fos","Jun","Jund","Atf3","Fosb","Fli1","Creb5","Batf3","Srebf2","Elf1","Bhlhe40","Nfil3","Zbtb18","Klf6"), cols = c("blue", "orange","red"), dot.scale = 8) + RotatedAxis()
plot_grid(p22, ptopinterestdot)
ptopinterest <- FeaturePlot(
  object = immune.combined2,
  features = ,
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
ptopinterest

DefaultAssay(immune.combined2) <- 'RNA'
gene_plot <- FeaturePlot(immune.combined2, features = "Il23a")
DefaultAssay(immune.combined2) <- 'chromvar'
motif_plot <- FeaturePlot(immune.combined2, features = "MA0466.3", min.cutoff = 0, cols = c("lightgrey", "darkred"))
gene_plot | motif_plot


plot_grid(plowzinterestdot,phighzinterestdot,ptopinterestdot, nrow =3)

ptopall <- DotPlot(immune.combined2, assay = "RNA",scale.min = 75, features = c("Il23a","Mertk","Cd68","Mrc1","Lilr4b","Adgre1","Itgax","Rbpj","Mef2c","Junb","Maf","Mef2a","Tcf12","Cebpb","Atf3","Jun", "Fos","Creb5","Mitf", "Mafb","Jund", "Klf2","Klf7","Klf4","Atf4","Nfe2l2","Fosb","Fli1", "Egr1","Arnt","Sp3","Batf3","Srebf2","Bhlhe40","Nfil3","Zbtb18","Nfkb1","Stat5b","Stat3","Elf1","Creb1","Nr4a2","Relb","Klf6","Zeb1"), cols = c("blue", "orange","red"), dot.scale = 8) + RotatedAxis()
pexon2tf <- DotPlot(immune.combined2, assay = "RNA", features = c("Il23a","Cebpb","Yy1","Arid3a","Cebpe","Cebpa","Cebpd","Tcf3"), cols = c("blue", "orange","red"), dot.scale = 8) + RotatedAxis()
pexon2tf <- DotPlot(immune.combined2, assay = "RNA", features = c("Il23a","Blvrb","Bsg","Capn2","Cd276","Cd47","Cd82","Cdh13","Ckb","Ckmt1a"), cols = c("blue", "orange","red"), dot.scale = 8) + RotatedAxis()


DefaultAssay(immune.combined2) <- 'peaks'
tile_plot <- TilePlot(
  object = immune.combined2,
  region = "chr10-128296139-128300084",
  idents = c("0","1")
)

cov_plot <- CoveragePlot(
  object = immune.combined2,
  region = "chr10-128296139-128300084",
  annotation = TRUE,
  peaks = TRUE,
  links = FALSE,
  tile = TRUE,
)

gene_plot <- AnnotationPlot(
  object = immune.combined2,
  region = "chr10-128296139-128300084"
)

link_plot <- LinkPlot(
  object = immune.combined2,
  region = "chr10-128296139-128300084"
)
128297490 - 128297550
CombineTracks(plotlist = list(cov_plot, tile_plot, gene_plot, link_plot), 
              heights = c(4,2,1,0)
              ) 


new.cluster.ids <-c("0","1","1")
names(new.cluster.ids) <- levels(immune.combined2)
immune.combined4 <- RenameIdents(immune.combined2, new.cluster.ids)


DefaultAssay(immune.combined4) <- "ATAC"
pwm_set <- getMatrixSet(x = JASPAR2022, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
motif.matrix <- CreateMotifMatrix(features = granges(immune.combined4), pwm = pwm_set, genome = 'mm10', use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
immune.combined4 <- SetAssayData(immune.combined4, assay = 'ATAC', slot = 'motifs', new.data = motif.object)



#devtools::install_github("immunogenomics/presto")
library('presto')
markers_rna <- presto:::wilcoxauc.Seurat(X = immune.combined4, assay = 'data', seurat_assay = 'SCT')
markers_motifs <- presto:::wilcoxauc.Seurat(X = immune.combined4, assay = 'data', seurat_assay = 'chromvar')
motif.names <- markers_motifs$feature
colnames(markers_rna) <- paste0("RNA.", colnames(markers_rna))
colnames(markers_motifs) <- paste0("motif.", colnames(markers_motifs))
markers_rna$gene <- markers_rna$RNA.feature
markers_motifs$gene <- ConvertMotifID(immune.combined4[["ATAC"]], id = motif.names)

# Provide the path to the downloaded mapping file
mapping_file_path <- "/Users/rubenprins/Downloads/HOM_MouseHumanSequence.rpt.txt"

# Read the mapping file into a data frame
mapping_df <- read.delim(mapping_file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Define a function to convert mouse gene names to human gene names
convertMouseToHuman <- function(mouse_gene) {
  mouse_rows <- which(mapping_df$Symbol == mouse_gene)
  if (length(mouse_rows) > 0) {
    human_gene <- mapping_df$Symbol[mouse_rows + 1]
    if (toupper(mouse_gene) == mouse_gene) {
      return(mouse_gene)  # Return mouse_gene if found in all caps
    } else {
      return(human_gene)  # Return human_gene for other cases
    }
  } else {
    return(mouse_gene)  # Return a specific value when gene is not found
  }
}

convertMouseToHuman("gm1025")
# Convert gene names in markers_rna and markers_motifs data frames
markers_rna$human_gene <- sapply(markers_rna$gene, convertMouseToHuman)
markers_motifs$human_gene <- sapply(markers_motifs$gene, convertMouseToHuman)

# Remove rows with NA in the human_gene column of markers_rna
markers_rna <- markers_rna[!is.na(markers_rna$human_gene), ]


# Define the topTFs function
topTFs <- function(stim, padj.cutoff = 5e-2) {
  ctmarkers_rna <- dplyr::filter(
    markers_rna, RNA.group %in% stim, RNA.padj < padj.cutoff, RNA.logFC > 0, RNA.auc >0.6) %>% 
    arrange(desc(RNA.auc))
  
  ctmarkers_motif <- dplyr::filter(
    markers_motifs, motif.group %in% stim, motif.padj < padj.cutoff, motif.logFC > 0, motif.auc >0.6) %>% 
    arrange(desc(motif.auc))
  # Convert human_gene column in ctmarkers_motif to character
  ctmarkers_motif$human_gene <- as.character(ctmarkers_motif$human_gene)
  ctmarkers_rna$human_gene <- as.character(ctmarkers_rna$human_gene)
  top_tfs <- inner_join(
    x = ctmarkers_rna[, c(2, 12, 6, 7)], 
    y = ctmarkers_motif[, c(2, 1, 12, 6, 7)], by = "human_gene"
  )
  top_tfs$avg_auc <- (top_tfs$RNA.auc + top_tfs$motif.auc) / 2
  top_tfs <- arrange(top_tfs, desc(avg_auc))
  top_tfs <- na.omit(top_tfs)
  return(top_tfs)
}

head(topTFs("0"), 100)
mostwanted <- topTFs("0")
mostwanted

#scatterplot(1:25,mostwanted$avg_auc)


mostwanted$motif.neglogpval <- -log(mostwanted$motif.pval, base = 10)
mostwanted <- mostwanted[order(mostwanted$motif.neglogpval, decreasing = TRUE), ]
view(mostwanted)

pTFil23 <- plot(1:nrow(mostwanted),mostwanted$motif.neglogpval, type = "p", xlab = "Rank", ylab = "-10Log(P_Value)")+
  with(mostwanted, text(x=1:3, y=mostwanted$motif.neglogpval[1:3], labels = mostwanted$human_gene[1:3], pos = 4, cex = .5))

pTFil23


# Define the topTFs function
topTFsuppress <- function(stim, padj.cutoff = 5e-2) {
  ctmarkers_rna <- dplyr::filter(
    markers_rna, RNA.group %in% stim, RNA.padj < padj.cutoff, RNA.logFC < 0, RNA.auc <0.4) %>% 
    arrange(desc(RNA.auc))
  
  ctmarkers_motif <- dplyr::filter(
    markers_motifs, motif.group %in% stim, motif.padj < padj.cutoff, motif.logFC < 0, motif.auc <0.4) %>% 
    arrange(desc(motif.auc))
  # Convert human_gene column in ctmarkers_motif to character
  ctmarkers_motif$human_gene <- as.character(ctmarkers_motif$human_gene)
  ctmarkers_rna$human_gene <- as.character(ctmarkers_rna$human_gene)
  top_tfs <- inner_join(
    x = ctmarkers_rna[, c(2, 12, 6, 7)], 
    y = ctmarkers_motif[, c(2, 1, 12, 6, 7)], by = "human_gene"
  )
  top_tfs$avg_auc <- (top_tfs$RNA.auc + top_tfs$motif.auc) / 2
  top_tfs <- arrange(top_tfs, desc(avg_auc))
  top_tfs <- na.omit(top_tfs)
  return(top_tfs)
}

head(topTFsuppress("0"), 50)
mostwantedsuppress <- topTFsuppress("0")
mostwantedsuppress



mostwantedsuppress$motif.neglogpval <- -log(mostwantedsuppress$motif.pval, base = 10)
mostwantedsuppress <- mostwantedsuppress[order(mostwantedsuppress$motif.neglogpval, decreasing = TRUE), ]
view(mostwantedsuppress)

pTFsuppressil23 <- plot(1:nrow(mostwantedsuppress),mostwantedsuppress$motif.neglogpval, type = "p", xlab = "Rank", ylab = "-10Log(P_Value)")+
  with(mostwantedsuppress[1:3,], text(x=1:3, y=mostwantedsuppress$motif.neglogpval[1:3], labels = mostwantedsuppress$human_gene[1:3], pos = 1, cex = .5))

write.csv(mostwanted, "/Users/rubenprins/Library/Mobile Documents/com~apple~CloudDocs/USC stuff/R01 edits/Figures R01/TFassociatedwithIL23.csv", row.names=FALSE)

write.csv(mostwantedsuppress, "/Users/rubenprins/Library/Mobile Documents/com~apple~CloudDocs/USC stuff/R01 edits/Figures R01/TFassociatedwithIL23suppress.csv", row.names=FALSE)


plot_grid(pTFil23,pTFsuppressil23)


library(ggplot2)
library(ggrepel)
library(scales)
# Create a vector of x-axis values with increments of 5 up to the maximum value or 30
x_values <- c(seq(0, max(nrow(mostwanted), 30), by = 5), 30)

# Convert the x-axis values to character format
x_labels <- as.character(x_values)

# Modify the aes mapping to use the x_labels instead of 1:nrow(mostwanted)
pTFil23 <- ggplot() +
  geom_point(data = mostwanted, aes(x = 1:nrow(mostwanted), y = motif.neglogpval), color = "blue") +
  geom_text_repel(data = mostwanted, aes(x = 1:nrow(mostwanted), y = motif.neglogpval, label = ifelse(motif.neglogpval > 1, human_gene, "")), 
                  vjust = -1, size = 5, color = "blue") +
  geom_point(data = mostwantedsuppress, aes(x = 1:nrow(mostwantedsuppress), y = -mostwantedsuppress$motif.neglogpval), color = "red") +
  geom_text_repel(data = mostwantedsuppress, aes(x = 1:nrow(mostwantedsuppress), y = -mostwantedsuppress$motif.neglogpval, label = ifelse(mostwantedsuppress$motif.neglogpval > 1, mostwantedsuppress$human_gene, "")), 
                  vjust = 1, size = 5, color = "red") +
  labs(x = "Rank", y = "-10Log(P_Value)") +
  scale_color_manual(values = c("blue", "red")) +
  coord_cartesian(ylim = c(-40, 40)) +
  scale_x_continuous(breaks = seq(0, 30, by = 5), limits = c(0, 30)) +
  scale_y_continuous(labels = function(x) comma(abs(x), accuracy = 1, trim = TRUE), 
                     breaks = seq(-40, 40, by = 10),
                     trans = "reverse", expand = c(0, 0)) +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10)) +
  theme_classic2() +
  geom_hline(yintercept = 0, aes(xmin = 0, xmax = 27, color = "black", linetype = "solid")) +
  theme(axis.text = element_text(size = 16),  # Increase the size of axis text
        axis.title = element_text(size = 16))  # Increase the size of axis titles

pTFil23


#RNA visualize

convertHumanToMouse <- function(human_gene) {
  human_rows <- which(mapping_df$Symbol == human_gene)
  if (length(human_rows) > 0) {
    mouse_gene <- mapping_df$Symbol[human_rows - 1]
    return(mouse_gene)  # Return mouse_gene 
  } else {
    return(human_gene)  # Return a specific value when gene is not found
  }
}


mostwanted$interest_gene <- sapply(mostwanted$human_gene, convertHumanToMouse)
listTFil23 <- mostwanted$interest_gene
DefaultAssay(immune.combined2) <- "RNA"
p50 <- FeaturePlot(object = immune.combined2, 
                  features = c("Il23a","Stat3","Nr4a2","Fos"),
                  cols = c("grey", "green", "red","purple"), 
                  reduction = "umap",
                  order = TRUE,
                  pt.size = .5) 

DefaultAssay(immune.combined2) <- "chromvar"
p51 <- FeaturePlot(
  object = immune.combined2,
  cols = c("grey", "green", "red","purple"), 
  features = mostwanted$motif.feature[1:10],
  reduction = "umap",
  order = TRUE,
  pt.size = 0.5,
  min.cutoff = "q10",
  max.cutoff = "q90")

p52<-FeaturePlot(object = immune.combined2, 
                 cols = c("grey", "green", "red","purple"), 
                features = c("Klf15","Klf3","Zbtb14","Tcfl5", "Il23a"), 
                reduction = "umap",
                order = TRUE,
                pt.size = 0.5,
                min.cutoff = "q10",
                max.cutoff = "q90")

plot_grid(p50,p51)
patchwork::wrap_plots(p50, p51, ncol = 1)

f1 = FeaturePlot(immune.combined2,reduction = "umap", features = c("MA0102.4", "MA0101.1"),pt.size = .5, blend = TRUE,order = TRUE, cols = c("black","yellow", "blue"), blend.threshold = 1, max.cutoff = 1)
f1 & DarkTheme()
 














stim.data22 <-Read10X_h5(filename = '/Users/rubenprins/Downloads/HPV-1-GEX_out/filtered_feature_bc_matrix.h5',
                         use.names = TRUE,
                         unique.features = TRUE)
view(stim.data22)
stim.datarna1 <-  stim.data22

stim22 <- CreateSeuratObject(counts = stim.datarna1, project = "HPV", min.cells = 5)


#QC
stim22[["percent.mt"]] <- PercentageFeatureSet(stim22, pattern = "^mt-")
VlnPlot(stim22, features = c("nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0)
stim22 <- subset(stim22, subset = nFeature_RNA > 1000 & nFeature_RNA < 3500 & percent.mt < 30 & nCount_RNA >2000)


ctrl.data11 <- Read10X_h5(filename = '/Users/rubenprins/Library/Mobile Documents/com~apple~CloudDocs/USC stuff/IL-17 and IL-23 Project/10xGenomics Multiome sequencing/GFP-2-GEX_out/filtered_feature_bc_matrix.h5', use.names = TRUE, unique.features = TRUE)


ctrl.datarna1 <-  ctrl.data11

ctrl11 <- CreateSeuratObject(counts = ctrl.datarna1, project = "GFP", min.cells = 5)


#QC
ctrl11[["percent.mt"]] <- PercentageFeatureSet(ctrl11, pattern = "^mt-")
VlnPlot(ctrl11, features = c( "nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0)
ctrl11 <- subset(ctrl11, subset =  nFeature_RNA > 1000 & nFeature_RNA < 3500 & percent.mt < 30 & nCount_RNA >2000)


stim22$stim <- "HPV1"

stim22 <- NormalizeData(stim22, assay = 'RNA', verbose = FALSE)
stim22 <- FindVariableFeatures(stim22, assay = "RNA", selection.method = "vst", nfeatures = 2000)


ctrl11$stim <- "GFP1"
ctrl11 <- NormalizeData(ctrl11, assay = "RNA", verbose = FALSE)
ctrl11 <- FindVariableFeatures(ctrl11, assay = "RNA", selection.method = "vst", nfeatures = 2000)

DefaultAssay(immune.combined) <-"integrated"
apc.query <- ctrl11
apc.anchors <- FindTransferAnchors(reference = immune.combined, query = apc.query,
                                        dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = apc.anchors, refdata = immune.combined$seurat_clusters,
                            dims = 1:30)
apc.query <- AddMetaData(apc.query, metadata = predictions)

immune.combined <- RunUMAP(immune.combined, dims = 1:30, reduction = "pca", return.model = TRUE)
apc.query <- MapQuery(anchorset = apc.anchors, reference = immune.combined, query = apc.query,
                            reference.reduction = "pca", reduction.model = "umap")


pextra <- DimPlot(apc.query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE)

apc.query1 <- stim22
apc.anchors1 <- FindTransferAnchors(reference = immune.combined, query = apc.query1,
                                   dims = 1:30, reference.reduction = "pca")
predictions1 <- TransferData(anchorset = apc.anchors1, refdata = immune.combined$seurat_clusters,
                            dims = 1:30)
apc.query1 <- AddMetaData(apc.query1, metadata = predictions1)

immune.combined <- RunUMAP(immune.combined, dims = 1:30, reduction = "pca", return.model = TRUE)
apc.query1 <- MapQuery(anchorset = apc.anchors1, reference = immune.combined, query = apc.query1,
                      reference.reduction = "pca", reduction.model = "umap")


pextra1 <- DimPlot(apc.query1, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
                  label.size = 3, repel = TRUE)




immune.anchors <- FindIntegrationAnchors(object.list = list(ctrl11, stim22), dims = 1:10)
#to_integrate <- Reduce(intersect, lapply(immune.anchors@object.list, rownames))
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)






#Subset only Il23a cells
ctrl33 <- subset(ctrl12, subset = Il23a >0| Btla >0)#|Tmsb4x <5
stim33 <- subset(stim21, subset = Il23a >0| Btla >0)#|Tmsb4x <5
ctrl33 <- FindVariableFeatures(ctrl33, selection.method = "vst", nfeatures = 2000)
stim33 <- FindVariableFeatures(stim33, selection.method = "vst", nfeatures = 2000)

immune.anchors3 <- FindIntegrationAnchors(object.list = list(ctrl33, stim33), dims = 1:10, k.filter = 71)
#to_integrate <- Reduce(intersect, lapply(immune.anchors@object.list, rownames))
immune.combined3 <- IntegrateData(anchorset = immune.anchors3, dims = 1:20, k.weight = 70)
DefaultAssay(immune.combined3) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined3 <- ScaleData(immune.combined3, verbose = FALSE)
immune.combined3 <- RunPCA(immune.combined3, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
immune.combined3 <- RunUMAP(immune.combined3, reduction = "pca", dims = 1:20)
immune.combined3 <- FindNeighbors(immune.combined3, reduction = "pca", dims = 1:20)
immune.combined3 <- FindClusters(immune.combined3, resolution = 0.5)


#dotplot difference between GFP and HPV genes
DotPlot(immune.combined,idents = c("0","1","2","3","4"), features = c("Cxcl2", "Cd38","Cxcl1", "Nfkb1", "Cd80","Il1b","Mrc1","Cd44", "Mertk", "Csf1r"), split.by = "stim", cols = c("blue", "red","green"))


















# call peaks using MACS2
peaks1 <- Signac::CallPeaks(
  immune.combined,
  group.by = "stim",
  assay = "ATAC",
  macs2.path = "/usr/local/Caskroom/miniconda/base/envs/py27/bin/macs2"
)



# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks1 <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks1 <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)


DefaultAssay(immune.combined) <- "ATAC"

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(immune.combined),
  features = peaks1,
  cells = colnames(immune.combined)
)


# create a new assay using the MACS2 peak set and add it to the Seurat object
immune.combined[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = Fragments(immune.combined),
  annotation = annotations
)


DefaultAssay(immune.combined) <- "peaks"





#Alternative matif analysis
library(Signac)
library(Seurat)
library(JASPAR2022)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
#library(BSgenome.Mmusculus.UCSC.mm39)
library(patchwork)
set.seed(1234)
library(motifmatchr)
library(ggseqlogo)
library(chromVAR)
#BiocManager::install('motifmatchr')
#BiocManager::install('TFBSTools')
#BiocManager::install('BSgenome.Mmusculus.UCSC.mm10')
#BiocManager::install('chromVAR')


# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2022,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

DefaultAssay(immune.combined) <- "peaks"
# add motif information
immune.combined <- AddMotifs(
  object = immune.combined,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

immune.combined <- RunChromVAR(
object = immune.combined,
genome = BSgenome.Mmusculus.UCSC.mm10
)

DefaultAssay(immune.combined) <- 'chromvar'

p52 <- FeaturePlot(
  object = immune.combined,
  cols = c("grey", "green", "red","purple"), 
  features = c("MA1515.1"),
  reduction = "umap",
  order = TRUE,
  pt.size = 0.2,
  min.cutoff = "q10",
  max.cutoff = "q90")
p52


DefaultAssay(immune.combined.sct) <- "integrated"
immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30, verbose = FALSE)
immune.combined.sct <- FindNeighbors(immune.combined.sct, reduction = "pca", dims = 1:30)
immune.combined.sct <- FindClusters(immune.combined.sct, resolution = 0.3)

immune.combined.sct <- immune.combined
DefaultAssay(immune.combined.sct) <- "ATAC"
immune.combined.sct <- RunTFIDF(immune.combined.sct)
immune.combined.sct <- FindTopFeatures(immune.combined.sct, min.cutoff = 'q0')
immune.combined.sct <- RunSVD(immune.combined.sct)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

immune.combined.sct <- FindMultiModalNeighbors(immune.combined.sct, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:50))
immune.combined.sct <- RunUMAP(immune.combined.sct, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
immune.combined.sct <- FindClusters(immune.combined.sct, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

DefaultAssay(immune.combined) <- "RNA"
all_markers <-FindAllMarkers(immune.combined,
                             min.pct =  0.65)

View(all_markers)
DefaultAssay(immune.combined.sct) <- "RNA"
all_markersrna <-FindAllMarkers(immune.combined,
                             min.pct =  0.65)

View(all_markersrna)

FeaturePlot(immune.combined.sct, features = c("MA0004.1","Il23a"), reduction = "wnn.umap", cols = c("gray","green","red","purple"), pt.size = .3, order = T, blend = T, max.cutoff = .5)

plot3 <- DimPlot(immune.combined.sct,
                     reduction = "wnn.umap",
                     label = TRUE,
                     pt.size = 1) 

plot4 <- FeaturePlot(immune.combined.sct, 
            features = "Mef2c",
            cols = c("black","yellow"), 
            reduction = "wnn.umap",
            order = TRUE,
            pt.size = 1) & DarkTheme()
