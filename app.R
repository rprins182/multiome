library(shiny)
library(bslib)
library(Seurat)
library(tidyverse)

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

# Increase the maximum upload size to 300 MB
options(shiny.maxRequestSize = 300*1024^2)

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

# Define UI for app that draws a histogram
ui <- fluidPage(
  # App title
  titlePanel("Multiome Analysis 2024"),
  
  textInput("sample_project_label", "Enter sample project label:", ""),
  textInput("control_project_label", "Enter control project label:", ""),
  
  # Sidebar panel for inputs
  fileInput("sample_h5", "Upload sample scRNA.h5 file", accept = ".h5"),
  fileInput("sample_fragment", "Upload sample fragment.tsv file", accept = ".h5"),
  fileInput("control_h5", "Upload sample scRNA.h5 file", accept = ".h5"),
  fileInput("control_fragment", "Upload sample fragment.tsv file", accept = ".h5"),
  textOutput("upload_status"),
  
  #Output
  titlePanel("Seurat Violin plot QC")
  plotlyOutput("violin_plots"),
  
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  # Reactive expression to read the input file
  stim_data <- reactive({
    req(input$sample_h5)
    req(input$control_h5)
    req(input$sample_fragment)
    req(input$control_fragment)
    output$upload_status <- renderText("Uploading and processing the file, please wait...")
    

    sample_s4_matrices <- read_10x(input$sample_h5)
    control_s4_matrices <- read_10(input$control_h5)
    
    output$upload_status <- renderText("File successfully uploaded and start the analyze.")

    sample_s4_matrix_rna <- columns_selection(sample_s4_matrices, `Gene Expression`)
    sample_s4_matrix_atac <- columns_selection(sample_s4_matrices, `Peaks`)

    control_s4_matrix_rna <- columns_selection(control_s4_matrices, `Gene Expression`)
    control_s4_matrix_atac <- columns_selection(control_s4_matrices, `Peaks`)

    sample_seurat_obj <- CreateSeuratObject(counts = sample_s4_matrix_rna, project = sample_project_label, min.cells = 5)
    control_seurat_obj <- CreateSeuratObject(counts = control_s4_matrix_rna, project = control_project_label, min.cells = 5)

    # Setting the table
    sample_s4_matrix_atac_named_coord <- atac_table_settings(sample_s4_matrix_atac)
    control_s4_matrix_atac_named_coord <- atac_table_settings(control_s4_matrix_atac)

    # Annotation
    genome = "mm39"
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotations) <- "UCSC"
    genome(annotations) <- genome
    
    sample_atac_fragments <- sample_fragment
    control_atac_fragments <- control_fragment

    sample_chrom_assay <- (sample_s4_matrix_atac_named_coord, annotation, genome, sample_atac_fragments)
    control_chrom_assay <- (control_s4_matrix_atac_named_coord, annotation, genome, control_atac_fragments)

    # FILTERING
    # Add new column with chrom_assay object
    sample_seurat_obj <- add_to_seurat(sample_seurat_obj, "ATAC", sample_chrom_assay)
    control_seurat_obj <- add_to_seurat(control_seurat_obj, "ATAC", control_chrom_assay)

    search_pattern = "^mt-"
    sample_seurat_obj <- add_to_seurat(sample_seurat_obj, "percent.mt", pattern=search_pattern)
    control_seurat_obj <- add_to_seurat(control_seurat_obj, "percent.mt", pattern=search_pattern)


    #Violin plot visualization for the percentage
    sample_qc_vlnplot <- VlnPlot(sample_seaurat_obj, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
            log = TRUE, pt.size = 0)
    control_qc_vlnplot <- VlnPlot(control_seaurat_obj, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
            log = TRUE, pt.size = 0)

    return c(sample_qc_vlnplot, control_qc_vlnplot)
  })

  output$violin_plots <- renderPlotly({
    stim_data(input$sample_h5, input$control_h5, input$sample_fragment, input$control_fragment)
  })
}

shinyApp(ui = ui, server = server)

