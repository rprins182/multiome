library(Seurat)
library(tidyverse)
library(cowplot)
library(ggplot2)
library(patchwork)
library(dplyr)
library(plotly)
#library(fcoex)
library(ggpubr)
#library(SingleR)
#library(glmGamPoi)
#library(sctransform)
#library("pheatmap")

#multiome
#library(Signac)
#library(EnsDb.Mmusculus.v79)
#library("LOLA")
#library(Biostrings)
#library(JASPAR2020)
#library(MotifDb)
#library(BSgenome.Mmusculus.UCSC.mm10)

#library(rtracklayer)
#library(BSgenome)
#library(GenomeInfoDb)
#library(PWMEnrich)
#library(PWMEnrich.Mmusculus.background)
#library('glmGamPoi')
#library("DESeq2")

#library("ggrepel")

# Increase the maximum upload size to 10 GB
options(shiny.maxRequestSize = 10 * 1024^3)  # 10 GB in bytes

# Define function to read 10x data from .h5 file
read_10x <- function(file_path) {
  data_file <- Read10X_h5(filename = file_path,
                          use.names = TRUE,
                          unique.features = TRUE)
  return(data_file)
}

# Function to select specific columns from S4 matrices
columns_selection <- function(s4_matrices, column_name) {
  selected_column <- s4_matrices[[column_name]]
  return(selected_column)
}

# Function to set ATAC table settings
atac_table_settings <- function(s4_matrix) {
  table_chr_coord <- StringToGRanges(rownames(s4_matrix), sep = c(":", "-"))
  table_chr_named_coord <- seqnames(table_chr_coord) %in% standardChromosomes(table_chr_coord)
  s4_matrix_atac_named_coord <- s4_matrix[as.vector(table_chr_named_coord), ]
  return(s4_matrix_atac_named_coord)
}

# Function to create chromatin assay
chrom_assay <- function(s4_matrix_atac_named_coord, annotations, genome, atac_fragments) {
  chrom_assay <- CreateChromatinAssay(
    counts = s4_matrix_atac_named_coord,
    sep = c(":", "-"),
    genome = genome,
    fragments = atac_fragments,
    min.cells = 10,
    annotation = annotations
  )
  return(chrom_assay)
}
# Function to add data to Seurat object
add_to_seurat <- function(seurat_object, column_to_add, data_to_add = "", pattern = FALSE) {
  if (pattern) {
    data_to_add <- PercentageFeatureSet(seurat_object, pattern = pattern)  # Assuming this function call is correct
  }
  seurat_object[[column_to_add]] <- data_to_add
  return(seurat_object)
}

# Define UI for app
ui <- fluidPage(
  titlePanel("Multiome Analysis 2024"),
  
  # Input fields
  textInput("sample_project_label", "Enter sample project label:", ""),
  textInput("control_project_label", "Enter control project label:", ""),
  
  fileInput("sample_h5", "Upload sample scRNA.h5 file", accept = ".h5"),
  fileInput("sample_fragment", "Upload sample fragment.tsv file", accept = ".tsv"),
  
  fileInput("control_h5", "Upload control scRNA.h5 file", accept = ".h5"),
  fileInput("control_fragment", "Upload control fragment.tsv file", accept = ".tsv"),
  
  # Output elements
  textOutput("upload_status"),
  tableOutput("rna"),
  
  # Violin plot section
  titlePanel("Seurat Violin plot QC"),
  plotlyOutput("violin_plots")
)

# Define server logic
server <- function(input, output) {
  
  # Reactive expression to process data files
  stim_data <- reactive({
    req(input$sample_h5)
    req(input$control_h5)
    req(input$sample_fragment)
    req(input$control_fragment)
    
    output$upload_status <- renderText("Uploading and processing the files, please wait...")
    
    # Read data files
    sample_s4_matrices <- read_10x(input$sample_h5$datapath)
    control_s4_matrices <- read_10x(input$control_h5$datapath)
    
    # Select columns
    sample_s4_matrix_rna <- columns_selection(sample_s4_matrices, "Gene Expression")
    sample_s4_matrix_atac <- columns_selection(sample_s4_matrices, "Peaks")
    
    control_s4_matrix_rna <- columns_selection(control_s4_matrices, "Gene Expression")
    control_s4_matrix_atac <- columns_selection(control_s4_matrices, "Peaks")
    
    # Create Seurat objects
    sample_seurat_obj <- CreateSeuratObject(counts = sample_s4_matrix_rna, project = isolate(input$sample_project_label), min.cells = 5)
    control_seurat_obj <- CreateSeuratObject(counts = control_s4_matrix_rna, project = isolate(input$control_project_label), min.cells = 5)
    
    # Settings for ATAC table
    sample_s4_matrix_atac_named_coord <- atac_table_settings(sample_s4_matrix_atac)
    control_s4_matrix_atac_named_coord <- atac_table_settings(control_s4_matrix_atac)
    
    # Annotation
    genome <- "mm39"
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotations) <- "UCSC"
    genome(annotations) <- genome
    
    sample_atac_fragments <- input$sample_fragment$datapath
    control_atac_fragments <- input$control_fragment$datapath
    
    # Chromatin assay
    sample_chrom_assay <- chrom_assay(sample_s4_matrix_atac_named_coord, annotations, genome, sample_atac_fragments)
    control_chrom_assay <- chrom_assay(control_s4_matrix_atac_named_coord, annotations, genome, control_atac_fragments)
    
    # Add data to Seurat objects
    sample_seurat_obj <- add_to_seurat(sample_seurat_obj, "ATAC", sample_chrom_assay)
    control_seurat_obj <- add_to_seurat(control_seurat_obj, "ATAC", control_chrom_assay)
    
    # Add percent.mt column
    search_pattern <- "^mt-"
    sample_seurat_obj <- add_to_seurat(sample_seurat_obj, "percent.mt", pattern = search_pattern)
    control_seurat_obj <- add_to_seurat(control_seurat_obj, "percent.mt", pattern = search_pattern)
    
    # Violin plot data
    sample_qc_vlnplot <- VlnPlot(sample_seurat_obj, features = c("nCount_ATAC", "nCount_RNA", "percent.mt"), ncol = 3,
                                 log = TRUE, pt.size = 0)
    control_qc_vlnplot <- VlnPlot(control_seurat_obj, features = c("nCount_ATAC", "nCount_RNA", "percent.mt"), ncol = 3,
                                  log = TRUE, pt.size = 0)
    
    # Return plots
    list(sample_qc_vlnplot, control_qc_vlnplot)
  })
  
  # Render violin plots
  output$violin_plots <- renderPlotly({
    stim_data()
  })
  
  # Observer to update upload status
  observe({
    req(input$sample_h5)
    req(input$control_h5)
    req(input$sample_fragment)
    req(input$control_fragment)
    output$upload_status <- renderText("Files successfully uploaded and analysis started.")
  })
}

# Run the application
shinyApp(ui = ui, server = server)

