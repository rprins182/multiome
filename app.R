library(shiny)
library(bslib)
library(Seurat)
library(tidyverse)

# Increase the maximum upload size to 300 MB
options(shiny.maxRequestSize = 300*1024^2)

# Define UI for app that draws a histogram
ui <- fluidPage(
  # App title
  titlePanel("Multiome Analysis 2024"),
  # Sidebar panel for inputs
  fileInput("upload", "Upload scRNA file", accept = ".h5"),
  textOutput("upload_status"),
  tableOutput("rna")
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  # Reactive expression to read the input file
  stim_data <- reactive({
    req(input$upload)
    output$upload_status <- renderText("Uploading and processing the file, please wait...")
    
    # Read data from the uploaded file
    data <- Read10X_h5(filename = input$upload$datapath,
                       use.names = TRUE,
                       unique.features = TRUE)
    
    # Extract RNA expression data
    rna_exprs <- data$`Gene Expression`
    
    # Convert to data frame for rendering
    exprs <- as.data.frame(as.matrix(rna_exprs))  # Convert to data frame
    
    # Update upload status
    output$upload_status <- renderText("File successfully uploaded and processed.")
    
    return(exprs)
  })
  
  output$rna <- renderTable({
    stim_data()
  })
}

shinyApp(ui = ui, server = server)
