library(shiny)
library(bslib)

# Define UI for app that draws a histogram ----
ui <- page_sidebar(
  # App title ----
  title = "Multiome Analysis 2024",
  # Sidebar panel for inputs ----
 fileInput("upload", "scRNAfile", placeholder = "scRNAfile"),
 tableOutput("files")
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  output$files <- renderTable(input$upload)
}

shinyApp(ui = ui, server = server)
