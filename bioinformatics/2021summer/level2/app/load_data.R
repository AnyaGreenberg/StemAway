library(shiny)
library(affy)

options(shiny.maxRequestSize=600*1024^2)

ui <- navbarPage(
  theme = bslib::bs_theme(bootswatch = 'cosmo'),
  'Bioinformatics \n Workflow',
  
  
  tabPanel('Load Data',
           sidebarLayout(
             sidebarPanel(
               fileInput('affy_file', 'Upload an RDS of an affyBatch object', multiple=FALSE, accept=c('.rds')),
               fileInput('meta_file', 'Upload a metadata file', multiple=FALSE, accept=c('.csv', '.tsv')),
               h3('Fun Facts'),
               p("RDS files are R data files which don't conform to normal file types. In this case, affybatch objects represent Affymetrix GeneChip probe level data. And the object has a special structure that allows this data to be accessible in R while keeping memory and space costs minimal."),
               p("CSVs and TSV are standard data files and can be viewed (if small enough) in a normal text editor. The difference between the 2 file types is how the data is separated. CSV stands for comma separated values and TSV stands for tab separated values.")
             ),
             
             mainPanel(
               h1('affyBatch data'),
               verbatimTextOutput('affybatch'),
               h1('Metadata'),
               tableOutput('metadata')
             )
           )
   ),
  
  
  tabPanel('Quality Control',
           'Here there will be QC output'
  ),
  
  
  tabPanel('DGE Analysis',
           'Here there will be GSE output'
  ),
  
  
  tabPanel('FunctionalAnalysis',
           'Here there will be FA output and it will probably contain subtabs'
  )
  
)

server <- function(input, output, session){
  # Verify .rds file
  affy_data <- reactive({
    req(input$affy_file)

    ext <- tools::file_ext(input$affy_file$name)

    switch(ext,
           rds = readRDS(input$affy_file$datapath),
           validate('Invalid file. Please upload an .rds file'))
  })

  # Load file contents and print object
  output$affybatch <- renderPrint(print(affy_data()))
  
  
  # Verify .csv or .tsv file
  meta_data <- reactive({
    req(input$meta_file)
    
    ext <- tools::file_ext(input$meta_file$name)
    
    switch(ext,
           csv = vroom::vroom(input$meta_file$datapath, delim = ','),
           tsv = vroom::vroom(input$meta_file$datapath, delim = '\t'),
           validate('Invalid file. Please upload .csv, or .tsv file'))
  })
  
  # Load file contents and print head
  output$metadata <- renderTable(meta_data())
}

shinyApp(ui, server)


### FUNCTIONS
# normalize <- (data){
#   return mas5()
# }