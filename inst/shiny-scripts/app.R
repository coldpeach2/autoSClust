# Define UI for random distribution app ----


ui <- fluidPage(

  # App title ----
  titlePanel(tags$h1(tags$b("autoSClust:"),"Automated Clustering from user defined inputes to assess metrics")),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      tags$p("Welcome to Shiny App of autSClust R package."),
      # br() element to introduce extra vertical spacing ----
      br(),
      tags$b("Description: "),

      # br() element to introduce extra vertical spacing ----
      br(),
      br(),

      # input
      tags$p("Instructions: Below, enter or select values required to
              perform the analysis. Default values are shown. Then
              press 'Run'. Navigate through the different tabs to the
              right to explore the results."),

      # br() element to introduce extra vertical spacing ----
      br(),

      # input
      uiOutput("tab"),
      shinyDirButton("dir", "Input directory", "Upload"),
      verbatimTextOutput("dir", placeholder = TRUE),

      tags$p("Enter or select values required for clustering.
             Default values are shown."),
      selectInput(inputId = "qcMetric",
                label = "Enter the quality control metrics required for filtering the dataset", "all",
                choices = c("all",
                            "percent_mitochondrial",
                            "percent_ribosomal",
                            "percent_dissociation")),
      textInput(inputId = "nFeaturemax",
                label = "Enter the maximum number of genes per cell
                This should be a positive integer:", "1"),
      textInput(inputId = 'nCountsmax',
                  label = "Enter the maximum number of molecule counts per cell
                  This should be a positive integer:", "1"),
      textInput(inputId = "features",
                label = "Enter the number of features in the dataset
                  This should be a positive integer:", "1"),
    sliderInput(inputId = "resolution",
                label = "Enter the resolution that sets the 'granularity; of the clustering. Increased values correspond to a higher number of clusters
                  This should be a positive float between 0.4-1.2 :",
                value = 0.5, min = 0.4, max = 1.2),
    textInput(inputId = "dimensions",
              label = "Enter the dimensionality of the dataset where the majority of true signals are captured
                  This should be an int from range 1:int", "1:10"),
    selectInput(inputId = "metric",
                label = "Enter the metric you would like to asses the output of the clusters", "all",
                choices = c("all",
                            "dunn",
                            "rand")),

      # actionButton
      actionButton(inputId = "button2",
                   label = "Run"),

      # br() element to introduce extra vertical spacing -
      br(),

    ), # End of side pannel


    # Main panel for displaying outputs
    mainPanel(

      # Output: Tabet
      tabsetPanel(type = "tabs",
                  tabPanel("Input Summary",
                           verbatimTextOutput("textOut")),
                  tabPanel("Violin Plot",
                           h4("Instructions: Enter values and click 'Run' at the bottom left side."),
                           h4("Violin Plot of Initial Dataset:"),
                           br(),
                           plotOutput("vlnplot")),
                  tabPanel("Plots",
                           uiOutput("plots")),

                  tabPanel("Normalization Density Plots",
                           plotOutput('normplot')),
                  tabPanel("Model Selection",
                           fluidRow(
                             splitLayout(cellWidths = c("50%", "50%"), plotOutput('BICvalues'), plotOutput('ICLvalues')),
                             splitLayout(cellWidths = c("50%", "50%"), plotOutput('AIC3values'), plotOutput('AICvalues')),
                           )),

      )
    )
  )
)

# Define server logic for random distribution app ----
server <- function(input, output) {
  shinyDirChoose(
    input,
    'dir',
    roots = c(home = '~'),
    filetypes = c('', 'txt', 'mtx', "tsv", "csv")
  )
  global <- reactiveValues(datapath = getwd())

  dir <- reactive(input$dir)

  output$dir <- renderText({
    global$datapath
  })
  observeEvent(ignoreNULL = TRUE,
               eventExpr = {
                 input$dir
               },
               handlerExpr = {
                 if (!"path" %in% names(dir())) return()
                 home <- normalizePath("~")
                 global$datapath <-
                   file.path(home, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
               })

  output$seuratSummary <- renderTable({
    req(input$dir)
    data_matrix <- Read10X(data.dir = input$dir)
  },)

  # reactive expression to process data
  runPreOut <- reactive({
    req(input$dir)
    req(input$qcMetric)
    plot_list <- runPreProcess(counts_matrix_path=Read10X(data.dir = input$dir$datapath), filterMetrics=input$qcMetric)
    plot_list
    print(plot_list)
  })

  # render UI
  output$plots <- renderUI({
    req(input$dir)
    req(input$qcMetric)
    #print(length(runPreOut()))
    lapply(1:2, function(i) {
      # creates a unique ID for each plotOutput
      id <- paste0("plot_", i)
      plotOutput(outputId = id)
      # render each plot
      output$id <- renderPlot({
        x <- runPreOut()[i]
        patchwork(x + NoLegend() +
                    plot_annotation(theme=theme(plot.title = element_text(hjust = 0.5, face="bold"))))
        })
    })
    })
}
shinyApp(ui, server)
