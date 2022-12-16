# Define UI for random distribution app ----

library(shiny)
library(shinyFiles)
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
      selectInput(inputId = "norm",
                  label = "Enter the normalization method you would like to use on the dataset", "all",
                  choices = c("all",
                              "log_norm",
                              "SCTransform",
                              "counts_per_million")),
      textInput(inputId = "features",
                label = "Enter the number of features in the dataset
                  This should be a positive integer:", "2000"),
    sliderInput(inputId = "resolution",
                label = "Enter the resolution that sets the 'granularity; of the clustering. Increased values correspond to a higher number of clusters
                  This should be a positive float between 0.4-1.2 :",
                value = 0.5, min = 0.4, max = 1.2),
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
                  tabPanel("Dataset Summary",
                           verbatimTextOutput("textOut")),
                  tabPanel("Plots",
                           uiOutput("plots")),
                  tabPanel("Normalized Density Plots",
                           h4("Visually inspect the results of the normalization method(s) selected."),
                           p("Warning: This step takes time to load"),
                           uiOutput('normplots')),
                  tabPanel("Clustering Results",
                           h4("Visually inspect the clustering outputs according to the normalization method selected."),
                           p("If normalization selection was 'all', default normalization is log normalized"),
                           p("Warning: This step takes time to load"),
                           plotOutput('clustplots')),
      )
    )
  )
)

# Define server logic for random distribution app ----
server <- function(input, output) {
  volumes <- (c(Home = fs::path_home(), getVolumes()()))
  shinyDirChoose(
    input,
    'dir',
    roots = volumes,
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
    req(input$norm)
    path_10x <- isolate(as.character(parseDirPath(volumes, input$dir)))
    plot_list <- autoSClust::runPreProcess(counts_matrix_path=path_10x, filterMetrics=input$qcMetric)
    plot_list
  })

  # render UI
  output$plots <- renderUI({
    req(input$dir)
    req(input$qcMetric)
    #print(length(runPreOut()))
    lapply(1:(length(runPreOut())-1), function(i) {
      # creates a unique ID for each plotOutput
      id <- paste0("plot_", i)
      plotOutput(outputId = id)
      # render each plot
      output$id <- renderPlot({
        x <- runPreOut()
        if (class(x[i]) != "Seurat") {
        x[[i]]
           }
        })
      })
    })

  # reactive expression to process data
  runNormOut <- reactive({
    req(input$dir)
    req(input$qcMetric)
    req(input$norm)
    path_10x <- isolate(as.character(parseDirPath(volumes, input$dir)))
    pre_list <- autoSClust::runPreProcess(counts_matrix_path=path_10x, filterMetrics="percent_mitochondrial")
    srt.lst <-  pre_list[sapply(pre_list, class) == "Seurat"]
    srt.obj <- srt.lst[[1]]
    normResult <- autoSClust::runNormalizarion(srt.data = srt.obj, norm=input$norm)
    normResult
  })

  # render UI
  output$normplots <- renderUI({
    req(input$dir)
    req(input$qcMetric)
    req(input$norm)
    if (input$norm != "all") {
      id <- "plot_1"
      plotOutput(outputId = id)
      output$id <- renderPlot({
        x <- runNormOut()
        x
      })
    }
    else {
      lapply(1:3, function(i) {
        # creates a unique ID for each plotOutput
        id <- paste0("plot_", i)
        plotOutput(outputId = id)
        # render each plot
        output$id <- renderPlot({
          x <- runNormOut()
          x[[i]]
        })
      })
    }
  })
  # reactive expression to process data
  runClustOut <- reactive({
    req(input$features)
    req(input$resolution)
    #path_10x <- isolate(as.character(parseDirPath(volumes, input$dir)))
    #pre_list <- autoSClust::runPreProcess(counts_matrix_path=path_10x, filterMetrics="percent_mitochondrial")
    #srt.lst <-  pre_list[sapply(pre_list, class) == "Seurat"]
    #srt.obj <- srt.lst[[1]]
    #normResult <- autoSClust::runNormalizarion(srt.data = srt.obj, norm=input$norm)
    normResult <- runNormOut()
    lst <-  normResult[sapply(normResult, class) == "Seurat"]
    obj <- lst[[1]]
    clustResult <- autoSClust::runClust(srt.data = obj, features=input$features, res=input$res)
    clustResult
  })

  # render UI
  output$clustplots <- renderPlot({
    req(input$dir)
    req(input$qcMetric)
    req(input$norm)
    req(input$features)
    req(input$resolution)
    clust <- runClustOut()
    DimPlot(clust, reduction = "umap")

  })
}
shinyApp(ui, server)
