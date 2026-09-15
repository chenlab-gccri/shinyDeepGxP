#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
library(htmltools)
library(htmlwidgets)
library(rintrojs)
library(shiny)
library(shinyBS)
library(shinydashboardPlus)
library(shinycssloaders)
library(shinydashboard)
library(shinyWidgets)
library(plotly)
library(DT)
library(VennDiagram)
library(tidyr)
library(stringr)
library(dplyr)
library(shinyjs)
library(igraph)
library(zip)
library(stats)

ui <- htmlTemplate(
  "index.html",
  privacytag1 = actionLink(
    inputId = "privacytag1",
    label = "Privacy Policy",
    style = "color:white;border-bottom: 1px solid #a9529e;"
  ),
  
  help_explore_full = includeHTML("www/helptab2_ExploreModel.html"),
  help_predict_full = includeHTML("www/helptab3_PredictProteins.html"),
  
  #Input Reset Buttons
  # Reset.p = actionBttn(
  #   inputId = "Reset.p",
  #   label = "Reset File Upload",
  #   size = "sm",
  #   style = "material-flat",
  # ),
  
  Reset.p = actionButton(
    inputId = "Reset.p",
    label = "Reset Upload",
    icon = icon("trash"),
    class = "upload-reset-btn"
  ),
  
  #Show description for each plot
  # Table1-1: nn.table.q
  # Figure1-3: nn.network.q
  nnboth1 = circleButton(
    inputId = "nnboth1",
    icon = icon("circle-question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  nnboth3 = circleButton(
    inputId = "nnboth3",
    icon = icon("circle-question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  
  # help_input1 = actionLink(
  #   inputId = "help_input1",
  #   label = "Download example data.",
  #   style = "color:blueviolet; border-bottom: 1px solid #a9529e;"
  # ),

  # help_input1 = tagList(
  #   actionLink(
  #     inputId = "help_input1",
  #     label = "Download example data.",
  #     style = "color:blueviolet; border-bottom: 1px solid #a9529e;"
  #   ),
  #   downloadButton(
  #     outputId = "download_example",
  #     label = NULL,  # no label since we trigger it manually
  #     style = "visibility: hidden; height: 1px;"  # hide the actual button
  #   )
  # ),
  
  # download_prediction_metadata = downloadButton(
  #   outputId = "download_prediction_metadata",
  #   label = "Download predicted proteins with clusters and cell types"
  # ),
  
  download_prediction_metadata = tagList(
    conditionalPanel(
      condition = "output.prediction_ready == false",
      actionButton(
        inputId = "download_prediction_metadata_check",
        label = "Download predicted proteins with clusters and cell types",
        icon = icon("download")
      )
    ),
    conditionalPanel(
      condition = "output.prediction_ready == true",
      downloadButton(
        outputId = "download_prediction_metadata",
        label = "Download predicted proteins with clusters and cell types"
      )
    )
  ),
  
  download_all_umaps = downloadButton(
    outputId = "download_all_umaps",
    label = "Download all UMAP plots"
  ),
  
  
  # download_prediction_metadata = tagList(
  #   actionButton(
  #     inputId = "download_prediction_metadata_check",
  #     label = "Download predicted proteins with clusters and cell types",
  #     icon = icon("download"),
  #     class = "upload-example-btn"
  #   ),
  #   downloadButton(
  #     outputId = "download_prediction_metadata",
  #     label = NULL,
  #     style = "display: none;"
  #   )
  # ),
  
  # question marks around each panel
  help1_1 = circleButton(
    inputId = "help1_1",
    icon = icon("circle-question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  
  help1_2 = circleButton(
    inputId = "help1_2",
    icon = icon("circle-question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  
  help1_3 = circleButton(
    inputId = "help1_3",
    icon = icon("circle-question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  
  help1_4 = circleButton(
    inputId = "help1_4",
    icon = icon("circle-question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  
  help2_1 = circleButton(
    inputId = "help2_1",
    icon = icon("circle-question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  
  help2_2 = circleButton(
    inputId = "help2_2",
    icon = icon("circle-question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  
  help2_3 = circleButton(
    inputId = "help2_3",
    icon = icon("circle-question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  
  help2_4 = circleButton(
    inputId = "help2_4",
    icon = icon("circle-question"),
    style = "color:darkblue; margin-bottom: 0.9vh;"
  ),
  
  helptab2p1 = actionLink(inputId = "helptab2p1", label = "Explore Model"),
  helptab3p1 = actionLink(inputId = "helptab3p1", label = "Predict Proteins"),

  help_input1 = actionLink(
    inputId = "helptab3p1",
    label = "Learn more about input data format and download examples",
    style = "color:blueviolet;
                                         border-bottom: 1px solid #a9529e;"
  ),
  # 
  ################## Module 1:Explore Model ##################
  run_interpret = actionButton(
    inputId = "run_interpret",
    label = "GO",
    style = "
    background-color: teal;
    border-color: teal;
    color: white;
    border-radius: 50%;
    font-weight: bold;
    width: 60px;
    height: 60px;
    font-size: 18px;
    padding: 0;
  "
  ),
  
  

  ################## Module 2:Predict Protein Abundance ##################
  #####Upload both files#####
  submitmodel = actionBttn(
    inputId = "submitmodel",
    label = "Submit",
    size = "md",
    icon = icon("play-circle"),
    style = "material-flat",
    color = "royal"
  ),
  
  submitexample = actionBttn(
    inputId = "submitexample",
    label = "Run Example",
    size = "md",
    icon = icon("rocket"),
    style = "material-flat",
    color = "primary"
  ),
  
  ####################################################################################################
  #Tab2:Prediction:
  #####Upload both files#####
  # submitmodel.p = actionBttn(
  #   inputId = "submitmodel.p",
  #   label = "Submit",
  #   size = "md",
  #   icon = icon("play-circle"),
  #   style = "material-flat",
  #   color = "royal"
  # ),
  
  # submitexample.p = actionBttn(
  #   inputId = "submitexample.p",
  #   label = "Run Example",
  #   size = "md",
  #   icon = icon("rocket"),
  #   style = "material-flat",
  #   color = "primary"
  # ),
  #example ID--> example.pd.both
  
  submitmodel.p = actionButton(
    inputId = "submitmodel.p",
    label = "GO",
    class = "go-circle-btn"
  ),
  
  submitexample.p = actionButton(
    inputId = "submitexample.p",
    label = "Run Example",
    icon = icon("flask"),
    class = "upload-example-btn"
  ),
  

  
  exp_data.p = fileInput(
    "exp_data.p",
    NULL,
    multiple = FALSE,
    accept = c("h5seurat", ".rds", ".txt"),
    width = "100%"
  ),

  #Show results:
  # Add output for the new prediction table
  # runtime_text = textOutput("runtime_text"),
  # gene_match_text = htmlOutput("gene_match_text"),
  prediction_summary = htmlOutput("prediction_summary"),
  
  prediction_table = shinycssloaders::withSpinner(dataTableOutput("prediction_table")),
  umap_clusters = shinycssloaders::withSpinner(plotOutput("umap_clusters")),
  umap_celltype1 = shinycssloaders::withSpinner(plotOutput("umap_celltype1")),
  umap_celltype2 = shinycssloaders::withSpinner(plotOutput("umap_celltype2")),

  ###### Module 1: Interpret Model #######
  
  cell_type = selectizeInput(inputId = "cell_type",label = NULL, choices = NULL, selected = NULL,
                                options = list(placeholder = "e.g. B", maxItems = 1),
                                multiple = FALSE ,width = "100%"),
  
  protein_name = selectizeInput(inputId = "protein_name",label = NULL,  choices = NULL, selected = NULL,
                                options = list(placeholder = "e.g. CD20", maxItems = 1),
                                multiple = FALSE, width = "100%"),

  geneset_name = selectizeInput(inputId = "geneset_name",label = NULL,  choices = NULL, selected = NULL,
                                options = list(placeholder = "e.g. GOBP", maxItems = 1),
                                multiple = FALSE, width = "100%"),
  
  #Show Panel 3-1 title
  text.gradient = textOutput("text.gradient"),
  
  filtered_gradient_table = shinycssloaders::withSpinner(dataTableOutput("filtered_gradient_table")),
  filtered_enrichment_table_pos = shinycssloaders::withSpinner(dataTableOutput("filtered_enrichment_table_pos")),
  filtered_enrichment_table_neg = shinycssloaders::withSpinner(dataTableOutput("filtered_enrichment_table_neg")),
  protein_violin_plot = shinycssloaders::withSpinner(plotlyOutput("protein_violin_plot")),
  gradient_density_plot = shinycssloaders::withSpinner(plotlyOutput("gradient_density_plot")),
  filtered_enrichment_plot_pos = shinycssloaders::withSpinner(plotlyOutput("filtered_enrichment_plot_pos")),
  filtered_enrichment_plot_neg = shinycssloaders::withSpinner(plotlyOutput("filtered_enrichment_plot_neg")),

  
  
  
)



