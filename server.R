#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
app_start_time <- Sys.time()

library(shiny) 
library(shinyjs) 
library(shinyWidgets)
library(Seurat)
library(SeuratDisk)
library(Rmagic)
library(reticulate) 
library(future)
library(keras)
library(ggplot2)
library(DOSE)
library(rintrojs)
library(shinyBS)
library(shinydashboardPlus)
library(shinycssloaders)
library(shinydashboard)
library(plotly)
library(DT)
library(VennDiagram)
library(tidyr)
library(stringr)
library(dplyr)
library(igraph)
library(data.table)
library(tidyverse)
library(tibble)
library(methods)
library(zip)
library(stats)
library(arrow)
# library(Matrix)

# options(shiny.sanitize.errors = FALSE)
options(shiny.maxRequestSize = 1024*1024^2)
options(future.globals.maxSize = 10 * 1024^3)  # 10 GB
use_condaenv("r-DL_env", required = TRUE) 

## Source functions
source("function/utils_processing.R")
source("function/utils_io.R")

# Load model data
#trained_model <- load_model_hdf5("Predictions/CNN_X_magic_Y_protein_norm_renorm_NOgeneset_model_3.h5")
trained_model <- NULL
model_metadata <- readRDS("data/model_metadata.rds")
model_gene_list <- model_metadata$model_gene_list
model_protein_list <- model_metadata$model_protein_list

alias_table_all <- data.table::fread("data/gene_alias_table.txt")
hgnc <- data.table::fread("data/hgnc_complete_set.txt")

# Load interpret module's dropdown data first
interpret_module <- readRDS("data/interpret_metadata_small.rds")
choices_protein <- interpret_module$choices_protein
choices_celltype <- interpret_module$choices_celltype
choices_geneset <- interpret_module$choices_geneset

# Load interpret module large data later
protein_abundance_pbmc <- NULL
cite_rna_protein_pair <- NULL
merged_gradient_df <- NULL
merged_gsea_df <- NULL

#---------------------------------------------------
# Reference object (lazy loading)
#---------------------------------------------------
reference <- NULL
get_reference <- function() {
  if (is.null(reference)) {
    message("Loading PBMC reference...")
    reference <<- readRDS(
      "data/pbmc_multimodal_2023.rds"
    )
    message("PBMC reference loaded.")
  }
  return(reference)
}

# Load help page and privacy info
helpInfo <- read.delim("www/Figure_explain.csv", sep = ",")

# Define Shiny server (mutation-related logic removed)
shinyServer(function(input, output, session) {
  
  # Measure app startup time
  session$onFlushed(function() {
    
    ready_time <- Sys.time()
    
    startup_seconds <- round(
      as.numeric(difftime(
        ready_time,
        app_start_time,
        units = "secs"
      )),
      2
    )
    
    print(paste(
      "App became ready in",
      startup_seconds,
      "seconds"
    ))
    
  }, once = TRUE)
  
  
  #### Privacy policy ####
  observeEvent(input$privacytag1, {
    shiny::showModal(shiny::modalDialog(
      size = "l",
      includeHTML("www/PrivacyPolicy.html"),
      easyClose = TRUE
    ))
  })
  
  ######help tab content#####
  #help page for input data format
  
  observeEvent(input$helptab2p1
               , {
                 showModal(modalDialog(
                   size = "l",
                   includeHTML("www/helptab2_ExploreModel.html"),
                   easyClose = TRUE
                 ))
               })
  
  observeEvent(input$helptab3p1
               , {
                 showModal(modalDialog(
                   size = "l",
                   includeHTML("www/helptab3_PredictProteins.html"),
                   
                   easyClose = TRUE,
                 ))
               })
  
  observeEvent(input$helptab3p2
               , {
                 showModal(modalDialog(
                   size = "l",
                   includeHTML("www/helptab3.html"),
                   
                   easyClose = TRUE,
                 ))
               })
  
  ######help tab content#####
  
  observeEvent(input$download_prediction_metadata_check, {
    showModal(modalDialog(
      title = "Prediction not ready",
      HTML("Please run the prediction first. The download file will be available after the prediction and UMAPs are completed."),
      easyClose = TRUE
    ))
  })
  
  observeEvent(input$help1_1, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$panel_description[1]),
      easyClose = TRUE
    ))
  })
  observeEvent(input$help1_2, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$panel_description[2]),
      easyClose = TRUE
    ))
  })
  observeEvent(input$help1_3, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$panel_description[3]),
      easyClose = TRUE
    ))
  })
  observeEvent(input$help1_4, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$panel_description[4]),
      easyClose = TRUE
    ))
  })
  observeEvent(input$help2_1, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$panel_description[5]),
      easyClose = TRUE
    ))
  })
  observeEvent(input$help2_2, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$panel_description[6]),
      easyClose = TRUE
    ))
  })
  observeEvent(input$help2_3, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$panel_description[7]),
      easyClose = TRUE
    ))
  })
  observeEvent(input$help2_4, {
    showModal(modalDialog(
      size = "l",
      HTML(helpInfo$panel_description[8]),
      easyClose = TRUE
    ))
  })
  
  #------------------------------------------------
  #------------------------------------------------
  ##### Predict Protein #####
  #------------------------------------------------
  #------------------------------------------------

  # track data from example or user
  dat.source <- reactiveValues(source = NULL)
  run_trigger <- reactiveVal(0)
  processing_flag <- reactiveVal(FALSE)
  
  # track prediction progress
  processing_stage <- reactiveVal("Waiting to start...")
  processing_percent <- reactiveVal(0)
  
  output$processing_stage <- renderText({
    processing_stage()
  })
  
  output$processing_percent <- renderText({
    paste0(processing_percent(), "%")
  })

  update_processing <- function(stage, percent) {
    processing_stage(stage)
    processing_percent(percent)
    
    shinyjs::runjs(
      paste0(
        # "$('#prediction_stage_text').text('", stage, "');",
        "$('#prediction_stage_text').html('", stage, "');",
        "$('#prediction_progress_bar').css('width', '", percent, "%');",
        "$('#prediction_percent_text').text('", percent, "%');"
      )
    )
  }
  
  
  # Store prediction results
  pred_results <- reactiveValues(data = NULL)
  
  observeEvent(input$submitexample.p, {
    if (processing_flag()) return()
    
    showModal(modalDialog(
      size = "m",
      title = "Run example data?",
      HTML("The example PBMC dataset will be processed. This may take around 2 minutes."),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_run_example", "Run Example", class = "btn btn-primary")
      ),
      easyClose = TRUE
    ))
  })
  
  observeEvent(input$confirm_run_example, {
    removeModal()
    
    processing_flag(TRUE)
    
    processing_stage("Starting example analysis...")
    processing_percent(0)
    
    showModal(modalDialog(
      size = "m",
      title = "Processing example data",
      tagList(
        div(
          id = "prediction_stage_text",
          style = "font-size:16px; margin-bottom:12px; color:#333;",
          "Starting analysis..."
        ),
        div(
          style = "width:100%; background:#e9ecef; border-radius:20px; height:18px;",
          div(
            id = "prediction_progress_bar",
            style = "width:0%; background:teal; height:18px; border-radius:20px;"
          )
        ),
        p(
          id = "prediction_percent_text",
          style = "text-align:center; margin-top:8px;",
          "0%"
        )
      ),
      easyClose = FALSE,
      footer = NULL
    ))
    
    dat.source$source <- "example"
    run_trigger(run_trigger() + 1)
  })
  
  observeEvent(input$Reset.p, {
    # Clear file input field
    shinyjs::reset("exp_data.p")
    
    # Clear reactive values
    dat.source$source <- NULL
    saved_results$seurat <- NULL
    saved_results$prediction <- NULL
    saved_results$umap_clusters <- NULL
    saved_results$umap_celltype1 <- NULL
    saved_results$umap_celltype2 <- NULL
    pred_results$data <- NULL
    
    showModal(modalDialog(
      title = "Reset complete",
      HTML("All inputs and results have been cleared. You can now upload new data or rerun the example."),
      easyClose = TRUE
    ))
  })
  
    #Message for submit
  # observeEvent(input$submitmodel.p, {
  #   if (processing_flag()) return()
  #   
  #   if (is.null(input$exp_data.p)) {
  #     showModal(modalDialog(
  #       size = "m",
  #       title = "No file uploaded",
  #       HTML("Please upload your data before clicking submit."),
  #       easyClose = TRUE
  #     ))
  #   } else {
  #     processing_flag(TRUE)
  #     dat.source$source <- "user"
  #     run_trigger(run_trigger() + 1)
  #     showModal(modalDialog(
  #       size = "m",
  #       title = "Processing uploaded data",
  #       HTML("Your data is being processed. This window will close automatically when the prediction and UMAPs are ready."),
  #       easyClose = FALSE,
  #       footer = NULL
  #     ))
  #   }
  # })
  
  observeEvent(input$submitmodel.p, {
    if (processing_flag()) return()
    
    if (is.null(input$exp_data.p)) {
      showModal(modalDialog(
        size = "m",
        title = "No file uploaded",
        HTML("Please upload your data before clicking GO."),
        easyClose = TRUE
      ))
    } else {
      showModal(modalDialog(
        size = "m",
        title = "Start analysis?",
        HTML("Your dataset will be processed. Large datasets may take several minutes. See the Help page for estimated runtimes."),
        footer = tagList(
          modalButton("Cancel"),
          actionButton("confirm_run_uploaded", "Run Analysis", class = "btn btn-primary")
        ),
        easyClose = TRUE
      ))
    }
  })
  
  observeEvent(input$confirm_run_uploaded, {
    removeModal()
    
    processing_flag(TRUE)
    
    
    # showModal(modalDialog(
    #   size = "m",
    #   title = "Processing uploaded data",
    #   HTML("Your data is being processed. This window will close automatically when the prediction and UMAPs are ready."),
    #   easyClose = FALSE,
    #   footer = NULL
    # ))
    
    processing_stage("Starting uploaded data analysis...")
    processing_percent(0)
    
    showModal(modalDialog(
      size = "m",
      title = "Processing uploaded data",
      tagList(
        div(
          id = "prediction_stage_text",
          style = "font-size:16px; margin-bottom:12px; color:#333;",
          "Starting analysis..."
        ),
        div(
          style = "width:100%; background:#e9ecef; border-radius:20px; height:18px;",
          div(
            id = "prediction_progress_bar",
            style = "width:0%; background:teal; height:18px; border-radius:20px;"
          )
        ),
        p(
          id = "prediction_percent_text",
          style = "text-align:center; margin-top:8px;",
          "0%"
        )
      ),
      easyClose = FALSE,
      footer = NULL
    ))
    dat.source$source <- "user"
    run_trigger(run_trigger() + 1)
  })
  
  # Returns Seurat object based on source
  # rna_seurat_input <- reactive({
  #   req(dat.source$source)
  #   
  #   if (dat.source$source == "user") {
  #     req(input$exp_data.p)
  #     if (!file.exists(input$exp_data.p$datapath)) return(NULL)
  #     load_rna_file(input$exp_data.p$datapath)
  #   } else if (dat.source$source == "example") {
  #     LoadH5Seurat("Predictions/samples/gene_data_subset_500_cells.h5seurat")
  #   } else {
  #     NULL
  #   }
  # })
  
  rna_seurat_input <- reactive({
    req(dat.source$source)
    
    if (dat.source$source == "user") {
      req(input$exp_data.p)
      if (!file.exists(input$exp_data.p$datapath)) return(NULL)
      
      obj <- load_rna_file(
        file_path = input$exp_data.p$datapath,
        original_filename = input$exp_data.p$name
      )
      
    } else if (dat.source$source == "example") {
      obj <- LoadH5Seurat("Predictions/samples/gene_data_subset_500_cells.h5seurat")
      
    } else {
      return(NULL)
    }
    
    obj <- tryCatch({
      validate_dataset_size(seurat_obj = obj, max_cells = 50000, auto_downsample = TRUE)
    }, error = function(e) {
      stop(paste("Dataset validation failed:", e$message))
    })
    
    
    return(obj)
  })
  
  saved_results <- reactiveValues(
    seurat = NULL,
    prediction = NULL,
    prediction_metadata = NULL,
    runtime = NULL,
    gene_match_summary = NULL,
    umap_clusters = NULL,
    umap_celltype1 = NULL,
    umap_celltype2 = NULL
  )
  
  output$prediction_ready <- reactive({
    !is.null(saved_results$prediction_metadata)
  })
  
  outputOptions(output, "prediction_ready", suspendWhenHidden = FALSE)
  
  # Full prediction and Seurat workflow
  prediction_matrix <- reactive({
    
    if (run_trigger() == 0) {
      return(NULL)
    }
    
    tryCatch({t_start <- Sys.time()
    
    # req(run_trigger())
    req(rna_seurat_input())
    
    # Load model
    if (is.null(trained_model)) {
      trained_model <<- load_model_hdf5(
        "Predictions/CNN_X_magic_Y_protein_norm_renorm_NOgeneset_model_3.h5"
      )
      print("Model loaded!!!")
    }
    
    update_processing("Mapping uploaded genes to DeepGxP model genes...", 8)
    
    rna_mapped <- map_counts_to_model_genes(
      seurat_obj = rna_seurat_input(),
      model_gene_list = model_gene_list,
      alias_table_all = alias_table_all,
      hgnc = hgnc
    )
    
    gene_summary <- attr(rna_mapped, "gene_match_summary")
    
    saved_results$gene_match_summary <- paste0(
      format(gene_summary$matched_model_genes, big.mark = ","),
      " of ",
      format(gene_summary$total_model_genes, big.mark = ","),
      " model genes were detected after gene-symbol mapping. ",
      "Exact matches: ",
      format(gene_summary$exact, big.mark = ","),
      "; alias matches: ",
      format(gene_summary$alias, big.mark = ","),
      ". Missing genes were filled with zeros."
    )
    
    update_processing("Preprocessing scRNA-seq data...", 10)
    
    rna_seurat <- SCTransform(
      rna_mapped,
      assay = "RNA",
      new.assay.name = "SCT",
      verbose = FALSE
    )
    
    magic_imputed <- magic(
      t(GetAssayData(rna_seurat, assay = "SCT", layer = "data")),
      n.jobs = -1
    )
    
    magic_df <- as.data.frame(magic_imputed)
    
    update_processing("Predicting 224 surface proteins...", 20)
    
    prediction_df <- make_predictions(
      trained_model,
      magic_df,
      model_gene_list,
      model_protein_list
    )
    
    
    # Remove large data 
    rm(magic_imputed, magic_df)
    gc()
    
    prediction_df <- prediction_df[, order(names(prediction_df))]
    prediction_df <- t(prediction_df)
    
    # Add predicted protein
    rna_seurat[["ADT"]] <- CreateAssayObject(counts = prediction_df)
    
    update_processing("Clustering cells...", 40)
    
    # RNA analysis
    #rna_seurat <- FindVariableFeatures(rna_seurat, assay = "SCT", selection.method = "vst", nfeatures = 2000)
    rna_seurat <- RunPCA(rna_seurat, assay = "SCT", reduction.name = "pca", verbose = FALSE)
    rna_seurat <- RunUMAP(rna_seurat, reduction = "pca", dims = 1:20, reduction.name = "rna.umap")
    rna_seurat <- FindNeighbors(rna_seurat, reduction = "pca", dims = 1:20)
    rna_seurat <- FindClusters(rna_seurat, graph.name = "SCT_snn")
    
    # Protein analysis
    rna_seurat <- NormalizeData(rna_seurat, assay = "ADT", normalization.method = "CLR", margin = 2)
    rna_seurat <- ScaleData(rna_seurat, assay = "ADT", features = rownames(rna_seurat), verbose = FALSE)
    rna_seurat <- RunPCA(rna_seurat, assay = "ADT", features = rownames(rna_seurat), reduction.name = "apca", verbose = FALSE)
    rna_seurat <- RunUMAP(rna_seurat, reduction = "apca", dims = 1:10, reduction.name = "adt.umap")
    rna_seurat <- FindNeighbors(rna_seurat, reduction = "apca", dims = 1:10)
    rna_seurat <- FindClusters(rna_seurat, graph.name = "ADT_snn")
    
    # Combined WNN
    rna_seurat <- FindMultiModalNeighbors(rna_seurat,
                                          reduction.list = list("pca", "apca"),
                                          dims.list = list(1:20, 1:10),
                                          modality.weight.name = c("SCT.weight", "ADT.weight")
    )
    rna_seurat <- RunUMAP(rna_seurat, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
    rna_seurat <- FindClusters(rna_seurat, graph.name = "wsnn")
    
    update_processing("Annotating cell types...", 75)
    
    # Cell type label transfer
    reference_obj <- tryCatch({
      get_reference()
    }, error = function(e) {
      stop(paste("Failed to load reference:", e$message))
    })
    anchors <- FindTransferAnchors(reference = get_reference(), query = rna_seurat,
                                   normalization.method = "SCT", reference.reduction = "spca", dims = 1:50)
    rna_seurat <- TransferData(anchorset = anchors, reference = get_reference(), query = rna_seurat,
                               refdata = list(
                                 celltype.l1 = "celltype.l1",
                                 celltype.l2 = "celltype.l2",
                                 predicted_ADT = "ADT"
                               )
    )
    
    # rm(reference, anchors)
    rm(anchors)
    gc()
    
    update_processing("Generating visualizations...", 90)
    
    # Generate plots
    umap_clusters <- DimPlot(
      rna_seurat, reduction = 'adt.umap', group.by = "ADT_snn_res.0.8", label = TRUE, repel = TRUE, label.size = 2.5) + ggtitle("Protein") |
      DimPlot(rna_seurat, reduction = 'wnn.umap', group.by = "wsnn_res.0.8", label = TRUE, repel = TRUE, label.size = 2.5) + ggtitle("Combined") |
      DimPlot(rna_seurat, reduction = 'rna.umap', group.by = "SCT_snn_res.0.8", label = TRUE, repel = TRUE, label.size = 2.5) + ggtitle("RNA")
    
    # umap_celltype1 <- DimPlot(rna_seurat, reduction = 'adt.umap', group.by = "predicted.celltype.l1", label = TRUE, repel = TRUE, label.size = 2.5) + ggtitle("Protein") |
    #   DimPlot(rna_seurat, reduction = 'wnn.umap', group.by = "predicted.celltype.l1", label = TRUE, repel = TRUE, label.size = 2.5) + ggtitle("Combined") |
    #   DimPlot(rna_seurat, reduction = 'rna.umap', group.by = "predicted.celltype.l1", label = TRUE, repel = TRUE, label.size = 2.5) + ggtitle("RNA")
    
    p_adt_l1 <- DimPlot(
      rna_seurat,
      reduction = "adt.umap",
      group.by = "predicted.celltype.l1",
      label = TRUE,
      repel = TRUE,
      label.size = 2.5
    ) +
      ggtitle("Protein") +
      NoLegend()
    
    p_wnn_l1 <- DimPlot(
      rna_seurat,
      reduction = "wnn.umap",
      group.by = "predicted.celltype.l1",
      label = TRUE,
      repel = TRUE,
      label.size = 2.5
    ) +
      ggtitle("Combined") +
      NoLegend()
    
    p_rna_l1 <- DimPlot(
      rna_seurat,
      reduction = "rna.umap",
      group.by = "predicted.celltype.l1",
      label = TRUE,
      repel = TRUE,
      label.size = 2.5
    ) +
      ggtitle("RNA") +
      theme(
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 10)
      )
    
    umap_celltype1 <- p_adt_l1 | p_wnn_l1 | p_rna_l1
    
    celltype_l2_colors <- c(
      "B memory" = "#1f77b4",
      "B naive" = "#ff7f0e",
      "CD14 Mono" = "#2ca02c",
      "CD16 Mono" = "#d62728",
      "CD4 CTL" = "#9467bd",
      "CD4 Naive" = "#8c564b",
      "CD4 Proliferating" = "#e377c2",
      "CD4 TCM" = "#7f7f7f",
      "CD4 TEM" = "#bcbd22",
      "CD8 Naive" = "#17becf",
      "CD8 TCM" = "#aec7e8",
      "CD8 TEM" = "#ffbb78",
      "dnT" = "#98df8a",
      "gdT" = "#ff9896",
      "MAIT" = "#c5b0d5",
      "NK" = "#c49c94",
      "pDC" = "#f7b6d2",
      "Plasmablast" = "#c7c7c7",
      "Treg" = "#dbdb8d"
    )
    
    p_adt_l2 <- DimPlot(
      rna_seurat,
      reduction = "adt.umap",
      group.by = "predicted.celltype.l2",
      label = TRUE,
      repel = TRUE,
      label.size = 2.5,
      cols = celltype_l2_colors
    ) +
      ggtitle("Protein") +
      NoLegend()
    
    p_wnn_l2 <- DimPlot(
      rna_seurat,
      reduction = "wnn.umap",
      group.by = "predicted.celltype.l2",
      label = TRUE,
      repel = TRUE,
      label.size = 2.5,
      cols = celltype_l2_colors
    ) +
      ggtitle("Combined") +
      NoLegend()
    
    p_rna_l2 <- DimPlot(
      rna_seurat,
      reduction = "rna.umap",
      group.by = "predicted.celltype.l2",
      label = TRUE,
      repel = TRUE,
      label.size = 2.5,
      cols = celltype_l2_colors
    ) +
      ggtitle("RNA") +
      guides(color = guide_legend(ncol = 2, override.aes = list(size = 4))) +
      theme(
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 9)
      )
    
    umap_celltype2 <- p_adt_l2 | p_wnn_l2 | p_rna_l2
    
    prediction_metadata <- as.data.frame(t(prediction_df))
    prediction_metadata$Cell_ID <- rownames(prediction_metadata)
    
    meta_df <- rna_seurat@meta.data
    
    prediction_metadata$RNA_cluster <- meta_df[rownames(prediction_metadata), "SCT_snn_res.0.8"]
    prediction_metadata$Protein_cluster <- meta_df[rownames(prediction_metadata), "ADT_snn_res.0.8"]
    prediction_metadata$WNN_cluster <- meta_df[rownames(prediction_metadata), "wsnn_res.0.8"]
    prediction_metadata$Celltype_Level1 <- meta_df[rownames(prediction_metadata), "predicted.celltype.l1"]
    prediction_metadata$Celltype_Level2 <- meta_df[rownames(prediction_metadata), "predicted.celltype.l2"]
    
    prediction_metadata <- prediction_metadata[, c(
      "Cell_ID",
      setdiff(colnames(prediction_metadata), "Cell_ID")
    )]
    
    update_processing("Finalizing results...", 100)
    
    # Save for later
    #saved_results$seurat <- rna_seurat
    saved_results$prediction <- prediction_df
    saved_results$prediction_metadata <- prediction_metadata
    saved_results$umap_clusters <- umap_clusters
    saved_results$umap_celltype1 <- umap_celltype1
    saved_results$umap_celltype2 <- umap_celltype2
    
    # Remove large object after plots are created
    rm(rna_seurat)
    gc()
    
    # ----------------------------------------
    # Evaluate prediction time (START)
    #-----------------------------------------
    t_end <- Sys.time()
    
    runtime_minutes <- round(
      as.numeric(difftime(
        t_end,
        t_start,
        units = "mins"
      )),
      2
    )
    
    saved_results$runtime <- runtime_minutes
    
    print(paste(
      "Prediction workflow completed in",
      runtime_minutes,
      "minutes"
    ))
    
    
    # ----------------------------------------
    # Evaluate prediction time (END)
    #-----------------------------------------
    processing_flag(FALSE)
    removeModal()
    
    prediction_df
    
    }, error = function(e) {
      
      processing_stage("Analysis failed")
      processing_percent(0)
      
      processing_flag(FALSE)
      run_trigger(0)
      removeModal()
      
      showModal(
        modalDialog(
          title = "Analysis failed",
          HTML(paste0("<b>Error:</b><br>", e$message)),
          easyClose = TRUE
        )
      )
      
      return(NULL)
      
    })
  })

  
  # output$prediction_table <- renderDataTable({
  #   req(prediction_matrix())
  #   pred_results$data <- prediction_matrix()  # store for export
  #   DT::datatable(
  #     pred_results$data,
  #     extensions = c("Scroller", "Buttons"),
  #     rownames = TRUE,
  #     options = list(
  #       scrollX = TRUE,
  #       paging = TRUE,
  #       dom = "frtip",
  #       lengthMenu = list(c(10, 15, 20), c(10, 15, 20)),
  #       pageLength = 10
  #     ),
  #     selection = list(mode = "single")
  #   )
  # }, server = TRUE)
  
  output$prediction_table <- renderDataTable({
    req(prediction_matrix())
    pred_results$data <- prediction_matrix()
    preview_df <- pred_results$data[, 1:min(10, ncol(pred_results$data))]
    DT::datatable(
      preview_df,
      extensions = c("Scroller", "Buttons"),
      rownames = TRUE,
      options = list(
        scrollX = TRUE,
        paging = TRUE,
        dom = "frtip",
        lengthMenu = list(c(10, 15, 20), c(10, 15, 20)),
        pageLength = 10
      )
    )
  }, server = TRUE)
  
  
  output$data_source_info <- renderText({
    if (is.null(dat.source$source)) return("No data submitted yet.")
    if (dat.source$source == "user") return("User uploaded data.")
    if (dat.source$source == "example") return("Using example PBMC data.")
  })
  
  output$runtime_text <- renderText({
    req(saved_results$runtime)
    
    paste(
      "Analysis completed in",
      saved_results$runtime,
      "minutes."
    )
  })
  
  output$gene_match_text <- renderUI({
    
    req(saved_results$gene_match_summary)
    
    HTML(
      paste0(
        "<div style='margin-top:8px; color:#444;'>",
        saved_results$gene_match_summary,
        "</div>"
      )
    )
  })
  
  output$prediction_summary <- renderUI({
    req(saved_results$prediction_metadata)
    req(saved_results$runtime)
    req(saved_results$gene_match_summary)
    
    n_cells <- nrow(saved_results$prediction_metadata)
    n_proteins <- length(model_protein_list)
    
    HTML(paste0(
      "<div style='background:white; border-radius:18px; padding:22px 26px;
                box-shadow:0 4px 18px rgba(0,0,0,0.08); color:#333;'>",
      
      "<h4 style='color:teal; margin-bottom:14px;'>
      Prediction completed successfully
    </h4>",
      
      "<p style='margin-bottom:6px;'><b>Runtime:</b> ",
      saved_results$runtime, " minutes</p>",
      
      "<p style='margin-bottom:6px;'><b>Cells processed:</b> ",
      format(n_cells, big.mark = ","), "</p>",
      
      "<p style='margin-bottom:6px;'><b>Proteins predicted:</b> ",
      format(n_proteins, big.mark = ","), "</p>",
      
      "<p style='margin-bottom:0px;'><b>Gene matching:</b> ",
      saved_results$gene_match_summary, "</p>",
      
      "</div>"
    ))
  })
  
  # UMAPS
  # Cluster-based UMAP plot
  output$umap_clusters <- renderPlot({
    req(saved_results$umap_clusters)
    saved_results$umap_clusters
  })
  
  # Cell type level 1 UMAP plot
  output$umap_celltype1 <- renderPlot({
    req(saved_results$umap_celltype1)
    saved_results$umap_celltype1
  })
  
  # Cell type level 2 UMAP plot
  output$umap_celltype2 <- renderPlot({
    req(saved_results$umap_celltype2)
    saved_results$umap_celltype2
  })
  
  # output$download_prediction_txt <- downloadHandler(
  #   filename = function() "prediction_matrix.txt",
  #   content = function(file) {
  #     req(saved_results$prediction)
  #     write.table(saved_results$prediction, file, sep = "\t", quote = FALSE)
  #   }
  # )
  
  # output$download_prediction_metadata <- downloadHandler(
  #   filename = function() {
  #     
  #     if (is.null(saved_results$prediction_metadata)) {
  #       return("prediction_not_ready.txt")
  #     }
  #     
  #     "predicted_proteins_with_clusters_celltypes.txt"
  #   },
  #   
  #   content = function(file) {
  #     
  #     if (is.null(saved_results$prediction_metadata)) {
  #       
  #       writeLines(
  #         "Prediction has not been run yet.",
  #         con = file
  #       )
  #       
  #       showModal(modalDialog(
  #         title = "Prediction not ready",
  #         HTML("Please run the prediction first before downloading."),
  #         easyClose = TRUE
  #       ))
  #       
  #       return()
  #     }
  #     
  #     write.table(
  #       saved_results$prediction_metadata,
  #       file,
  #       sep = "\t",
  #       quote = FALSE,
  #       row.names = FALSE
  #     )
  #   }
  # )
  
  output$download_prediction_metadata <- downloadHandler(
    filename = function() {
      "Predictions.txt"
    },
    content = function(file) {
      req(saved_results$prediction_metadata)
      
      write.table(
        saved_results$prediction_metadata,
        file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
    }
  )
  
  # output$download_clusters_pdf <- downloadHandler(
  #   filename = function() "plot_clusters.pdf",
  #   content = function(file) {
  #     ggsave(file, plot = saved_results$umap_clusters, width = 10, height = 4)
  #   }
  # )
  # 
  # output$download_celltype1_pdf <- downloadHandler(
  #   filename = function() "plot_celltype1.pdf",
  #   content = function(file) {
  #     ggsave(file, plot = saved_results$umap_celltype1, width = 10, height = 4)
  #   }
  # )
  # 
  # output$download_celltype2_pdf <- downloadHandler(
  #   filename = function() "plot_celltype2.pdf",
  #   content = function(file) {
  #     ggsave(file, plot = saved_results$umap_celltype2, width = 10, height = 4)
  #   }
  # )
  
  output$download_all_umaps <- downloadHandler(
    filename = function() {
      "UMAP_plots.zip"
    },
    content = function(file) {
      req(saved_results$umap_clusters)
      req(saved_results$umap_celltype1)
      req(saved_results$umap_celltype2)
      
      temp_dir <- tempdir()
      
      clusters_pdf <- file.path(temp_dir, "UMAP_clusters.pdf")
      celltype1_pdf <- file.path(temp_dir, "UMAP_celltype_level1.pdf")
      celltype2_pdf <- file.path(temp_dir, "UMAP_celltype_level2.pdf")
      
      ggsave(clusters_pdf, plot = saved_results$umap_clusters, width = 10, height = 4)
      ggsave(celltype1_pdf, plot = saved_results$umap_celltype1, width = 10, height = 4)
      ggsave(celltype2_pdf, plot = saved_results$umap_celltype2, width = 10, height = 4)
      
      zip::zipr(
        zipfile = file,
        files = c(clusters_pdf, celltype1_pdf, celltype2_pdf)
      )
    }
  )
  
  #------------------------------------------------
  #------------------------------------------------
  ##### Model Interpret #####
  #------------------------------------------------
  #------------------------------------------------
  
  load_interpret_data <- function() {
    
    if (is.null(protein_abundance_pbmc) || is.null(cite_rna_protein_pair)) {
      message("Loading interpret heavy RDS data...")
      interpret_heavy_data <- readRDS("data/interpret_heavy_data.rds")
      
      protein_abundance_pbmc <<- interpret_heavy_data$protein_abundance_pbmc
      cite_rna_protein_pair <<- interpret_heavy_data$cite_rna_protein_pair
    }
    
    if (is.null(merged_gradient_df)) {
      message("Loading gradient parquet...")
      merged_gradient_df <<- arrow::read_parquet(
        "data/merged_gradientRank_df.parquet"
      )
    }
    
    if (is.null(merged_gsea_df)) {
      message("Loading enrichment parquet...")
      merged_gsea_df <<- arrow::read_parquet(
        "data/merged_enrichmentORA.parquet"
      )
    }
  }
  
  # Update dropdowns on startup
  print("Initializing cell_type choices")
  #print(choices_celltype)
  updateSelectizeInput(session = session,
                       inputId = "cell_type",
                       choices = c(NA, choices_celltype),
                       selected = NULL,
                       options = list(placeholder ="e.g. B", maxItems = 1),
                       server = TRUE)
  #Gene selection#
  print("Initializing protein_name choices")
  #print(choices_protein)
  updateSelectizeInput(session = session,
                       inputId = "protein_name",
                       choices = c(NA, choices_protein),
                       selected = NULL,
                       options = list(placeholder = "e.g. CD20", maxItems = 1),
                       server = TRUE)
  
  #Gene selection#
  print("Initializing geneset choices")
  #print(choices_geneset)
  updateSelectizeInput(session = session,
                       inputId = "geneset_name",
                       choices = c(NA, choices_geneset),
                       selected = NULL,
                       options = list(placeholder = "e.g. GOBP", maxItems = 1),
                       server = TRUE)
  
  # observeEvent(input$cell_type, {load_interpret_data()}, ignoreInit = TRUE)
  # observeEvent(input$protein_name, {load_interpret_data()}, ignoreInit = TRUE)
  # observeEvent(input$geneset_name, {load_interpret_data()}, ignoreInit = TRUE)
  
  interpret_status <- reactiveVal(
    "Please select a cell type and protein, then click Run Analysis."
  )
  
  interpret_processing <- reactiveVal(FALSE)
  
  interpret_selection_changed <- reactiveVal(FALSE)
  
  # interpret_trigger <- eventReactive(input$run_interpret, {
  #   req(input$cell_type, input$protein_name)
  #   load_interpret_data()
  #   interpret_selection_changed(FALSE)
  #   interpret_status(
  #     paste(
  #       "Selected cell type:", input$cell_type,
  #       "| Selected protein:", input$protein_name
  #     )
  #   )
  #   list(cell_type = input$cell_type,
  #        protein_name = input$protein_name)
  # })
  
  # interpret_trigger <- eventReactive(input$run_interpret, {
  #   req(input$cell_type, input$protein_name)
  #   
  #   # showModal(modalDialog(
  #   #   size = "m",
  #   #   title = "Loading Explore Model results",
  #   #   tagList(
  #   #     div(
  #   #       style = "font-size:16px; color:#333; margin-bottom:10px;",
  #   #       "Retrieving protein abundance and gene importance results for the selected cell type and protein."
  #   #     ),
  #   #     div(
  #   #       style = "font-size:14px; color:#666;",
  #   #       "This window will close automatically when the plots are ready."
  #   #     )
  #   #   ),
  #   #   easyClose = FALSE,
  #   #   footer = NULL
  #   # ))
  #   # 
  #   # on.exit(removeModal(), add = TRUE)
  #   
  #   load_interpret_data()
  #   
  #   interpret_selection_changed(FALSE)
  #   
  #   # interpret_status(
  #   #   paste(
  #   #     "Selected cell type:", input$cell_type,
  #   #     "| Selected protein:", input$protein_name
  #   #   )
  #   # )
  #   
  #   interpret_status(
  #     "Loading protein abundance and gene importance results..."
  #   )
  #   
  #   list(
  #     cell_type = input$cell_type,
  #     protein_name = input$protein_name
  #   )
  # })
  
  interpret_trigger <- eventReactive(input$run_interpret, {
    
    # check if interpret data still not ready
    if (is.null(choices_protein) || is.null(choices_celltype)) {
      
      showModal(modalDialog(
        title = "Application initializing",
        HTML("The Explore Model module is still loading required data. Please wait a few moments and try again."),
        easyClose = TRUE
      ))
      
      return(NULL)
    }
    
    req(input$cell_type, input$protein_name)
    
    interpret_processing(TRUE)
    
    showModal(modalDialog(
      size = "m",
      title = "Preparing Explore Model results",
      tagList(
        
        div(
          style = "font-size:16px; color:#333; margin-bottom:12px;",
          "Retrieving protein abundance and gene importance information..."
        ),
      ),
      easyClose = FALSE,
      footer = NULL
    ))
    
    load_interpret_data()
    
    interpret_selection_changed(FALSE)
    
    interpret_status(
      paste(
        "Selected cell type:", input$cell_type,
        "| Selected protein:", input$protein_name
      )
    )
    
    list(
      cell_type = input$cell_type,
      protein_name = input$protein_name
    )
    
  })
  
  
  observeEvent(
    list(input$cell_type, input$protein_name),
    {interpret_selection_changed(TRUE)
      interpret_status("Selection changed. Please click Run Analysis to update results.")
      updateSelectizeInput(
        session = session,
        inputId = "geneset_name",
        selected = character(0)
      )
      }, ignoreInit = TRUE)
  
  
  output$text.gradient <- renderText({
    interpret_status()
  })
  
  
  # --- Reactive filtered data ---
  output$protein_violin_plot <- renderPlotly({
    # req(input$cell_type, input$protein_name)
    selected <- interpret_trigger()
    
    plotViolin(df = protein_abundance_pbmc, 
               celltype = selected$cell_type, 
               protein = selected$protein_name)
  })

  filtered_gradient <- reactive({
    # req(input$cell_type, input$protein_name)
    selected <- interpret_trigger()

    df <- merged_gradient_df %>%
      filter(CellType == selected$cell_type, Protein == selected$protein_name)
    df$high_zscore <- round(df$high_zscore, 4)
    df <- df %>%
        select(GeneName, high_zscore, high_rank)
      
    df
  })

  
  # --- Output filtered table (optional for UI inspection) ---
  output$filtered_gradient_table <- renderDataTable({
    req(filtered_gradient())
    datatable(filtered_gradient(), options = list(scrollX = TRUE, pageLength = 10))
  })
  
  # output$gradient_density_plot <- renderPlotly({
  #   selected <- interpret_trigger()
  #   table.out <- filtered_gradient()
  #   plotDensity2.0(df = table.out, 
  #                  selected_protein = selected$protein_name, 
  #                  match_pair = cite_rna_protein_pair)
  # })
  
  output$gradient_density_plot <- renderPlotly({
    
    selected <- interpret_trigger()
    
    table.out <- filtered_gradient()
    
    p <- plotDensity2.0(
      df = table.out,
      selected_protein = selected$protein_name,
      match_pair = cite_rna_protein_pair
    )
    
    if (interpret_processing()) {
      removeModal()
      interpret_processing(FALSE)
    }
    
    p
  })
  
  # --- filtered enrichment data ---
  filtered_enrichment <- reactive({
    selected <- interpret_trigger()
    req(input$geneset_name)
    
    df <- merged_gsea_df %>%
      filter(celltype == selected$cell_type,
             proName == selected$protein_name,
             Cluster == input$geneset_name )
    
    common_cols <- c(
      "Description", "GeneRatio", "BgRatio", "RichFactor", 
      "FoldEnrichment", "zScore", "pvalue", "p.adjust", "qvalue", 
      "geneID", "Count"
    )
    
    selected_cols <- if (input$geneset_name %in% c("GOBP", "GOCC", "GOMF")) {
      c("ID", common_cols)
    } else {
      common_cols
    }
    
    df_pos <- df %>%
      filter(direction == "pos") %>%
      select(all_of(selected_cols))
    
    df_neg <- df %>%
      filter(direction == "neg") %>%
      select(all_of(selected_cols))
    
    list(pos = df_pos, neg = df_neg)
  })
  
  
  # --- Output filtered table (optional for UI inspection) ---
  output$filtered_enrichment_table_pos <- renderDataTable({
    req(filtered_enrichment()$pos)
    if (is.null(merged_gsea_df)) {
      return(NULL)
    }
    datatable(filtered_enrichment()$pos, options = list(scrollX = TRUE, pageLength = 5)) %>%
      formatRound(columns = c('RichFactor', 'FoldEnrichment', "zScore"), digits = 3) %>%
      formatSignif(columns = c("pvalue", "p.adjust", "qvalue"), digits = 3)
  })
  
  output$filtered_enrichment_table_neg <- renderDataTable({
    req(filtered_enrichment()$neg)
    datatable(filtered_enrichment()$neg, options = list(scrollX = TRUE, pageLength = 5)) %>%
      formatRound(columns = c('RichFactor', 'FoldEnrichment', "zScore"), digits = 3) %>%
      formatSignif(columns = c("pvalue", "p.adjust", "qvalue"), digits = 3)
  })
  
  output$filtered_enrichment_plot_pos <- renderPlotly({
    df <- filtered_enrichment()$pos
    plotEnrichment(df)
  })
  output$filtered_enrichment_plot_neg <- renderPlotly({
    table.out <- filtered_enrichment()$neg
    plotEnrichment(table.out)
  })

})


