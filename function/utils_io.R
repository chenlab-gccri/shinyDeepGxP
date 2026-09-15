# utils_io.R for file handling

# Function to load RNA file and create Seurat object
# load_rna_file <- function(file_path) {
#   rna_data <- LoadH5Seurat(file_path)
#   CreateSeuratObject(counts = rna_data@assays$RNA@counts)
# }

load_rna_file <- function(file_path, original_filename = NULL) {
  if (!file.exists(file_path)) {
    stop("Uploaded file was not found on the server.")
  }
  
  file_label <- if (!is.null(original_filename)) original_filename else basename(file_path)
  file_label_lower <- tolower(file_label)
  
  get_ext <- function(x) {
    x <- tolower(x)
    if (grepl("\\.h5seurat$", x)) return("h5seurat")
    if (grepl("\\.rds$", x)) return("rds")
    if (grepl("\\.txt(\\.gz)?$", x)) return("txt")
    
    return(tools::file_ext(x))
  }
  
  make_seurat <- function(counts, project = "query") {
    if (is.list(counts)) {
      # Read10X returns a list for multi-modal data.
      # Use Gene Expression/RNA when available.
      if ("Gene Expression" %in% names(counts)) {
        counts <- counts[["Gene Expression"]]
      } else if ("RNA" %in% names(counts)) {
        counts <- counts[["RNA"]]
      } else {
        counts <- counts[[1]]
      }
    }
    Seurat::CreateSeuratObject(counts = counts, project = project)
  }
  
  ext <- get_ext(file_label_lower)
  
  if (ext == "h5seurat") {
    message("Loading h5Seurat file: ", file_label)
    
    # Shiny upload temp paths often lose the original extension.
    # Copy to a temporary path with the expected .h5seurat extension.
    tmp_dir <- tempfile("h5seurat_load_")
    dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
    h5seurat_path <- file.path(tmp_dir, "query.h5seurat")
    file.copy(file_path, h5seurat_path, overwrite = TRUE)
    
    obj <- SeuratDisk::LoadH5Seurat(h5seurat_path)
    return(obj)
  }
  
  if (ext == "rds") {
    obj <- readRDS(file_path)
    if (!inherits(obj, "Seurat")) {
      stop("Uploaded RDS file is not a Seurat object.")
    }
    return(obj)
  }
  
  if (ext == "txt") {
    message("Loading TXT gene-by-cell matrix: ", file_label)
    
    counts_df <- utils::read.table(
      file_path,
      sep = "\t",
      header = TRUE,
      row.names = 1,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      comment.char = "",
      quote = "\""
    )
    
    if (nrow(counts_df) == 0 || ncol(counts_df) == 0) {
      stop("The uploaded expression matrix is empty. Please upload a TXT gene-by-cell matrix with genes in rows and cells in columns.")
    }
    
    counts_df[] <- lapply(counts_df, function(x) {
      as.numeric(as.character(x))
    })
    
    if (anyNA(counts_df)) {
      stop("The uploaded TXT file contains non-numeric expression values. Please check that all expression columns contain numeric counts.")
    }
    
    counts_mat <- as.matrix(counts_df)
    rownames(counts_mat) <- rownames(counts_df)
    
    return(make_seurat(counts_mat, project = "txt_query"))
  }
  
  stop(paste0(
    "Unsupported input format: ", file_label,
    ". Supported formats are .h5Seurat, .rds, and .txt."
  ))
}

validate_dataset_size <- function(seurat_obj,
                                  max_cells = 10000,
                                  auto_downsample = TRUE) {
  
  n_cells <- ncol(seurat_obj)
  n_genes <- nrow(seurat_obj)
  
  message("Uploaded dataset: ", n_cells, " cells x ", n_genes, " genes")
  
  if (n_cells > max_cells) {
    
    if (auto_downsample) {
      
      warning(
        paste0(
          "Dataset contains ", n_cells,
          " cells. Randomly downsampling to ",
          max_cells,
          " cells for web-server processing."
        )
      )
      
      set.seed(123)
      
      selected_cells <- sample(
        colnames(seurat_obj),
        max_cells
      )
      
      seurat_obj <- subset(
        seurat_obj,
        cells = selected_cells
      )
      
      attr(seurat_obj, "downsample_message") <-
        paste0(
          "Input dataset contained ",
          n_cells,
          " cells and was downsampled to ",
          max_cells,
          " cells for web-server analysis."
        )
      
    } else {
      
      stop(
        paste0(
          "Dataset contains ",
          n_cells,
          " cells. The web server currently supports up to ",
          max_cells,
          " cells."
        )
      )
    }
  }
  
  return(seurat_obj)
}

map_counts_to_model_genes <- function(seurat_obj,
                                      model_gene_list,
                                      alias_table_all = NULL,
                                      hgnc = NULL) {
  
  counts_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  counts_mat <- as.matrix(counts_mat)
  # if (!inherits(counts_mat, "dgCMatrix")) {
  #   counts_mat <- as(counts_mat, "dgCMatrix")
  # }
  
  original_genes <- trimws(rownames(counts_mat))
  original_genes <- toupper(original_genes)
  
  # Fix common Excel-converted SEPT/MARCH genes
  original_genes <- gsub("^([0-9]+)-SEP$", "SEPT\\1", original_genes)
  original_genes <- gsub("^([0-9]+)-MAR$", "MARCH\\1", original_genes)
  
  rownames(counts_mat) <- original_genes
  model_gene_list <- toupper(model_gene_list)
  
  gene_map <- data.frame(
    input_gene = original_genes,
    model_gene = original_genes,
    match_type = ifelse(original_genes %in% model_gene_list, "exact", NA),
    stringsAsFactors = FALSE
  )
  
  # Optional alias mapping
  if (!is.null(alias_table_all) && !is.null(hgnc)) {
    
    approved_symbols <- unique(toupper(hgnc$symbol))
    
    alias_long <- alias_table_all[
      ,
      .(alias_gene = unlist(strsplit(aliases, "\\|"))),
      by = model_gene
    ]
    
    alias_long[, model_gene := toupper(trimws(model_gene))]
    alias_long[, alias_gene := toupper(trimws(alias_gene))]
    alias_long <- alias_long[alias_gene != "" & !is.na(alias_gene)]
    
    # protect approved HGNC symbols from being wrongly reassigned
    alias_long <- alias_long[
      !(alias_gene %in% approved_symbols & !(alias_gene %in% model_gene_list))
    ]
    
    alias_long <- unique(alias_long)
    
    idx <- is.na(gene_map$match_type) &
      gene_map$input_gene %in% alias_long$alias_gene
    
    gene_map$model_gene[idx] <- alias_long$model_gene[
      match(gene_map$input_gene[idx], alias_long$alias_gene)
    ]
    gene_map$match_type[idx] <- "alias"
  }
  
  keep <- gene_map$model_gene %in% model_gene_list
  counts_keep <- counts_mat[keep, , drop = FALSE]
  mapped_genes <- gene_map$model_gene[keep]
  
  # Collapse duplicated genes after alias mapping
  counts_collapsed <- rowsum(counts_keep, group = mapped_genes)
  
  # Fill missing model genes with zeros
  missing_genes <- setdiff(model_gene_list, rownames(counts_collapsed))
  zero_mat <- matrix(
    0,
    nrow = length(missing_genes),
    ncol = ncol(counts_collapsed),
    dimnames = list(missing_genes, colnames(counts_collapsed))
  )
  
  counts_model <- rbind(counts_collapsed, zero_mat)
  counts_model <- counts_model[model_gene_list, , drop = FALSE]
  
  mapped_obj <- CreateSeuratObject(counts = counts_model, meta.data = seurat_obj@meta.data)
  
  attr(mapped_obj, "gene_match_summary") <- list(
    exact = sum(gene_map$match_type == "exact", na.rm = TRUE),
    alias = sum(gene_map$match_type == "alias", na.rm = TRUE),
    matched_model_genes = sum(model_gene_list %in% rownames(counts_collapsed)),
    total_model_genes = length(model_gene_list),
    missing_model_genes = length(missing_genes)
  )
  
  return(mapped_obj)
}


# Function to handle predictions
make_predictions <- function(trained_model, magic_df, model_gene_list, model_protein_list) {
  missing_genes <- setdiff(model_gene_list, colnames(magic_df))
  magic_df[missing_genes] <- 0
  magic_df <- magic_df[, model_gene_list, drop = FALSE]
  x_test <- array(data.matrix(magic_df), dim = c(nrow(magic_df), ncol(magic_df), 1))
  prediction_matrix <- predict(trained_model, x_test)
  prediction_matrix <- round(prediction_matrix, 4)
  prediction_df <- data.frame(prediction_matrix, row.names = rownames(magic_df))
  colnames(prediction_df) <- model_protein_list
  return(prediction_df)
}

process_seurat_and_predict <- function(seurat_obj) {
  seurat_obj <- SCTransform(seurat_obj, assay = 'RNA', new.assay.name = 'SCT')
  magic_imputed <- magic(t(GetAssayData(seurat_obj, assay = "SCT", layer = "data")))
  magic_df <- as.data.frame(magic_imputed)
  prediction_df <- make_predictions(trained_model, magic_df, model_gene_list, model_protein_list)
  prediction_df <- prediction_df[, order(names(prediction_df))]
  t(prediction_df)  # Return transposed version directly
}

# plotDensity <- function(df){
#   # Compute density curve
#   dens <- density(df$high_zscore)
#   df_dens <- data.frame(x = dens$x, y = dens$y)
#   
#   # Assign color
#   df_sig <- df[df$high_zscore > 1.96 | df$high_zscore < -1.96, ]
#   df_sig$color <- factor(ifelse(df_sig$high_zscore > 1.96, "pink", "blue"),
#                          levels = c("blue", "pink"))
#   
#   # Start plot with density curve
#   p <- plot_ly(df_dens, x = ~x, y = ~y, type = 'scatter', mode = 'lines',
#                line = list(width = 2), name = "Density", hoverinfo = 'none')
#   
#   # Add significant genes with tooltips and color
#   p <- p %>% add_markers(data = df_sig, x = ~high_zscore, y = rep(0, nrow(df_sig)),
#                          marker = list(size = 6, color = df_sig$color),
#                          text = ~paste("Gene:", GeneName,
#                                        "<br>Z-score:", round(high_zscore, 2),
#                                        "<br>Rank:", high_rank),
#                          #"<br>Protein:", Protein,
#                          #"<br>CellType:", CellType),
#                          hoverinfo = 'text', name = "Significant Genes")
#   
#   # Layout
#   p <- p %>% layout(title = '',
#                xaxis = list(title = "normalized IG (nIG)"),
#                yaxis = list(title = "Density"))
#   
#   return(p)
# }

plotDensity2.0 <- function(df, selected_protein, match_pair){
  target_gene <- match_pair %>%
    filter(protName == selected_protein) %>%
    slice_max(corr, n = 1, with_ties = FALSE) %>%
    pull(geneName)
  
  # Compute density curve
  dens <- density(df$high_zscore)
  df_dens <- data.frame(x = dens$x, y = dens$y)
  
  # Assign color to significant genes
  # df_sig <- df %>%
  #   filter(high_zscore > 1.96 | high_zscore < -1.96) %>%
  #   mutate(color = factor(ifelse(high_zscore > 1.96, "pink", "blue"),
  #                         levels = c("blue", "pink")))
  
  df_sig_pos <- df %>%
    filter(high_zscore > 1.96)
  
  df_sig_neg <- df %>%
    filter(high_zscore < -1.96)
  
  # Subset target gene info if available
  if (length(target_gene) > 0 && target_gene %in% df$GeneName) {
    df_target <- df %>%
      filter(GeneName == target_gene) %>%
      mutate(color = case_when(
        high_zscore > 1.96 ~ "pink",
        high_zscore < -1.96 ~ "blue",
        TRUE ~ "gray"
      ))
    self_gene_label <- "Self-Gene"
  } else {
    df_target <- NULL
    self_gene_label <- "No Self-Gene Found"
  }
  
  # Start plot with density curve
  p <- plot_ly(df_dens, x = ~x, y = ~y, type = 'scatter', mode = 'lines',
               line = list(width = 2), name = "Density", hoverinfo = 'none')
  
  # Add significant gene points
  # p <- p %>% add_markers(data = df_sig, x = ~high_zscore, y = rep(0, nrow(df_sig)),
  #                        marker = list(size = 6, color = df_sig$color),
  #                        text = ~paste("Gene:", GeneName,
  #                                      "<br>Z-score:", round(high_zscore, 2),
  #                                      "<br>Rank:", high_rank),
  #                        hoverinfo = 'text', name = "Significant Genes")
  
  p <- p %>% add_markers(
    data = df_sig_pos,
    x = ~high_zscore,
    y = rep(0, nrow(df_sig_pos)),
    marker = list(size = 6, color = "pink"),
    text = ~paste(
      "Gene:", GeneName,
      "<br>Z-score:", round(high_zscore, 2),
      "<br>Rank:", high_rank
    ),
    hoverinfo = "text",
    name = "Positive Significant Genes"
  )
  
  p <- p %>% add_markers(
    data = df_sig_neg,
    x = ~high_zscore,
    y = rep(0, nrow(df_sig_neg)),
    marker = list(size = 6, color = "blue"),
    text = ~paste(
      "Gene:", GeneName,
      "<br>Z-score:", round(high_zscore, 2),
      "<br>Rank:", high_rank
    ),
    hoverinfo = "text",
    name = "Negative Significant Genes"
  )
  
  # Add one dummy trace for the self-gene legend (always shown)
  p <- p %>% add_segments(
    x = 0, xend = 0, y = 0, yend = 0.00001,
    line = list(color = "black", dash = "dot", width = 2),
    name = self_gene_label,
    showlegend = TRUE,
    hoverinfo = "none"
  )
  
  # Conditionally add vertical line only if target gene exists
  if (!is.null(df_target)) {
    # p <- p %>% add_segments(
    #   data = df_target,
    #   x = ~high_zscore, xend = ~high_zscore,
    #   y = 0, yend = max(df_dens$y),
    #   line = list(color = df_target$color, dash = "dot", width = 2),
    #   hoverinfo = "text",
    #   text = ~paste("Target Gene:", GeneName,
    #                 "<br>Z-score:", round(high_zscore, 2),
    #                 "<br>Rank:", high_rank),
    #   showlegend = FALSE
    # )
    
    p <- p %>% add_segments(
      data = df_target,
      x = ~high_zscore,
      xend = ~high_zscore,
      y = 0,
      yend = max(df_dens$y),
      
      line = list(
        color = "black",
        dash = "dot",
        width = 4
      ),
      
      hoverinfo = "text",
      
      text = ~paste(
        "<b>Self-Gene:</b> ", GeneName,
        "<br><b>nIG:</b> ", round(high_zscore, 2),
        "<br><b>Rank:</b> ", high_rank
      ),
      
      showlegend = FALSE
    )
    
  }
  
  # Layout
  p <- p %>% layout(
    title = "",
    xaxis = list(title = "normalized IG (nIG)"),
    yaxis = list(title = "Density")
  )
  
  return(p)
}


plotEnrichment <- function(df){
  # Prepare data
  df_top <- df %>%
    slice_min(pvalue, n = 5) %>%
    mutate(GeneRatio_num = parse_ratio(GeneRatio)) %>%
    mutate(Description = factor(Description, levels = Description[order(GeneRatio_num)]))
  
  
  # Plotly dot plot
  p <- plot_ly(data = df_top, x = ~GeneRatio_num, y = ~Description,
               type = 'scatter', mode = 'markers',
               marker = list(
                 sizemode = 'diameter',
                 sizeref = 1,        # tweak this to control dot size scale
                 size = ~Count,
                 color = ~pvalue,
                 colorscale = list(
                   c(0, "#e06663"),
                   c(0.5, "#c04cbc"),
                   c(1, "#327eba")
                 ),
                 reversescale = FALSE,
                 showscale = TRUE,
                 colorbar = list(title = "pvalue", tickformat = ".1e")

               ),
               text = ~paste("Term:", Description,
                             "<br>GeneRatio:", round(GeneRatio_num, 3),
                             "<br>Count:", Count,
                             "<br>p.adjust:", signif(p.adjust, 3)),
               hoverinfo = 'text'
  ) %>%
    layout(
      title = "",
      xaxis = list(title = "GeneRatio"),
      yaxis = list(title = "", categoryorder = "array", categoryarray = levels(df_top$Description))
    )
  return(p)

}

# plotEnrichment <- function(df) {
#   required_cols <- c("pvalue", "GeneRatio", "Description", "Count", "p.adjust")
#   if (!all(required_cols %in% colnames(df))) {
#     stop("Input data frame is missing required columns: ",
#          paste(setdiff(required_cols, colnames(df)), collapse = ", "))
#   }
# 
#   df_top <- df %>%
#     slice_min(pvalue, n = 5) %>%
#     mutate(
#       GeneRatio_num = parse_ratio(GeneRatio),
#       Description = factor(Description, levels = Description[order(GeneRatio_num)])
#     )
# 
#   # Ensure colorbar is a named list
#   colorbar_settings <- list(title = "p.adjust")
# 
#   p <- plot_ly(
#     data = df_top,
#     x = ~GeneRatio_num,
#     y = ~Description,
#     type = 'scatter',
#     mode = 'markers',
#     marker = list(
#       sizemode = 'diameter',
#       sizeref = 2,
#       size = ~Count,
#       color = ~p.adjust,
#       colorscale = 'RdBu',
#       reversescale = TRUE,
#       showscale = TRUE,
#       colorbar = colorbar_settings
#     ),
#     text = ~paste("Term:", Description,
#                   "<br>GeneRatio:", round(GeneRatio_num, 3),
#                   "<br>Count:", Count,
#                   "<br>p.adjust:", signif(p.adjust, 3)),
#     hoverinfo = 'text'
#   ) %>%
#     layout(
#       title = "Top 5 Enriched Terms",
#       xaxis = list(title = "GeneRatio"),
#       yaxis = list(
#         title = "",
#         categoryorder = "array",
#         categoryarray = levels(df_top$Description)
#       )
#     )
#   return (p)
# }

plotViolin <- function(df, celltype, protein){
  df_plot <- df %>%
    mutate(
      group_raw = ifelse(celltype.l1 == celltype, celltype, "Others")
    )
  
  group_counts <- df_plot %>%
    count(group_raw) %>%
    mutate(label = paste0(group_raw, " (n=", n, ")"))
  
  # Create new labeled factor for plotting
  df_plot <- df_plot %>%
    left_join(group_counts, by = c("group_raw")) %>%
    mutate(group = factor(label, levels = group_counts$label))
  
  # Define color map using new labels
  color_map <- setNames(c("steelblue", "gray"), group_counts$label)
  
  # Create plot
  p <- df_plot %>%
    plot_ly(
      x = ~group,
      y = ~.data[[protein]],
      type = 'violin',
      color = ~group,
      colors = color_map,
      box = list(visible = TRUE),
      meanline = list(visible = TRUE)
    ) %>%
    layout(
      title = "",
      yaxis = list(title = "Expression Level"),
      xaxis = list(title = ""),
      showlegend = F
    )
  return(p)
  
}
