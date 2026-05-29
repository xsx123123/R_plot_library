#' Plot GO Enrichment Dot Plot
#'
#' Visualize GO enrichment analysis results as a dot plot.
#' The size of dots represents gene counts, and color represents adjusted p-value.
#' Results are grouped and faceted by GO ONTOLOGY (BP/CC/MF).
#'
#' @param data A data frame containing GO enrichment results. Required columns:
#'   Description, p.adjust, Count, ONTOLOGY, and the column specified by x_var.
#' @param top_n Integer. Number of top enriched terms to show per ONTOLOGY group. Default: 10.
#' @param x_var Character. Column name for x-axis values. Default: "FoldEnrichment".
#' @param wrap_length Integer. Maximum width for wrapping Description text labels. Default: 50.
#' @param highlight_pathways Character vector. Specific pathway names to highlight in red and bold.
#'   Default: NULL (no highlighting).
#' @param facet Logical. Whether to facet the plot by ONTOLOGY. Default: TRUE.
#'
#' @return A ggplot object.
#'
#' @examples
#' plot_go_dotplot(go_result, top_n = 15, highlight_pathways = c("cell cycle", "apoptosis"))
#' plot_go_dotplot(go_result, facet = FALSE, x_var = "GeneRatio")
plot_go_dotplot <- function(data, top_n = 10, x_var = "FoldEnrichment", 
                            wrap_length = 50, highlight_pathways = NULL,
                            facet = TRUE) {
  require(dplyr)
  require(stringr)
  require(ggplot2)
  require(ggtext)
  require(rlang)
  
  plot_df <- data %>%
    mutate(p.adjust = as.numeric(p.adjust)) %>%
    group_by(ONTOLOGY) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    mutate(
      Description = str_wrap(Description, width = wrap_length),
      Description = gsub("\n", "<br>", Description)
    )
  
  if (!is.null(highlight_pathways)) {
    highlight_wrapped <- str_wrap(highlight_pathways, width = wrap_length)
    highlight_wrapped <- gsub("\n", "<br>", highlight_wrapped)
    
    plot_df <- plot_df %>%
      mutate(
        Description_label = ifelse(
          Description %in% highlight_wrapped,
          paste0("<span style='color:red;'><b>", Description, "</b></span>"),
          Description
        )
      )
  } else {
    plot_df <- plot_df %>%
      mutate(Description_label = Description)
  }
  
  if (facet) {
    plot_df <- plot_df %>%
      arrange(ONTOLOGY, !!sym(x_var)) %>%
      mutate(Description_label = factor(Description_label, levels = unique(Description_label)))
  } else {
    plot_df <- plot_df %>%
      arrange(!!sym(x_var)) %>%
      mutate(Description_label = factor(Description_label, levels = unique(Description_label)))
  }
  
  p <- ggplot(plot_df, aes(x = !!sym(x_var), y = Description_label)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "red", high = "blue") +
    theme_pubclean() +
    theme(
      legend.position = "right",
      axis.text.y = ggtext::element_markdown(size = 10),
      axis.text.x = element_text(size = 10),
      strip.text = element_text(size = 12, face = "bold"),
      strip.background = element_rect(fill = "grey90"),
      panel.grid.major.y = element_line(linetype = "dashed", color = "grey80"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    ) +
    labs(
      title = "GO Enrichment Analysis",
      x = x_var,
      y = "",
      color = "p.adjust",
      size = "Gene Count"
    )
  
  if (facet) {
    p <- p + facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y")
  }
  
  return(p)
}

#' Plot KEGG Enrichment Dot Plot
#'
#' Visualize KEGG pathway enrichment analysis results as a dot plot.
#' The size of dots represents gene counts, and color represents adjusted p-value.
#'
#' @param data A data frame containing KEGG enrichment results. Required columns:
#'   Description, p.adjust, Count, and the column specified by x_var.
#' @param top_n Integer. Number of top enriched pathways to show. Default: 10.
#' @param x_var Character. Column name for x-axis values. Default: "FoldEnrichment".
#' @param wrap_length Integer. Maximum width for wrapping Description text labels. Default: 50.
#' @param highlight_pathways Character vector. Specific pathway names to highlight in red and bold.
#'   Default: NULL (no highlighting).
#'
#' @return A ggplot object.
#'
#' @examples
#' plot_kegg_dotplot(kegg_result, top_n = 15, highlight_pathways = c("MAPK signaling pathway"))
#' plot_kegg_dotplot(kegg_result, x_var = "GeneRatio")
plot_kegg_dotplot <- function(data, top_n = 10, x_var = "FoldEnrichment", 
                              wrap_length = 50, highlight_pathways = NULL) {
  require(dplyr)
  require(stringr)
  require(ggplot2)
  require(ggtext)
  require(rlang)
  
  plot_df <- data %>%
    mutate(p.adjust = as.numeric(p.adjust)) %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n) %>%
    mutate(
      Description = str_wrap(Description, width = wrap_length),
      Description = gsub("\n", "<br>", Description)
    )
  
  if (!is.null(highlight_pathways)) {
    highlight_wrapped <- str_wrap(highlight_pathways, width = wrap_length)
    highlight_wrapped <- gsub("\n", "<br>", highlight_wrapped)
    
    plot_df <- plot_df %>%
      mutate(
        Description_label = ifelse(
          Description %in% highlight_wrapped,
          paste0("<span style='color:red;'><b>", Description, "</b></span>"),
          Description
        )
      )
  } else {
    plot_df <- plot_df %>%
      mutate(Description_label = Description)
  }
  
  plot_df <- plot_df %>%
    arrange(!!sym(x_var)) %>%
    mutate(Description_label = factor(Description_label, levels = unique(Description_label)))
  
  p <- ggplot(plot_df, aes(x = !!sym(x_var), y = Description_label)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "red", high = "blue") +
    theme_pubclean() +
    theme(
      legend.position = "right",
      axis.text.y = ggtext::element_markdown(size = 10),
      axis.text.x = element_text(size = 10),
      strip.text = element_text(size = 12, face = "bold"),
      strip.background = element_rect(fill = "grey90"),
      panel.grid.major.y = element_line(linetype = "dashed", color = "grey80"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    ) +
    labs(
      title = "KEGG Enrichment Analysis",
      x = x_var,
      y = "",
      color = "p.adjust",
      size = "Gene Count"
    )
  
  return(p)
}