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