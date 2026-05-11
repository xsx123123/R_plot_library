plot_go_dotplot <- function(data, top_n = 10, x_var = "FoldEnrichment", 
                            wrap_length = 50, highlight_pathways = NULL) {
  # =========================================================
  # Function:
  #   Generate a dot plot for GO enrichment results.
  #
  # Features:
  #   1. Select the top_n enriched terms within each ONTOLOGY group
  #   2. Wrap long Description text for better display
  #   3. Highlight selected pathways on the y-axis labels
  #      using red color and bold font
  #   4. Display results in separate facets by ONTOLOGY
  #
  # Arguments:
  #   data               : A data frame containing GO enrichment results.
  #                        It should include columns such as ONTOLOGY,
  #                        Description, p.adjust, Count, and the variable
  #                        specified by x_var.
  #   top_n              : Number of top terms to keep for each ONTOLOGY.
  #   x_var              : Column name used for the x-axis
  #                        (e.g. "FoldEnrichment").
  #   wrap_length        : Maximum width used to wrap Description text.
  #   highlight_pathways : A character vector of pathway names to highlight.
  #
  # Returns:
  #   A ggplot object.
  # =========================================================
  
  # Load required packages
  require(dplyr)
  require(stringr)
  require(ggplot2)
  require(ggtext)
  require(rlang)
  
  # ---------------------------------------------------------
  # Step 1. Data preprocessing
  #   - Convert p.adjust to numeric
  #   - Group data by ONTOLOGY
  #   - Sort terms by adjusted p-value within each group
  #   - Keep the top_n terms for each ONTOLOGY
  #   - Wrap Description text to avoid overly long axis labels
  # ---------------------------------------------------------
  plot_df <- data %>%
    mutate(p.adjust = as.numeric(p.adjust)) %>%
    group_by(ONTOLOGY) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    mutate(Description = str_wrap(Description, width = wrap_length))
  
  # ---------------------------------------------------------
  # Step 2. Process highlighted pathways
  #   - Since Description has already been wrapped,
  #     highlight_pathways must also be wrapped in the same way
  #     to ensure correct matching
  #   - Matching pathway labels will be converted into
  #     HTML-styled text for rendering by ggtext
  # ---------------------------------------------------------
  if (!is.null(highlight_pathways)) {
    highlight_wrapped <- str_wrap(highlight_pathways, width = wrap_length)
    
    plot_df <- plot_df %>%
      mutate(
        Description_label = ifelse(
          Description %in% highlight_wrapped,
          paste0("<span style='color:red;'><b>", Description, "</b></span>"),
          Description
        )
      )
  } else {
    # If no highlighted pathways are provided,
    # use the original wrapped Description as label
    plot_df <- plot_df %>%
      mutate(Description_label = Description)
  }
  
  # ---------------------------------------------------------
  # Step 3. Set the order of y-axis labels
  #   - Sort by ONTOLOGY and the selected x-axis variable
  #   - Convert Description_label to a factor so that
  #     the plotting order can be controlled
  # ---------------------------------------------------------
  plot_df <- plot_df %>%
    arrange(ONTOLOGY, !!sym(x_var)) %>%
    mutate(Description_label = factor(Description_label, levels = unique(Description_label)))
  
  # ---------------------------------------------------------
  # Step 4. Create the dot plot
  #   - x-axis  : selected variable (x_var)
  #   - y-axis  : pathway descriptions
  #   - size    : gene count
  #   - color   : adjusted p-value
  #   - facets  : separated by ONTOLOGY
  # ---------------------------------------------------------
  p <- ggplot(plot_df, aes(x = !!sym(x_var), y = Description_label)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    
    # Define color gradient for adjusted p-values
    scale_color_gradient(low = "red", high = "blue") +
    
    # Create one facet for each ONTOLOGY
    facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
    
    # Apply clean publication-style theme
    theme_pubclean() +
    theme(
      legend.position = "right",
      
      # Use markdown rendering so highlighted labels can display
      # custom HTML styles such as color and bold font
      axis.text.y = ggtext::element_markdown(size = 10),
      axis.text.x = element_text(size = 10),
      
      # Style facet labels
      strip.text = element_text(size = 12, face = "bold"),
      strip.background = element_rect(fill = "grey90"),
      
      # Add horizontal guide lines for readability
      panel.grid.major.y = element_line(linetype = "dashed", color = "grey80"),
      
      # Center and bold the plot title
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    ) +
    
    # Set labels and legend titles
    labs(
      title = "GO Enrichment Analysis",
      x = x_var,
      y = "",
      color = "p.adjust",
      size = "Gene Count"
    )
  
  # Return the ggplot object
  return(p)
}
