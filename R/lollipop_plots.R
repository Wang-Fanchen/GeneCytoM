#' @importFrom magrittr %>%
#' @importFrom dplyr slice_max mutate arrange
#'
#' Plot Lollipop Charts for Cytoband Analysis
#' @importFrom magrittr %>%
#' Plot Lollipop Charts for Cytoband Analysis
#'
#' Creates lollipop charts showing difference and ratio between attribution values
#' for top cytobands.
#'
#' @param cytoband_summary Data frame with cytoband summary statistics
#' @param top_n Numeric number of top cytobands to show (default: 10)
#' @param output_dir Character string specifying output directory (default: ".")
#' @param plot_type Character string specifying which plots to create:
#'   "difference", "ratio", or "both" (default: "both")
#' @param width Numeric plot width in inches (default: 7)
#' @param height Numeric plot height in inches (default: 5)
#'
#' @return A list containing ggplot objects
#'
#' @details
#' The function creates two types of lollipop charts:
#' 1. Difference plot: Attribution.M_sum - Attribution.G_sum
#' 2. Ratio plot: Attribution.M_sum / Attribution.G_sum
#'
#' Both plots show the top N cytobands sorted by Attribution.M_sum.
#'
#' @examples
#' \dontrun{
#' # Create both difference and ratio plots
#' summary_data <- analyze_cytoband_location(candidates)
#' plots <- plot_lollipop_charts(summary_data, top_n = 10)
#'
#' # Create only difference plot
#' diff_plot <- plot_lollipop_charts(summary_data, plot_type = "difference")
#' }
#'
#' @export
plot_lollipop_charts <- function(cytoband_summary, top_n = 10, output_dir = ".",
                                plot_type = "both", width = 7, height = 5) {

  # Input validation
  if (!is.data.frame(cytoband_summary)) {
    stop("Input must be a data frame")
  }

  required_cols <- c("Cytoband", "Attribution.G_sum", "Attribution.M_sum")
  missing_cols <- setdiff(required_cols, colnames(cytoband_summary))

  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  valid_types <- c("difference", "ratio", "both")
  if (!plot_type %in% valid_types) {
    stop("plot_type must be one of: ", paste(valid_types, collapse = ", "))
  }

  if (!is.numeric(top_n) || top_n <= 0) {
    stop("top_n must be a positive number")
  }

  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Calculate difference and ratio
  plot_data <- cytoband_summary %>%
    mutate(
      Difference = Attribution.M_sum - Attribution.G_sum,
      Ratio = ifelse(Attribution.G_sum == 0,
                    NA_real_,
                    Attribution.M_sum / Attribution.G_sum)
    ) %>%
    arrange(desc(Attribution.M_sum))

  # Take top N cytobands
  top_data <- plot_data %>%
    slice_head(n = top_n)

  # Initialize results list
  plot_results <- list()

  # Plot 1: Difference lollipop chart
  if (plot_type %in% c("difference", "both")) {
    diff_plot <- ggplot(
      top_data,
      aes(x = reorder(Cytoband, Difference), y = Difference)
    ) +
      geom_segment(aes(xend = Cytoband, y = 0, yend = Difference), color = "gray") +
      geom_point(aes(color = Difference > 0), size = 4) +
      scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red"),
                        name = "Positive") +
      coord_flip() +
      labs(
        title = paste("Difference Between Attribution.M and Attribution.G\n(Top", top_n, "Cytobands by Attribution.M_sum)"),
        x = "Cytoband",
        y = "Difference (Attribution.M_sum - Attribution.G_sum)"
      ) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none")

    # Save difference plot
    diff_file <- file.path(output_dir, "Cytoband_Difference_lollipop.pdf")
    ggsave(diff_file, plot = diff_plot, width = width, height = height, device = cairo_pdf)

    plot_results$difference_plot <- diff_plot
    plot_results$difference_file <- diff_file

    cat("Difference lollipop plot saved to:", diff_file, "\n")
  }

  # Plot 2: Ratio lollipop chart
  if (plot_type %in% c("ratio", "both")) {
    # Remove NA ratios for plotting
    ratio_data <- top_data[!is.na(top_data$Ratio), ]

    if (nrow(ratio_data) > 0) {
      ratio_plot <- ggplot(
        ratio_data,
        aes(x = reorder(Cytoband, Ratio), y = Ratio)
      ) +
        geom_segment(aes(xend = Cytoband, y = 0, yend = Ratio), color = "gray") +
        geom_point(aes(color = Ratio > 1), size = 4) +
        scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red"),
                          name = "Ratio > 1") +
        coord_flip() +
        labs(
          title = paste("Ratio Between Attribution.M and Attribution.G\n(Top", top_n, "Cytobands by Attribution.M_sum)"),
          x = "Cytoband",
          y = "Ratio (Attribution.M_sum / Attribution.G_sum)"
        ) +
        theme_minimal(base_size = 12) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
        theme(legend.position = "none")

      # Save ratio plot
      ratio_file <- file.path(output_dir, "Cytoband_Ratio_lollipop.pdf")
      ggsave(ratio_file, plot = ratio_plot, width = width, height = height, device = cairo_pdf)

      plot_results$ratio_plot <- ratio_plot
      plot_results$ratio_file <- ratio_file

      cat("Ratio lollipop plot saved to:", ratio_file, "\n")
    } else {
      warning("No valid ratio data available for plotting")
      plot_results$ratio_plot <- NULL
    }
  }

  return(plot_results)
}

#' Create Publication-Ready Plot Theme
#'
#' Internal function to apply consistent theming across plots.
#'
#' @param base_size Numeric base font size
#' @param font_family Character string for font family
#'
#' @return theme object for ggplot2
#'
#' @keywords internal
create_plot_theme <- function(base_size = 12, font_family = "sans") {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(family = font_family),
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
      axis.title = element_text(size = base_size + 1, face = "bold"),
      axis.text = element_text(size = base_size),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95")
    )
}
