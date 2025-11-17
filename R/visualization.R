#' @importFrom magrittr %>%
#' @importFrom dplyr mutate filter group_by arrange
#' @importFrom karyoploteR plotKaryotype kpAddCytobands kpAddChromosomeNames kpPlotRegions kpPoints kpPlotMarkers
#' @importFrom regioneR toGRanges
#' @importFrom GenomicRanges GRanges
#'
#' Plot Chromosome Ideogram
#' @importFrom magrittr %>%
#' Plot Chromosome Ideogram
#'
#' Creates chromosome ideogram plots showing gene locations on chromosomes.
#' Generates both regional and point-based visualizations.
#'
#' @param candidate_data Data frame of candidate genes from screen_candidate_genes()
#' @param output_dir Character string specifying output directory (default: ".")
#' @param plot_type Character string specifying plot type:
#'   "regions", "points", or "both" (default: "both")
#' @param genome Character string specifying genome version (default: "hg38")
#' @param width Numeric width of PDF plots in inches (default: 7)
#' @param height Numeric height of PDF plots in inches (default: 9)
#'
#' @return Invisibly returns a list with plot objects and chromosome statistics
#'
#' @details
#' The function generates high-quality chromosome ideogram plots showing:
#' 1. Regional visualization: Genes shown as red regions on chromosomes
#' 2. Point visualization: Genes shown as red points at gene midpoints
#' 3. Chromosome distribution statistics
#'
#' @examples
#' \dontrun{
#' # Plot both regional and point ideograms
#' candidates <- screen_candidate_genes(processed_data)
#' plot_ideogram(candidates, output_dir = "plots")
#'
#' # Plot only regional ideogram
#' plot_ideogram(candidates, plot_type = "regions")
#' }
#'
#' @export
plot_ideogram <- function(candidate_data, output_dir = ".", plot_type = "both",
                         genome = "hg38", width = 7, height = 9) {

  # Input validation
  if (!is.data.frame(candidate_data)) {
    stop("Input must be a data frame")
  }

  required_cols <- c("chromosome", "start_position_on_the_genomic_accession",
                    "end_position_on_the_genomic_accession")
  missing_cols <- setdiff(required_cols, colnames(candidate_data))

  if (length(missing_cols) > 0) {
    stop("Missing required columns for plotting: ", paste(missing_cols, collapse = ", "))
  }

  valid_types <- c("regions", "points", "both")
  if (!plot_type %in% valid_types) {
    stop("plot_type must be one of: ", paste(valid_types, collapse = ", "))
  }

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Filter genes with valid position information
  plot_data <- candidate_data[!is.na(candidate_data$start_position_on_the_genomic_accession) &
                             !is.na(candidate_data$end_position_on_the_genomic_accession), ]

  if (nrow(plot_data) == 0) {
    warning("No genes with valid position information found for plotting")
    return(invisible(NULL))
  }

  # Build clean data frame for plotting
  gene_data <- data.frame(
    Chr = as.character(plot_data$chromosome),
    Start = as.numeric(plot_data$start_position_on_the_genomic_accession),
    End = as.numeric(plot_data$end_position_on_the_genomic_accession),
    Gene = plot_data$Gene,
    stringsAsFactors = FALSE
  )

  # Chromosome distribution statistics
  chr_gene_count <- as.data.frame(table(gene_data$Chr))
  colnames(chr_gene_count) <- c("Chr", "Count")

  cat("========== Chromosome Distribution Statistics ==========\n")
  print(chr_gene_count)
  cat("Total genes with position info:", nrow(gene_data), "\n")
  cat("=====================================================\n")

  # Build GRanges object
  gene_gr <- GRanges(
    seqnames = paste0("chr", gene_data$Chr),
    ranges = IRanges(start = gene_data$Start, end = gene_data$End),
    gene = gene_data$Gene
  )

  # Calculate midpoint positions for point plot
  mid_pos <- (start(gene_gr) + end(gene_gr)) / 2

  # Chromosome order: Y, X, 22, ..., 1
  chr.order <- rev(paste0("chr", c(1:22, "X", "Y")))

  # Initialize results list
  plot_results <- list()

  # Plot 1: Regional ideogram
  if (plot_type %in% c("regions", "both")) {
    regions_pdf <- file.path(output_dir, "chromosome_ideogram_regions.pdf")
    regions_png <- file.path(output_dir, "chromosome_ideogram_regions.png")

    pdf(regions_pdf, width = width, height = height)
    kp1 <- plotKaryotype(
      genome = genome,
      chromosomes = chr.order,
      plot.type = 1,
      ideogram.plotter = kpAddCytobands,
      labels.plotter = kpAddChromosomeNames
    )

    kpPlotRegions(
      kp1,
      data = gene_gr,
      r0 = 0.0,
      r1 = 0.5,
      col = "red",
      border = NA
    )
    dev.off()

    plot_results$regions_plot <- kp1
    plot_results$regions_files <- c(pdf = regions_pdf, png = regions_png)
  }

  # Plot 2: Point ideogram
  if (plot_type %in% c("points", "both")) {
    points_pdf <- file.path(output_dir, "chromosome_ideogram_points.pdf")
    points_png <- file.path(output_dir, "chromosome_ideogram_points.png")

    pdf(points_pdf, width = width, height = height)
    kp2 <- plotKaryotype(
      genome = genome,
      chromosomes = chr.order,
      plot.type = 1,
      ideogram.plotter = kpAddCytobands,
      labels.plotter = kpAddChromosomeNames
    )

    kpPoints(
      kp2,
      chr = seqnames(gene_gr),
      x = mid_pos,
      y = 0.25,
      pch = 16,
      cex = 0.3,
      col = "red"
    )
    dev.off()

    plot_results$points_plot <- kp2
    plot_results$points_files <- c(pdf = points_pdf, png = points_png)
  }

  plot_results$chromosome_stats <- chr_gene_count
  plot_results$genes_plotted <- nrow(gene_data)

  cat("\nGenerated plot files:\n")
  if (plot_type %in% c("regions", "both")) {
    cat("- chromosome_ideogram_regions.pdf\n")
  }
  if (plot_type %in% c("points", "both")) {
    cat("- chromosome_ideogram_points.pdf\n")
  }

  return(invisible(plot_results))
}

#' Plot Top Cytobands Line Plot
#'
#' Creates a line and scatter plot showing attribution values for top N cytobands.
#'
#' @param candidate_data Data frame of candidate genes
#' @param top_n Numeric number of top cytobands to show (default: 20)
#' @param output_dir Character string specifying output directory (default: ".")
#' @param width Numeric plot width in inches (default: 8)
#' @param height Numeric plot height in inches (default: 5)
#'
#' @return A ggplot object
#'
#' @details
#' The function creates a line plot with the following features:
#' 1. Shows absolute attribution sums for both conditions
#' 2. Displays top N cytobands sorted by Attribution.M_sum
#' 3. Includes both lines and points for better visualization
#' 4. Automatic saving to PDF file
#'
#' @examples
#' \dontrun{
#' # Create plot for top 20 cytobands
#' plot_obj <- plot_top_cytobands(candidates, top_n = 20)
#'
#' # Save custom plot
#' ggsave("custom_plot.pdf", plot_obj)
#' }
#'
#' @export
plot_top_cytobands <- function(candidate_data, top_n = 20, output_dir = ".",
                              width = 8, height = 5) {

  # Input validation
  if (!is.data.frame(candidate_data)) {
    stop("Input must be a data frame")
  }

  required_cols <- c("Cytoband", "Attribution.G", "Attribution.M")
  missing_cols <- setdiff(required_cols, colnames(candidate_data))

  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!is.numeric(top_n) || top_n <= 0) {
    stop("top_n must be a positive number")
  }

  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Add absolute value columns
  plot_data <- candidate_data %>%
    mutate(
      Attribution.G_abs = abs(Attribution.G),
      Attribution.M_abs = abs(Attribution.M)
    )

  # Aggregate by cytoband
  cytoband_summary <- plot_data %>%
    group_by(Cytoband) %>%
    summarize(
      Attribution.G_sum = sum(Attribution.G_abs, na.rm = TRUE),
      Attribution.M_sum = sum(Attribution.M_abs, na.rm = TRUE),
      n_genes = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(Attribution.M_sum))

  # Take top N cytobands
  top_cytobands <- cytoband_summary %>%
    slice_head(n = top_n)

  # Create the plot
  p <- ggplot(
    top_cytobands %>%
      mutate(Cytoband = fct_reorder(Cytoband, Attribution.M_sum, .desc = TRUE)),
    aes(x = Cytoband)
  ) +
    geom_line(aes(y = Attribution.G_sum, group = 1, color = "Gene-only (G)"), linewidth = 1) +
    geom_point(aes(y = Attribution.G_sum, color = "Gene-only (G)"), size = 3) +
    geom_line(aes(y = Attribution.M_sum, group = 1, color = "Multimodal (M)"), linewidth = 1) +
    geom_point(aes(y = Attribution.M_sum, color = "Multimodal (M)"), size = 3) +
    scale_color_manual(
    values = c("Gene-only (G)" = "blue", "Multimodal (M)" = "red"),
    name = "Model Type"
  ) +
  labs(
    title = paste("Top", top_n, "Cytobands by Multimodal Attribution"),
    subtitle = "Ranked by total multimodal attribution sum",
    x = "Cytoband",
    y = "Sum of Absolute Attribution Value",
    color = "Model Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.subtitle = element_text(size = 10, color = "gray40")
  )

  # Save the plot
  output_file <- file.path(output_dir, paste0("Top", top_n, "_Cytoband_lineplot.pdf"))
  ggsave(output_file, plot = p, width = width, height = height, device = cairo_pdf)

  cat("Top cytobands plot saved to:", output_file, "\n")

  return(p)
}

#' Plot Top Cytobands by Difference
#'
#' Creates a line and scatter plot showing attribution values for top N cytobands
#' sorted by the difference (Attribution.M_sum - Attribution.G_sum).
#'
#' @param candidate_data Data frame of candidate genes
#' @param top_n Numeric number of top cytobands to show (default: 20)
#' @param output_dir Character string specifying output directory (default: ".")
#' @param width Numeric plot width in inches (default: 8)
#' @param height Numeric plot height in inches (default: 5)
#'
#' @return A ggplot object
#'
#' @details
#' This function creates a line plot showing cytobands ranked by the difference
#' between multimodal and gene-only attribution values. Positive differences
#' indicate cytobands where the multimodal model shows greater importance.
#'
#' @examples
#' \dontrun{
#' # Create plot for top 20 cytobands by difference
#' plot_obj <- plot_top_cytobands_by_difference(candidates, top_n = 20)
#' }
#'
#' @export
plot_top_cytobands_by_difference <- function(candidate_data, top_n = 20, output_dir = ".",
                                            width = 8, height = 5) {

  # Input validation
  if (!is.data.frame(candidate_data)) {
    stop("Input must be a data frame")
  }

  required_cols <- c("Cytoband", "Attribution.G", "Attribution.M")
  missing_cols <- setdiff(required_cols, colnames(candidate_data))

  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!is.numeric(top_n) || top_n <= 0) {
    stop("top_n must be a positive number")
  }

  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Add absolute value columns
  plot_data <- candidate_data %>%
    mutate(
      Attribution.G_abs = abs(Attribution.G),
      Attribution.M_abs = abs(Attribution.M)
    )

  # Aggregate by cytoband and calculate difference
  cytoband_summary <- plot_data %>%
    group_by(Cytoband) %>%
    summarize(
      Attribution.G_sum = sum(Attribution.G_abs, na.rm = TRUE),
      Attribution.M_sum = sum(Attribution.M_abs, na.rm = TRUE),
      Difference = Attribution.M_sum - Attribution.G_sum,
      n_genes = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(Difference))  # Sort by difference

  # Take top N cytobands by difference
  top_cytobands <- cytoband_summary %>%
    slice_head(n = top_n)

  # Create the plot
  p <- ggplot(
    top_cytobands %>%
      mutate(Cytoband = fct_reorder(Cytoband, Difference, .desc = TRUE)),
    aes(x = Cytoband)
  ) +
    geom_line(aes(y = Attribution.G_sum, group = 1, color = "Gene-only (G)"), linewidth = 1) +
    geom_point(aes(y = Attribution.G_sum, color = "Gene-only (G)"), size = 3) +
    geom_line(aes(y = Attribution.M_sum, group = 1, color = "Multimodal (M)"), linewidth = 1) +
    geom_point(aes(y = Attribution.M_sum, color = "Multimodal (M)"), size = 3) +
    scale_color_manual(
      values = c("Gene-only (G)" = "blue", "Multimodal (M)" = "red"),
      name = "Model Type"
    ) +
    labs(
      title = paste("Top", top_n, "Cytobands by Difference (M - G)"),
      subtitle = "Cytobands showing greatest multimodal advantage",
      x = "Cytoband",
      y = "Sum of Absolute Attribution Value",
      color = "Model Type"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      plot.subtitle = element_text(size = 10, color = "gray40")
    )

  # Save the plot
  output_file <- file.path(output_dir, paste0("Top", top_n, "_Cytoband_by_Difference.pdf"))
  ggsave(output_file, plot = p, width = width, height = height, device = cairo_pdf)

  cat("Top cytobands by difference plot saved to:", output_file, "\n")

  return(p)
}

#' Plot Top Cytobands by Ratio
#'
#' Creates a line and scatter plot showing attribution values for top N cytobands
#' sorted by the ratio (Attribution.M_sum / Attribution.G_sum).
#'
#' @param candidate_data Data frame of candidate genes
#' @param top_n Numeric number of top cytobands to show (default: 20)
#' @param output_dir Character string specifying output directory (default: ".")
#' @param width Numeric plot width in inches (default: 8)
#' @param height Numeric plot height in inches (default: 5)
#'
#' @return A ggplot object
#'
#' @details
#' This function creates a line plot showing cytobands ranked by the ratio
#' of multimodal to gene-only attribution values. Ratios > 1 indicate
#' cytobands where the multimodal model shows proportionally greater importance.
#'
#' @examples
#' \dontrun{
#' # Create plot for top 20 cytobands by ratio
#' plot_obj <- plot_top_cytobands_by_ratio(candidates, top_n = 20)
#' }
#'
#' @export
plot_top_cytobands_by_ratio <- function(candidate_data, top_n = 20, output_dir = ".",
                                       width = 8, height = 5) {

  # Input validation
  if (!is.data.frame(candidate_data)) {
    stop("Input must be a data frame")
  }

  required_cols <- c("Cytoband", "Attribution.G", "Attribution.M")
  missing_cols <- setdiff(required_cols, colnames(candidate_data))

  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!is.numeric(top_n) || top_n <= 0) {
    stop("top_n must be a positive number")
  }

  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Add absolute value columns
  plot_data <- candidate_data %>%
    mutate(
      Attribution.G_abs = abs(Attribution.G),
      Attribution.M_abs = abs(Attribution.M)
    )

  # Aggregate by cytoband and calculate ratio
  cytoband_summary <- plot_data %>%
    group_by(Cytoband) %>%
    summarize(
      Attribution.G_sum = sum(Attribution.G_abs, na.rm = TRUE),
      Attribution.M_sum = sum(Attribution.M_abs, na.rm = TRUE),
      Ratio = ifelse(Attribution.G_sum == 0,
                     NA_real_,
                     Attribution.M_sum / Attribution.G_sum),
      n_genes = n(),
      .groups = "drop"
    ) %>%
    filter(!is.na(Ratio)) %>%  # Remove NA ratios
    arrange(desc(Ratio))  # Sort by ratio

  # Take top N cytobands by ratio
  top_cytobands <- cytoband_summary %>%
    slice_head(n = top_n)

  # Create the plot
  p <- ggplot(
    top_cytobands %>%
      mutate(Cytoband = fct_reorder(Cytoband, Ratio, .desc = TRUE)),
    aes(x = Cytoband)
  ) +
    geom_line(aes(y = Attribution.G_sum, group = 1, color = "Gene-only (G)"), linewidth = 1) +
    geom_point(aes(y = Attribution.G_sum, color = "Gene-only (G)"), size = 3) +
    geom_line(aes(y = Attribution.M_sum, group = 1, color = "Multimodal (M)"), linewidth = 1) +
    geom_point(aes(y = Attribution.M_sum, color = "Multimodal (M)"), size = 3) +
    scale_color_manual(
      values = c("Gene-only (G)" = "blue", "Multimodal (M)" = "red"),
      name = "Model Type"
    ) +
    labs(
      title = paste("Top", top_n, "Cytobands by Ratio (M / G)"),
      subtitle = "Cytobands showing highest multimodal-to-gene ratio",
      x = "Cytoband",
      y = "Sum of Absolute Attribution Value",
      color = "Model Type"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      plot.subtitle = element_text(size = 10, color = "gray40")
    )

  # Save the plot
  output_file <- file.path(output_dir, paste0("Top", top_n, "_Cytoband_by_Ratio.pdf"))
  ggsave(output_file, plot = p, width = width, height = height, device = cairo_pdf)

  cat("Top cytobands by ratio plot saved to:", output_file, "\n")

  return(p)
}
