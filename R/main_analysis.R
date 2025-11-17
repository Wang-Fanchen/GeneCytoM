#' Run Complete Cytoband Analysis Pipeline
#'
#' Executes the complete cytoband analysis pipeline from data processing to visualization.
#' This is the main function that integrates all analysis steps.
#'
#' @param input_file Character string specifying path to input CSV file
#' @param output_dir Character string specifying output directory (default: "results")
#' @param comparison_filter Character string for candidate gene filtering:
#'   "up", "down", "mid", or "all" (default: "up")
#' @param pvalue_threshold_g Numeric p-value threshold for P.Value.G filtering (default: NULL).
#'   If NULL or if P.Value.G column doesn't exist, this filter is skipped.
#' @param pvalue_threshold_m Numeric p-value threshold for P.Value.M filtering (default: NULL).
#'   If NULL or if P.Value.M column doesn't exist, this filter is skipped.
#' @param top_n Numeric number of top cytobands for visualization (default: 20)
#' @param create_plots Logical whether to create visualization plots (default: TRUE)
#' @param save_intermediate Logical whether to save intermediate results (default: TRUE)
#'
#' @return A list containing analysis results:
#'   - processed_data: Full processed gene data with location information
#'   - candidate_genes: Filtered candidate genes
#'   - cytoband_summary: Cytoband-level attribution summary
#'   - comprehensive_summary: 5-column comprehensive summary
#'   - plots: List of plot objects (if create_plots = TRUE)
#'
#' @details
#' This function executes the complete analysis pipeline:
#' 1. Data processing and location matching
#' 2. Candidate gene screening
#' 3. Cytoband attribution analysis
#' 4. Comprehensive summary generation
#' 5. Visualization (ideograms, line plots, lollipop charts)
#' 6. Automatic file saving
#'
#' The input CSV file must contain these required columns:
#' Gene, Local.Index.G, Local.Index.M, Attribution.G, Attribution.M
#' Optional columns: P.Value.G, P.Value.M
#'
#' @examples
#' \dontrun{
#' # Run complete analysis with default parameters
#' results <- run_cytoband_analysis(
#'   input_file = "gene_data.csv",
#'   output_dir = "analysis_results"
#' )
#'
#' # Run analysis with custom parameters
#' results <- run_cytoband_analysis(
#'   input_file = "my_data.csv",
#'   output_dir = "custom_results",
#'   comparison_filter = "all",
#'   pvalue_threshold_g = 0.05,
#'   pvalue_threshold_m = 0.01,
#'   top_n = 30
#' )
#'
#' # Run without P-value filtering
#' results <- run_cytoband_analysis(
#'   input_file = "my_data.csv",
#'   output_dir = "custom_results",
#'   comparison_filter = "up"
#' )
#'
#' # Access specific results
#' candidate_count <- nrow(results$candidate_genes)
#' top_cytobands <- head(results$cytoband_summary, 10)
#' }
#'
#' @export
run_cytoband_analysis <- function(input_file,
                                  output_dir = "results",
                                  comparison_filter = "up",
                                  pvalue_threshold_g = NULL,
                                  pvalue_threshold_m = NULL,
                                  top_n = 20,
                                  create_plots = TRUE,
                                  save_intermediate = TRUE) {

  # Input validation
  if (!is.character(input_file) || length(input_file) != 1) {
    stop("input_file must be a single character string")
  }

  if (!file.exists(input_file)) {
    stop("Input file does not exist: ", input_file)
  }

  if (!is.character(output_dir) || length(output_dir) != 1) {
    stop("output_dir must be a single character string")
  }

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  }

  cat("========== Starting Cytoband Analysis Pipeline ==========\n")
  cat("Input file:", input_file, "\n")
  cat("Output directory:", output_dir, "\n")
  cat("Parameters:\n")
  cat("  - comparison_filter =", comparison_filter, "\n")

  if (!is.null(pvalue_threshold_g)) {
    cat("  - pvalue_threshold_g =", pvalue_threshold_g, "\n")
  } else {
    cat("  - pvalue_threshold_g = NULL (no P.Value.G filtering)\n")
  }

  if (!is.null(pvalue_threshold_m)) {
    cat("  - pvalue_threshold_m =", pvalue_threshold_m, "\n")
  } else {
    cat("  - pvalue_threshold_m = NULL (no P.Value.M filtering)\n")
  }

  cat("  - top_n =", top_n, "\n")
  cat("========================================================\n\n")

  # Initialize results list
  results <- list()

  # Step 1: Read and process input data
  cat("Step 1: Processing input data...\n")
  gene_data <- read.csv(input_file, header = TRUE, check.names = FALSE)
  cat("Loaded", nrow(gene_data), "genes from input file\n")

  processed_data <- process_gene_data(gene_data)

  if (save_intermediate) {
    processed_file <- file.path(output_dir, "01_processed_data.csv")
    write.csv(processed_data, processed_file, row.names = FALSE)
    cat("Processed data saved to:", processed_file, "\n")
  }

  results$processed_data <- processed_data
  cat("✓ Data processing completed\n\n")

  # Step 2: Screen candidate genes
  cat("Step 2: Screening candidate genes...\n")
  candidate_genes <- screen_candidate_genes(
    processed_data,
    comparison_filter = comparison_filter,
    pvalue_threshold_g = pvalue_threshold_g,
    pvalue_threshold_m = pvalue_threshold_m
  )

  if (save_intermediate) {
    candidate_file <- file.path(output_dir, "02_candidate_genes.csv")
    write.csv(candidate_genes, candidate_file, row.names = FALSE)
    cat("Candidate genes saved to:", candidate_file, "\n")
  }

  results$candidate_genes <- candidate_genes
  cat("✓ Candidate screening completed\n\n")

  # Step 3: Cytoband attribution analysis
  cat("Step 3: Performing cytoband attribution analysis...\n")
  cytoband_summary <- analyze_cytoband_location(candidate_genes, output_dir)
  results$cytoband_summary <- cytoband_summary
  cat("✓ Cytoband analysis completed\n\n")

  # Step 4: Generate comprehensive summary
  cat("Step 4: Generating comprehensive cytoband summary...\n")
  comprehensive_summary <- generate_cytoband_summary(processed_data, output_dir)
  results$comprehensive_summary <- comprehensive_summary
  cat("✓ Comprehensive summary completed\n\n")

  # Step 5: Create visualizations
  if (create_plots) {
    cat("Step 5: Creating visualizations...\n")

    plots_dir <- file.path(output_dir, "plots")
    if (!dir.exists(plots_dir)) {
      dir.create(plots_dir)
    }

    plot_results <- list()

    # Ideogram plots
    if (nrow(candidate_genes) > 0) {
      cat("  - Creating chromosome ideograms...\n")
      ideogram_results <- plot_ideogram(candidate_genes, plots_dir)
      plot_results$ideogram <- ideogram_results
    }

    # Top cytobands line plots (4 different sortings)
    if (nrow(cytoband_summary) > 0) {
      cat("  - Creating top cytobands plots (4 sorting methods)...\n")

      # 1. By Multimodal attribution sum
      cat("    * By Multimodal attribution sum\n")
      line_plot_m <- plot_top_cytobands(candidate_genes, top_n, plots_dir)
      plot_results$line_plot_multimodal <- line_plot_m

      # 2. By Difference (M - G)
      cat("    * By Difference (Multimodal - Gene-only)\n")
      line_plot_diff <- plot_top_cytobands_by_difference(candidate_genes, top_n, plots_dir)
      plot_results$line_plot_difference <- line_plot_diff

      # 3. By Ratio (M / G)
      cat("    * By Ratio (Multimodal / Gene-only)\n")
      line_plot_ratio <- plot_top_cytobands_by_ratio(candidate_genes, top_n, plots_dir)
      plot_results$line_plot_ratio <- line_plot_ratio
    }

    # Lollipop charts
    if (nrow(cytoband_summary) > 0) {
      cat("  - Creating lollipop charts...\n")
      lollipop_results <- plot_lollipop_charts(cytoband_summary, top_n, plots_dir)
      plot_results$lollipop <- lollipop_results
    }

    results$plots <- plot_results
    cat("✓ All visualizations completed\n")
    cat("Plots saved to:", plots_dir, "\n\n")
  }

  # Final summary
  cat("========== Analysis Summary ==========\n")
  cat("Total input genes:", nrow(gene_data), "\n")
  cat("Processed genes with location info:", sum(!is.na(processed_data$chromosome)), "\n")
  cat("Candidate genes:", nrow(candidate_genes), "\n")
  cat("Cytobands with candidates:", nrow(cytoband_summary), "\n")
  cat("Output directory:", output_dir, "\n")

  if (create_plots) {
    cat("Plots created:", paste(list.files(plots_dir, pattern = "\\.pdf$"), collapse = ", "), "\n")
  }
  cat("===================================\n")

  cat("\n🎉 Analysis completed successfully!\n")
  cat("Check the output directory for all results.\n")

  return(results)
}

#' Quick Analysis with Example Data
#'
#' Runs a quick analysis using built-in example data for testing and demonstration.
#'
#' @param output_dir Character string specifying output directory (default: "example_results")
#'
#' @return Analysis results from run_cytoband_analysis()
#'
#' @details
#' This function creates example data and runs the complete analysis pipeline.
#' Useful for testing the package functionality and exploring output formats.
#'
#' @examples
#' \dontrun{
#' # Run example analysis
#' example_results <- run_example_analysis()
#'
#' # Run with custom output directory
#' example_results <- run_example_analysis("test_output")
#' }
#'
#' @export
run_example_analysis <- function(output_dir = "example_results") {

  cat("Creating example data for demonstration...\n")

  # Create example gene data
  example_data <- data.frame(
    Gene = c("TP53", "BRCA1", "BRCA2", "MYC", "EGFR", "KRAS", "PIK3CA", "PTEN",
             "AKT1", "ERBB2", "CDH1", "ATM", "CHEK2", "PALB2", "RAD51"),
    Local.Index.G = c(0.523, 0.445, 0.678, 0.234, 0.567, 0.389, 0.456, 0.612,
                      0.334, 0.478, 0.289, 0.567, 0.423, 0.345, 0.234),
    Local.Index.M = c(0.687, 0.512, 0.891, 0.456, 0.623, 0.567, 0.389, 0.734,
                      0.456, 0.512, 0.345, 0.234, 0.567, 0.123, 0.567),
    Attribution.G = c(-0.234, 0.123, -0.345, 0.089, 0.234, -0.156, 0.267, -0.123,
                      0.178, -0.089, 0.123, -0.234, 0.156, -0.178, 0.089),
    Attribution.M = c(0.456, 0.234, 0.567, 0.178, 0.289, 0.345, 0.189, 0.423,
                      0.267, 0.156, 0.234, 0.123, 0.289, 0.156, 0.234),
    P.Value.G = c(0.002, 0.034, 0.0005, 0.056, 0.078, 0.023, 0.089, 0.005,
                  0.045, 0.098, 0.067, 0.034, 0.078, 0.023, 0.045),
    P.Value.M = c(0.001, 0.023, 0.0001, 0.045, 0.067, 0.012, 0.078, 0.003,
                  0.034, 0.089, 0.056, 0.023, 0.067, 0.012, 0.034),
    stringsAsFactors = FALSE
  )

  # Save example data
  example_file <- file.path(output_dir, "example_input.csv")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  write.csv(example_data, example_file, row.names = FALSE)

  cat("Example data saved to:", example_file, "\n")

  # Run analysis
  results <- run_cytoband_analysis(
    input_file = example_file,
    output_dir = output_dir,
    create_plots = TRUE
  )

  return(results)
}