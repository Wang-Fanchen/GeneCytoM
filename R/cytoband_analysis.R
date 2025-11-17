#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarize mutate arrange left_join select desc n
#'
#' Analyze Cytoband Location Attribution (Multimodal vs Gene-Only)
#' @importFrom magrittr %>%
#' Analyze Cytoband Location Attribution (Multimodal vs Gene-Only)
#'
#' Performs comprehensive cytoband-level analysis including aggregation,
#' difference and ratio calculations comparing multimodal and gene-only
#' attribution values.
#'
#' @param candidate_data Data frame of candidate genes from screen_candidate_genes()
#' @param output_dir Character string specifying output directory (default: ".")
#'
#' @return A data frame with cytoband-level summary statistics including:
#'   - Cytoband: Cytoband identifier
#'   - Attribution.G_sum: Sum of absolute attribution values for gene-only model
#'   - Attribution.M_sum: Sum of absolute attribution values for multimodal model
#'   - Difference: Attribution.M_sum - Attribution.G_sum (positive = multimodal advantage)
#'   - Ratio: Attribution.M_sum / Attribution.G_sum (>1 = multimodal advantage)
#'   - n_genes: Number of genes in each cytoband
#'
#' @details
#' The function performs the following steps:
#' 1. Calculate absolute values for attribution columns
#' 2. Aggregate genes by cytoband
#' 3. Calculate sums for both gene-only and multimodal attribution
#' 4. Compute differences and ratios to quantify multimodal advantage
#' 5. Sort by Attribution.M_sum (multimodal attribution) in descending order
#' 6. Save results to CSV file
#'
#' Positive differences and ratios > 1 indicate cytobands where the multimodal
#' model shows systematic advantage over the gene-only model.
#'
#' @examples
#' \dontrun{
#' # Analyze cytoband attribution for candidate genes
#' candidates <- screen_candidate_genes(processed_data)
#' summary <- analyze_cytoband_location(candidates, output_dir = "results")
#'
#' # View top cytobands
#' head(summary, 10)
#' }
#'
#' @export
analyze_cytoband_location <- function(candidate_data, output_dir = ".") {

  # Input validation
  if (!is.data.frame(candidate_data)) {
    stop("Input must be a data frame")
  }

  required_cols <- c("Cytoband", "Attribution.G", "Attribution.M")
  missing_cols <- setdiff(required_cols, colnames(candidate_data))

  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  cat("Analyzing", nrow(candidate_data), "candidate genes for cytoband attribution\n")

  # Add absolute value columns
  analysis_data <- candidate_data %>%
    mutate(
      Attribution.G_abs = abs(Attribution.G),
      Attribution.M_abs = abs(Attribution.M)
    )

  # Aggregate by cytoband
  cytoband_summary <- analysis_data %>%
    group_by(Cytoband) %>%
    summarize(
      Attribution.G_sum = sum(Attribution.G_abs, na.rm = TRUE),
      Attribution.M_sum = sum(Attribution.M_abs, na.rm = TRUE),
      Difference = Attribution.M_sum - Attribution.G_sum,
      Ratio = ifelse(Attribution.G_sum == 0,
                     NA_real_,
                     Attribution.M_sum / Attribution.G_sum),
      n_genes = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(Attribution.M_sum))

  # Report statistics
  cat("\n========== Cytoband Summary Statistics ==========\n")
  cat("Total cytobands with candidate genes:", nrow(cytoband_summary), "\n")
  cat("Cytobands with positive difference:", sum(cytoband_summary$Difference > 0, na.rm = TRUE), "\n")
  cat("Cytobands with ratio > 1:", sum(cytoband_summary$Ratio > 1, na.rm = TRUE), "\n")
  cat("Mean difference:", round(mean(cytoband_summary$Difference, na.rm = TRUE), 4), "\n")
  cat("Mean ratio:", round(mean(cytoband_summary$Ratio, na.rm = TRUE), 2), "\n")

  # Show top 10 cytobands
  cat("\nTop 10 Cytobands by Attribution.M_sum:\n")
  print(head(cytoband_summary, 10))
  cat("===============================================\n")

  # Save to file
  output_file <- file.path(output_dir, "Candidate_cytoband_attribution_summary.csv")
  write.csv(cytoband_summary, output_file, row.names = FALSE)

  cat("\nCytoband attribution summary saved to:", output_file, "\n")

  return(cytoband_summary)
}

#' Generate Comprehensive Cytoband Summary
#'
#' Creates a comprehensive 5-column summary table including all genes and candidate genes
#' for each cytoband, along with attribution differences.
#'
#' @param allgene_data Data frame with all processed genes from process_gene_data()
#' @param output_dir Character string specifying output directory (default: ".")
#'
#' @return A data frame with comprehensive cytoband summary containing:
#'   - Cytoband: Cytoband identifier
#'   - Genes: Comma-separated list of all genes in the cytoband
#'   - Candidate_Genes: Comma-separated list of candidate genes (if any)
#'   - All_Diff: Attribution difference for all genes
#'   - Cand_Diff: Attribution difference for candidate genes only
#'
#' @details
#' This function provides a comprehensive view by:
#' 1. Listing all genes in each cytoband
#' 2. Identifying which genes are candidates
#' 3. Calculating attribution differences for all genes
#' 4. Calculating attribution differences for candidate genes only
#' 5. Presenting results in a 5-column format for easy comparison
#'
#' @examples
#' \dontrun{
#' # Generate comprehensive summary
#' processed <- process_gene_data(gene_data)
#' comprehensive_summary <- generate_cytoband_summary(processed, output_dir = "results")
#'
#' # View cytobands with most candidate genes
#' summary_df <- comprehensive_summary
#' summary_df$n_candidates <- sapply(strsplit(summary_df$Candidate_Genes, ", "), length)
#' summary_df[order(-summary_df$n_candidates), ]
#' }
#'
#' @export
generate_cytoband_summary <- function(allgene_data, output_dir = ".") {

  # Input validation
  if (!is.data.frame(allgene_data)) {
    stop("Input must be a data frame")
  }

  required_cols <- c("Cytoband", "Gene", "Attribution.G", "Attribution.M", "Comparison")
  missing_cols <- setdiff(required_cols, colnames(allgene_data))

  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Check for optional P-value columns
  has_pvalue_g <- "P.Value.G" %in% colnames(allgene_data)
  has_pvalue_m <- "P.Value.M" %in% colnames(allgene_data)

  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  cat("Generating comprehensive cytoband summary from", nrow(allgene_data), "genes\n")

  # Add absolute value columns
  allgene_data <- allgene_data %>%
    mutate(
      Attribution.G_abs = abs(Attribution.G),
      Attribution.M_abs = abs(Attribution.M)
    )

  # Filter out rows with missing cytoband information
  allgene_cb <- allgene_data %>%
    filter(!is.na(Cytoband), Cytoband != "")

  cat("Genes with valid cytoband information:", nrow(allgene_cb), "\n")

  # 1. All genes in each cytoband
  cytoband_genes <- allgene_cb %>%
    group_by(Cytoband) %>%
    summarize(
      Genes = paste(sort(unique(Gene)), collapse = ", "),
      .groups = "drop"
    )

  # 2. Candidate genes in each cytoband (conditional P-value filtering)
  if (has_pvalue_m) {
    cytoband_candidates <- allgene_cb %>%
      filter(Comparison == "up", P.Value.M < 0.05) %>%
      group_by(Cytoband) %>%
      summarize(
        Candidate_Genes = paste(sort(unique(Gene)), collapse = ", "),
        .groups = "drop"
      )
  } else {
    cytoband_candidates <- allgene_cb %>%
      filter(Comparison == "up") %>%
      group_by(Cytoband) %>%
      summarize(
        Candidate_Genes = paste(sort(unique(Gene)), collapse = ", "),
        .groups = "drop"
      )
  }

  # 3. Attribution differences for all genes
  cyto_all_attr <- allgene_cb %>%
    group_by(Cytoband) %>%
    summarize(
      All_AttrG_sum = sum(Attribution.G_abs, na.rm = TRUE),
      All_AttrM_sum = sum(Attribution.M_abs, na.rm = TRUE),
      All_Diff = All_AttrM_sum - All_AttrG_sum,
      .groups = "drop"
    )

  # 4. Attribution differences for candidate genes (conditional P-value filtering)
  if (has_pvalue_m) {
    cyto_cand_attr <- allgene_cb %>%
      filter(Comparison == "up", P.Value.M < 0.05) %>%
      group_by(Cytoband) %>%
      summarize(
        Cand_AttrG_sum = sum(Attribution.G_abs, na.rm = TRUE),
        Cand_AttrM_sum = sum(Attribution.M_abs, na.rm = TRUE),
        Cand_Diff = Cand_AttrM_sum - Cand_AttrG_sum,
        .groups = "drop"
      )
  } else {
    cyto_cand_attr <- allgene_cb %>%
      filter(Comparison == "up") %>%
      group_by(Cytoband) %>%
      summarize(
        Cand_AttrG_sum = sum(Attribution.G_abs, na.rm = TRUE),
        Cand_AttrM_sum = sum(Attribution.M_abs, na.rm = TRUE),
        Cand_Diff = Cand_AttrM_sum - Cand_AttrG_sum,
        .groups = "drop"
      )
  }

  # 5. Combine into final 5-column table
  cytoband_summary_5 <- cytoband_genes %>%
    left_join(cytoband_candidates, by = "Cytoband") %>%
    mutate(Candidate_Genes = ifelse(is.na(Candidate_Genes), "", Candidate_Genes)) %>%
    left_join(cyto_all_attr %>% select(Cytoband, All_Diff), by = "Cytoband") %>%
    left_join(cyto_cand_attr %>% select(Cytoband, Cand_Diff), by = "Cytoband") %>%
    arrange(Cytoband)

  # Add statistics
  cytoband_summary_5$n_all_genes <- sapply(strsplit(cytoband_summary_5$Genes, ", "), length)
  cytoband_summary_5$n_candidate_genes <- sapply(strsplit(cytoband_summary_5$Candidate_Genes, ", "), length)

  # Report summary statistics
  cat("\n========== Comprehensive Cytoband Summary Statistics ==========\n")
  cat("Total cytobands:", nrow(cytoband_summary_5), "\n")
  cat("Cytobands with candidate genes:", sum(cytoband_summary_5$n_candidate_genes > 0), "\n")
  cat("Total all genes:", sum(cytoband_summary_5$n_all_genes), "\n")
  cat("Total candidate genes:", sum(cytoband_summary_5$n_candidate_genes), "\n")

  # Show cytobands with most candidate genes
  top_candidates <- cytoband_summary_5[order(-cytoband_summary_5$n_candidate_genes), ]
  cat("\nTop 5 cytobands by number of candidate genes:\n")
  print(head(top_candidates[, c("Cytoband", "n_all_genes", "n_candidate_genes", "All_Diff", "Cand_Diff")], 5))
  cat("==========================================================\n")

  # Save to file
  output_file <- file.path(output_dir, "Cytoband_summary_5cols.csv")
  write.csv(cytoband_summary_5, output_file, row.names = FALSE)

  cat("\nComprehensive cytoband summary saved to:", output_file, "\n")

  return(cytoband_summary_5)
}
