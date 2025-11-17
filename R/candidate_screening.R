#' @importFrom magrittr %>%
#' @importFrom dplyr filter mutate
#'

#' Screen Candidate Genes (Multimodal-Advantaged)
#'
#' Filters genes based on comparison type and p-value thresholds to identify
#' candidate genes showing enhanced importance in the multimodal model compared
#' to the gene-only model.
#'
#' @param processed_data Data frame output from process_gene_data()
#' @param comparison_filter Character string specifying which comparison type to include:
#'   "up" (multimodal > gene-only), "down" (multimodal < gene-only),
#'   "mid" (equal), or "all" (default: "up")
#' @param pvalue_threshold_g Numeric p-value threshold for P.Value.G filtering (default: NULL).
#'   If NULL or if P.Value.G column doesn't exist, this filter is skipped.
#' @param pvalue_threshold_m Numeric p-value threshold for P.Value.M filtering (default: NULL).
#'   If NULL or if P.Value.M column doesn't exist, this filter is skipped.
#'
#' @return A data frame containing only the candidate genes that meet the filtering criteria.
#'   Includes all original columns plus any location information.
#'
#' @details
#' The function filters genes based on:
#' 1. Comparison type (whether multimodal shows advantage over gene-only)
#' 2. P-value thresholds for both models (optional)
#' 3. Reports filtering statistics
#'
#' Typical usage is comparison_filter = "up" to find genes where the multimodal
#' model (gene + image) shows greater importance than the gene-only model.
#'
#' @examples
#' \dontrun{
#' # Process data first
#' processed <- process_gene_data(gene_data)
#'
#' # Screen for multimodal-advantaged genes with P-value filtering
#' candidates <- screen_candidate_genes(
#'   processed,
#'   comparison_filter = "up",
#'   pvalue_threshold_g = 0.05,
#'   pvalue_threshold_m = 0.01
#' )
#'
#' # Screen without P-value filtering
#' all_up <- screen_candidate_genes(processed, comparison_filter = "up")
#'
#' # Screen with only multimodal P-value filtering
#' candidates_m <- screen_candidate_genes(processed, pvalue_threshold_m = 0.05)
#' }
#'
#' @export
screen_candidate_genes <- function(processed_data,
                                  comparison_filter = "up",
                                  pvalue_threshold_g = NULL,
                                  pvalue_threshold_m = NULL) {

  # Input validation
  if (!is.data.frame(processed_data)) {
    stop("Input must be a data frame")
  }

  if (!"Comparison" %in% colnames(processed_data)) {
    stop("Input data must contain 'Comparison' column. Run process_gene_data() first.")
  }

  valid_filters <- c("up", "down", "mid", "all")
  if (!comparison_filter %in% valid_filters) {
    stop("comparison_filter must be one of: ", paste(valid_filters, collapse = ", "))
  }

  # Check if P-value columns exist
  has_pvalue_g <- "P.Value.G" %in% colnames(processed_data)
  has_pvalue_m <- "P.Value.M" %in% colnames(processed_data)

  total_genes <- nrow(processed_data)
  filtered_data <- processed_data

  # Apply P.Value.G filter if column exists AND threshold provided
  if (has_pvalue_g && !is.null(pvalue_threshold_g)) {
    filtered_data <- filtered_data[!is.na(filtered_data$P.Value.G) &
                                   filtered_data$P.Value.G < pvalue_threshold_g, ]
    cat("Applied P.Value.G filter (<", pvalue_threshold_g, "):", nrow(filtered_data), "genes remaining\n")
  }

  # Apply P.Value.M filter if column exists AND threshold provided
  if (has_pvalue_m && !is.null(pvalue_threshold_m)) {
    filtered_data <- filtered_data[!is.na(filtered_data$P.Value.M) &
                                   filtered_data$P.Value.M < pvalue_threshold_m, ]
    cat("Applied P.Value.M filter (<", pvalue_threshold_m, "):", nrow(filtered_data), "genes remaining\n")
  }

  # Apply comparison filter
  if (comparison_filter == "all") {
    candidates <- filtered_data
  } else {
    candidates <- filtered_data[filtered_data$Comparison == comparison_filter, ]
  }

  # Report filtering statistics
  cat("========== Candidate Gene Screening ==========\n")
  cat("Total genes:", total_genes, "\n")

  if (has_pvalue_g && !is.null(pvalue_threshold_g)) {
    cat("P.Value.G filtering: APPLIED (threshold =", pvalue_threshold_g, ")\n")
  } else {
    cat("P.Value.G filtering: SKIPPED\n")
  }

  if (has_pvalue_m && !is.null(pvalue_threshold_m)) {
    cat("P.Value.M filtering: APPLIED (threshold =", pvalue_threshold_m, ")\n")
  } else {
    cat("P.Value.M filtering: SKIPPED\n")
  }

  if (comparison_filter != "all") {
    cat("Comparison filter:", comparison_filter, "\n")
  }

  cat("Final candidate genes:", nrow(candidates), "\n")
  cat("==============================================\n")

  return(candidates)
}

#' Calculate Filtering Statistics
#'
#' Internal function to calculate detailed filtering statistics.
#'
#' @param processed_data Original processed data
#' @param candidates Filtered candidate genes
#' @param comparison_filter Comparison filter used
#' @param pvalue_threshold_g P-value threshold used for P.Value.G
#' @param pvalue_threshold_m P-value threshold used for P.Value.M
#'
#' @return List with detailed statistics
#'
#' @keywords internal
calculate_filtering_stats <- function(processed_data, candidates, comparison_filter, pvalue_threshold_g = NULL, pvalue_threshold_m = NULL) {

  total_genes <- nrow(processed_data)
  candidate_count <- nrow(candidates)

  # Count by comparison type in original data
  comparison_counts <- table(processed_data$Comparison, useNA = "ifany")

  # Check if P-value columns exist
  has_pvalue_g <- "P.Value.G" %in% colnames(processed_data)
  has_pvalue_m <- "P.Value.M" %in% colnames(processed_data)

  # Count significant genes by comparison type for G if applicable
  sig_counts_g <- NULL
  if (has_pvalue_g && !is.null(pvalue_threshold_g)) {
    sig_genes_g <- processed_data[!is.na(processed_data$P.Value.G) &
                                processed_data$P.Value.G < pvalue_threshold_g, ]
    sig_counts_g <- table(sig_genes_g$Comparison, useNA = "ifany")
  }

  # Count significant genes by comparison type for M if applicable
  sig_counts_m <- NULL
  if (has_pvalue_m && !is.null(pvalue_threshold_m)) {
    sig_genes_m <- processed_data[!is.na(processed_data$P.Value.M) &
                                processed_data$P.Value.M < pvalue_threshold_m, ]
    sig_counts_m <- table(sig_genes_m$Comparison, useNA = "ifany")
  }

  stats <- list(
    total_genes = total_genes,
    candidate_genes = candidate_count,
    pvalue_threshold_g = pvalue_threshold_g,
    pvalue_threshold_m = pvalue_threshold_m,
    comparison_filter = comparison_filter,
    comparison_distribution = comparison_counts,
    significant_distribution_g = sig_counts_g,
    significant_distribution_m = sig_counts_m,
    has_pvalue_g = has_pvalue_g,
    has_pvalue_m = has_pvalue_m
  )

  return(stats)
}
