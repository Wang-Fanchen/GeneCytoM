#' @importFrom magrittr %>%
#' @importFrom dplyr mutate filter select left_join
#'
#' Process Gene Data with Location Information
#' @importFrom magrittr %>%
#' Process Gene Data with Location Information
#'
#' Integrates gene attribution data from multimodal and gene-only models with
#' chromosomal location information. Adds comparison column (multimodal vs gene-only)
#' and matches genes to their cytoband locations.
#'
#' @param gene_data A data frame containing gene analysis data with required columns:
#'   Gene, Local.Index.G, Local.Index.M, Attribution.G, Attribution.M
#'   - G suffix: Gene-only model metrics
#'   - M suffix: Multimodal model metrics (gene + image)
#'   Optional columns: P.Value.G, P.Value.M
#' @param location_data Optional data frame with gene location information.
#'   If NULL, uses built-in gene_location data.
#'
#' @return A data frame with added location information and comparison column.
#'   Includes new columns: Comparison, Cytoband, chromosome, start_position_on_the_genomic_accession,
#'   end_position_on_the_genomic_accession, orientation, exon_count, map_location
#'
#' @details
#' The function performs the following steps:
#' 1. Adds Comparison column based on Local.Index differences (multimodal vs gene-only)
#' 2. Matches genes to chromosome locations using Symbol names
#' 3. For unmatched genes, attempts matching through Aliases
#' 4. Reports matching statistics
#' 5. Creates map_location column preserving original cytoband information
#'
#' The Comparison column categorizes genes as:
#' - "up": Multimodal model shows greater importance than gene-only
#' - "down": Gene-only model shows greater importance than multimodal
#' - "mid": Equal importance between models
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data <- read.csv("your_gene_data.csv", header = TRUE, check.names = FALSE)
#' processed <- process_gene_data(data)
#' }
#'
#' @export
process_gene_data <- function(gene_data, location_data = NULL) {

  # Input validation
  required_cols <- c("Gene", "Local.Index.G", "Local.Index.M",
                    "Attribution.G", "Attribution.M")
  missing_cols <- setdiff(required_cols, colnames(gene_data))

  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Check for optional P-value columns
  has_pvalue_g <- "P.Value.G" %in% colnames(gene_data)
  has_pvalue_m <- "P.Value.M" %in% colnames(gene_data)

  if (has_pvalue_g) {
    cat("✓ P.Value.G column detected\n")
  } else {
    cat("ℹ P.Value.G column not found - P-value filtering for condition G will be skipped\n")
  }

  if (has_pvalue_m) {
    cat("✓ P.Value.M column detected\n")
  } else {
    cat("ℹ P.Value.M column not found - P-value filtering for condition M will be skipped\n")
  }

  # Use built-in location data if not provided
  if (is.null(location_data)) {
    # Try to load from package first, fallback to global environment
    tryCatch({
      data("gene_location", package = "GeneCytoM", envir = environment())
      location_data <- gene_location
    }, error = function(e) {
      # Fallback to global environment if package not properly installed
      if (exists("gene_location")) {
        location_data <- get("gene_location")
      } else {
        stop("gene_location data not found. Please ensure GeneCytoM package is properly installed.")
      }
    })
  }

  # Step 1: Add Comparison column
  processed_data <- gene_data %>%
    mutate(Comparison = ifelse(Local.Index.M - Local.Index.G > 0, "down",
                              ifelse(Local.Index.M - Local.Index.G == 0, "mid", "up")))

  # Step 2: Match by Symbol column
  processed_data <- merge(processed_data,
                          location_data[, c("Symbol", "Cytoband", "chromosome",
                                           "start_position_on_the_genomic_accession",
                                           "end_position_on_the_genomic_accession",
                                           "orientation", "exon_count")],
                          by.x = "Gene", by.y = "Symbol", all.x = TRUE)

  # Track Symbol matching
  matched_symbol <- sum(!is.na(processed_data$chromosome))

  # Step 3: Match unmatched genes by Aliases
  unmatched <- is.na(processed_data$chromosome)
  matched_aliases <- 0

  if (any(unmatched)) {
    for (i in which(unmatched)) {
      gene <- processed_data$Gene[i]
      match_idx <- grep(paste0("\\b", gene, "\\b"), location_data$Aliases)
      if (length(match_idx) > 0) {
        processed_data[i, c("Cytoband", "chromosome", "start_position_on_the_genomic_accession",
                            "end_position_on_the_genomic_accession", "orientation", "exon_count")] <-
          location_data[match_idx[1], c("Cytoband", "chromosome", "start_position_on_the_genomic_accession",
                                       "end_position_on_the_genomic_accession", "orientation", "exon_count")]
        matched_aliases <- matched_aliases + 1
      }
    }
  }

  # Step 4: Report matching statistics
  total_genes <- nrow(processed_data)
  final_matched <- sum(!is.na(processed_data$chromosome))
  final_unmatched <- sum(is.na(processed_data$chromosome))

  cat("========== Gene Matching Results ==========\n")
  cat("Total genes:", total_genes, "\n")
  cat("Matched by Symbol:", matched_symbol, "\n")
  cat("Matched by Aliases:", matched_aliases, "\n")
  cat("Total matched:", final_matched, "\n")
  cat("Unmatched:", final_unmatched, "\n")
  cat("Matching success rate:", round(final_matched/total_genes*100, 2), "%\n")
  cat("======================================\n")

  # Step 5: Create map_location column
  processed_data$map_location <- processed_data$Cytoband
  processed_data$Cytoband <- gsub("-.*", "", processed_data$Cytoband)

  return(processed_data)
}

#' Validate Gene Data Input
#'
#' Internal function to validate input data format and structure.
#'
#' @param gene_data Data frame to validate
#'
#' @return TRUE if valid, throws error if invalid
#'
#' @keywords internal
validate_gene_data <- function(gene_data) {
  if (!is.data.frame(gene_data)) {
    stop("Input must be a data frame")
  }

  required_cols <- c("Gene", "Local.Index.G", "Local.Index.M",
                    "Attribution.G", "Attribution.M")
  missing_cols <- setdiff(required_cols, colnames(gene_data))

  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Check for empty gene names
  if (any(is.na(gene_data$Gene) | gene_data$Gene == "")) {
    warning("Some rows have empty gene names and will be excluded from analysis")
  }

  return(TRUE)
}
