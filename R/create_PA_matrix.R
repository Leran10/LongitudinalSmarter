#' Create presence/absence matrix from count data
#'
#' Converts a count matrix to presence/absence (0/1) using TMM normalization
#' and CPM thresholding. A feature is considered "present" if its normalized
#' CPM >= target_cpm AND raw count >= min_reads_floor.
#'
#' @param count_matrix Matrix or data.frame of count data (features x samples).
#'   Rownames should be feature IDs, colnames should be sample IDs.
#' @param target_cpm CPM threshold for presence (default: 0.5)
#' @param min_reads_floor Minimum raw read count threshold (default: 3)
#' @param normalization_method edgeR normalization method (default: "TMMwsp" for sparse data)
#'
#' @return A numeric matrix (0/1) with same dimensions and dimnames as input
#'
#' @details
#' This function uses edgeR's TMM (trimmed mean of M-values) normalization,
#' which is robust for sparse microbiome data. The "TMMwsp" method is
#' specifically designed for data with many zeros.
#'
#' The presence rule requires BOTH conditions:
#' - Normalized CPM >= target_cpm
#' - Raw count >= min_reads_floor
#'
#' This guards against calling rare features "present" due to normalization
#' artifacts.
#'
#' @export
#'
#' @importFrom edgeR DGEList calcNormFactors cpm
#'
#' @examples
#' \dontrun{
#' # Create PA matrix from counts
#' PA_matrix <- create_PA_matrix(
#'   count_matrix = my_counts,
#'   target_cpm = 0.5,
#'   min_reads_floor = 3
#' )
#'
#' # More stringent threshold
#' PA_strict <- create_PA_matrix(
#'   count_matrix = my_counts,
#'   target_cpm = 1.0,
#'   min_reads_floor = 5
#' )
#' }
create_PA_matrix <- function(count_matrix,
                             target_cpm = 0.5,
                             min_reads_floor = 3L,
                             normalization_method = "TMMwsp") {

  # Input validation
  if (!is.matrix(count_matrix) && !is.data.frame(count_matrix)) {
    stop("count_matrix must be a matrix or data.frame")
  }

  if (target_cpm <= 0) {
    stop("target_cpm must be > 0")
  }

  if (min_reads_floor < 0) {
    stop("min_reads_floor must be >= 0")
  }

  # Convert to matrix if data.frame
  if (is.data.frame(count_matrix)) {
    rn <- rownames(count_matrix)
    cn <- colnames(count_matrix)
    count_matrix <- as.matrix(count_matrix)
    rownames(count_matrix) <- rn
    colnames(count_matrix) <- cn
  }

  # Check for negative values
  if (any(count_matrix < 0, na.rm = TRUE)) {
    stop("count_matrix contains negative values")
  }

  # Remove features with all zeros
  feature_sums <- rowSums(count_matrix)
  if (all(feature_sums == 0)) {
    stop("All features have zero counts")
  }

  n_zeros <- sum(feature_sums == 0)
  if (n_zeros > 0) {
    message("Removing ", n_zeros, " features with all zero counts")
    count_matrix <- count_matrix[feature_sums > 0, , drop = FALSE]
  }

  # Build DGEList
  y <- edgeR::DGEList(counts = count_matrix)

  # Normalize using TMM (or specified method)
  y <- edgeR::calcNormFactors(y, method = normalization_method)

  # Calculate normalized CPMs
  cpm_mat <- edgeR::cpm(y, normalized.lib.sizes = TRUE)

  # Apply presence rule: CPM >= threshold AND raw count >= floor
  present_logical <- (cpm_mat >= target_cpm) & (count_matrix >= min_reads_floor)

  # Convert to 0/1 matrix
  PA <- matrix(
    as.integer(present_logical),
    nrow = nrow(present_logical),
    ncol = ncol(present_logical),
    dimnames = dimnames(present_logical)
  )

  # Report statistics
  n_present <- sum(PA)
  n_total <- length(PA)
  pct_present <- round(100 * n_present / n_total, 2)

  message(
    "Created PA matrix: ",
    nrow(PA), " features x ", ncol(PA), " samples\n",
    "  Present cells: ", n_present, " / ", n_total, " (", pct_present, "%)\n",
    "  Thresholds: CPM >= ", target_cpm, " & count >= ", min_reads_floor
  )

  return(PA)
}
