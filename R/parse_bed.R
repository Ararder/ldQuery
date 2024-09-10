#' Read a bed file into a tibble
#'
#' @param path filepath
#' @param ... arguments to pass to [arrow::read_tsv_arrow()]
#'
#' @return a tibble
#' @export
#'
#' @examples \dontrun{
#' read_bed("file.bed")
#' }
read_bed <- function(path, ...) {
  arrow::read_tsv_arrow(path, col_names = FALSE, ...) |>
    dplyr::rename(
      chr = f0,
      start = f1,
      end = f2,
      n_snps = f3,
      top_pval = f4,
      top_rsid = f5
    )
}


merge_bed_with_gwas <- function(bed, gwas) {

  start <- dplyr::pull(bed, "start")
  end <- dplyr::pull(bed, "end")
  chr <- dplyr::pull(bed, "chr")
  chr_numeric <- as.integer(stringr::str_remove(chr, "chr"))

  # if multiple snps per loci, snps and pvals are entered with a "," separator
  top_pvals <- stringr::str_split(bed[[5]], ",")[[1]]
  top_snps <- stringr::str_split(bed[[6]], ",")[[1]]
  if(length(top_snps) > 1) {
    idx_top <- which.min((as.numeric(top_pvals)))
    snps <- top_snps[idx_top]
    pvals <- top_pvals[idx_top]

  } else {
    snps <- top_snps
    pvals <- top_pvals

  }

  pos <- gwas |>
    dplyr::filter(CHR == chr_numeric) |>
    dplyr::filter(RSID == snps) |>
    dplyr::pull(POS)

  list(
    start = start,
    end = end,
    chr = chr,
    chr_numeric = chr_numeric,
    top_pvals = top_pvals,
    pvals = pvals,
    top_snps = top_snps,
    snps = snps,
    pos = pos

  )
}
