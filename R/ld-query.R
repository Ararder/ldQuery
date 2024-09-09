query_sv <- function(tbl, type = "SNV") {
  stopifnot(all(c("POS", "CHR") %in% colnames(tbl)))
  stopifnot(type %in% c("SNV", "SV"))
  if(type == "SV") {
    dset <- arrow::open_dataset(Sys.getenv("ld_ref_sv"))
  } else{
    dset <- arrow::open_dataset(Sys.getenv("ld_ref_snv"))
  }

  split_df <- split(tbl, tbl$CHR)
  res <- purrr::imap(split_df, \(snps,chrom) single_query(chr = as.integer(chrom), pos = snps[["POS"]], dset = dset))

  purrr::list_rbind(res)


}


query_topmed <- function(chr, pos, ld_ref) {
  rlang::check_required(chr)
  rlang::check_required(pos)
  stopifnot(length(chr) == length(pos))
  # type <- rlang::arg_match(type)

  # check which LD data to query
  pos <- as.integer(pos)
  dset <- arrow::open_dataset(ld_ref)


  single_query(chr = chr , pos = pos, dset = dset)


}

single_query <- function(chr, pos, dset) {

  out <- vector("list", length = length(pos))
  if(length(pos) <= 10) {
    for(idx in seq_along(pos)) {
      out[[idx]] <-
        dset |>
        dplyr::filter(CHR == {{ chr }}) |>
        dplyr::filter(SNP1 == pos[idx] | SNP2 == pos[idx]) |>
        dplyr::collect()
    }
    purrr::list_rbind(out) |>
      dplyr::mutate(
        tmp = dplyr::if_else(SNP1 %in% pos, SNP1, SNP2),
        SNP2 = dplyr::if_else(SNP2 %in% pos, SNP1, SNP2),
        SNP1 = tmp
      ) |>
      dplyr::select(-tmp)
  } else {


    dset |>
      dplyr::filter(CHR == {{ chr }}) |>
      dplyr::filter(SNP1 %in% pos | SNP2 %in% pos) |>
      dplyr::mutate(
        tmp = dplyr::if_else(SNP1 %in% pos, SNP1, SNP2),
        SNP2 = dplyr::if_else(SNP2 %in% pos, SNP1, SNP2),
        SNP1 = tmp
      ) |>
      dplyr::select(-tmp) |>
      dplyr::collect()
  }
}
