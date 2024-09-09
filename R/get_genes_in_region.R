get_genes_in_region <- function(region_chr, region_start, region_end, genetype, build = c("38", "37"), gene_df) {
  build <- rlang::arg_match(build)
  stopifnot("chr should be provided in UCSC format: chr11 not 11" = is.character(region_chr))


  if(!missing(genetype)) {
    stopifnot("enter a character vector to filter for genetypes" = is.character(genetype))
    gene_matrix <- dplyr::filter(gene_matrix, gene_type %in% {{ genetype }})
  }

  if(build == "38") {
    gene_df <- dplyr::rename(
      gene_df,
      chr = chr_38,
      start = start_38,
      end = end_38,
      strand = strand_38
    )
  } else {
    gene_df <- dplyr::rename(
      gene_df,
      chr = chr_37,
      start = start_37,
      end = end_37,
      strand = strand_37
    )

  }

  # There are four possible scenarios for overlap:
  # ----------------------------------------------------------------------------
  # 1. A gene ends inside the region, and starts outside
  # 2. A gene starts inside the region, and ends outside
  # 3. A gene both starts outside the region, and and ends outside it it. It spans the full region
  # 4. A gene is starts AND ends inside the region. it is contained inside it.
  # Scenario four is automatically detected by filters for (1) and (2)



  # -------------------------------------------------------------------------
  # Scenario (1) - The gene ends inside the region
  #   |------------->
  # ------------|_____________________________________________|------
  #.      start|                                             |end

  right <- gene_df |>
    dplyr::filter(chr == {{ region_chr }}) |>
    dplyr::filter(end >= {{ region_start }} & end <={{ region_end }})

  # -------------------------------------------------------------------------
  # Scenario (2) - Starts inside gene boundary
  #                  |---------------------------------------------->
  # -----------|_____________________________________________|------
  #.      start|                                             |end
  left <- gene_df |>
    dplyr::filter(chr == {{ region_chr }}) |>
    dplyr::filter(start >= {{ region_start }} & start <= {{ region_end }} )





  # -------------------------------------------------------------------------
  # Scenario (3) - The gene overlaps the full region
  #   |---------------------------------------------------------------->
  # ------------|_____________________________________________|------
  #        start|                                             |end

  spanning <- gene_df |>
    dplyr::filter(chr == {{ region_chr }}) |>
    dplyr::filter(start <= {{ region_start }} & end >= {{ region_end }})


  dplyr::bind_rows(left, right, spanning) |>
    dplyr::distinct(gene_name, .keep_all = TRUE)

}
