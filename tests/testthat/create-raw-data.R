gene_matrix <- arrow::read_tsv_arrow("~/geneMatrix.tsv.gz")
gene_matrix |>
  dplyr::select(ensgid, gene_name, dup_gene_name, gene_type, gene_id,hg38h0, h1, h2, hstr, hg19g0, g1, g2, gstr) |>
  dplyr::tibble() |>
  dplyr::rename(
    chr_38 = hg38h0, start_38 = h1, end_38 = h2, strand_38 = hstr,
    chr_37 = hg19g0, start_37 = g1, end_37 = g2, strand_37 = gstr
  ) |>
  arrow::write_parquet("inst/extdata/gene_df.parquet")

