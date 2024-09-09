test_that("Get genes work", {
  skip()
  gene_df <- arrow::read_parquet("inst/extdata/gene_df.parquet")

  expect_no_error(

    get_genes_in_region(
      region_chr = "chr1",
      region_start = 1,
      region_end = 10^5L,
      build = "38",
      gene_df = gene_df
    )
  )
})


test_that("locus zoom works", {
  skip()

  bed_file <- read_bed("~/arvhar/nordic_mdd_gwas/results/clumping/broad_mdd/n_loci.bedtools") |>
    dplyr::slice(10)
  gwas <- arrow::read_tsv_arrow("~/arvhar/nordic_mdd_gwas/sumstats/broad_mdd.gz")
  gene_df <- arrow::read_parquet("inst/extdata/gene_df.parquet")


  zm <- locus_zoom(
    gwas = gwas,
    bed_file = bed_file,
    build = "38",
    gene_df = gene_df,
    ld_path = "~/arvhar/snp_level_annotatations/ld_ref_topmed/ld_single_nucleotide_variants/"
  )



})
