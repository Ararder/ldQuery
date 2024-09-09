test_that("multiplication works", {
  skip()
  bed_file = "~/arvhar/nordic_mdd_gwas/results/clumping/broad_mdd/n_loci.bedtools"
  gwas <- arrow::read_tsv_arrow("~/arvhar/nordic_mdd_gwas/sumstats/broad_mdd.gz")
  bed <- read_bed(bed_file)
  gene_df <- arrow::read_parquet("inst/extdata/gene_df.parquet")

  locus_zoom


})
