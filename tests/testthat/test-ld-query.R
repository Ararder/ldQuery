test_that("ld query ", {
  skip("Only runs on longleaf")

  path <- "~/arvhar/snp_level_annotatations/ld_ref_topmed/ld_single_nucleotide_variants/"
  ds <- arrow::open_dataset(path)



  chr <- 7
  pos <- 24398341
  query_topmed(chr, pos, path)

})
