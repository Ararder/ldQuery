utils::globalVariables(c(
  "start", "end", "RSID", "POS", "P",
  "chr_38", "start_38", "chr_37", "start_37",
  "chr", "CHR","f0", "f1", "f2", "f3", "f4", "f5", "tmp",
  "strand_37", "strand_38", "gene_name", "gene_type", "strand",
  "line_start","line_end", "top_snps", "R2", "end_37","end_38","SNP1", "SNP2",
  "gene_df"
))


#' Plot a GWAS loci colored by ld with nearby SNPs
#'
#' @param gwas a summary statistics file in tidyGWAS format
#' @param bed_file a bed file with six columns
#' @param build genome build of GWAS and bedfile
#' @param expand how much to expand sides of loci
#' @param ld_path filepath to LD reference data
#'
#' @return a [[ggplot2::ggplot()]] figure
#' @export
#'
#' @examples \dontrun{
#' locus_zoom(gwas, bedfile, build = "38", ld_path ="path/to/topmed_ld")
#' }
locus_zoom <- function(gwas, bed_file, build = c("38","37"), expand=100000L, ld_path) {

  build <- rlang::arg_match(build)
  rlang::check_required(gwas)
  rlang::check_required(expand)
  rlang::check_required(gene_df)
  rlang::check_required(ld_path)
  stopifnot("Can only plot 1 loci at a time" = nrow(bed_file) == 1)
  gene_df <- arrow::read_parquet(system.file("extdata", "gene_df.parquet", package = "ldQuery"))


  locus_data <- merge_bed_with_gwas(gwas = gwas, bed = bed_file)

  # extract required data
  # ----------------------------------------------------------------------------

  ld_info <- query_topmed(chr = locus_data[["chr_numeric"]], pos = locus_data[["pos"]], ld_ref = ld_path)
  plot_start_limit <- as.integer(locus_data[["start"]] - expand)
  plot_end_limit <- as.integer(locus_data[["end"]] + expand)

  pdf <- dplyr::filter(gwas, CHR == locus_data[["chr_numeric"]] & POS >= {{ plot_start_limit }} & POS <= {{ plot_end_limit }}) |>
    dplyr::filter(P < 0.005) |>
    dplyr::left_join(dplyr::select(ld_info, SNP2, R2, CHR), by = c("POS" = "SNP2", "CHR"))



  # plot SNP figure ---------------------------------------------------------

  x_scale <- ggplot2::scale_x_continuous(breaks = scales::breaks_extended(5), limits = c(plot_start_limit, plot_end_limit))

  log10_snp_plot <- create_snp_figure(
    pdf = pdf,
    plot_start_limit = plot_start_limit,
    plot_end_limit = plot_end_limit,
    locus_data = locus_data
  ) +
    x_scale +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())



  # plot gene figure --------------------------------------------------------

  plot_genes <- create_gene_figure(
    region_chr = locus_data[["chr"]],
    region_start = plot_start_limit,
    region_end = plot_end_limit,
    build = build,
    gene_df = gene_df
  ) +
    x_scale

  cowplot::plot_grid(log10_snp_plot, plot_genes, ncol = 1, align = "v", vjust= 0, axis="lr", hjust=-2)
}



create_snp_figure <- function(pdf, plot_start_limit,plot_end_limit, locus_data) {

  list2env(locus_data,  envir = environment())
  labels <- dplyr::filter(pdf, RSID %in% top_snps) |>
    dplyr::mutate(R2 = NA_real_)

  pdf |>
    ggplot2::ggplot(ggplot2::aes(POS, -log10(P), color = R2)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient(low = "blue", high = "red", limits=c(0.2,1)) +
    ggplot2::labs(
      x = "",
      title = glue::glue("Chromosome {chr_numeric}, {plot_start_limit} - {plot_end_limit}")
    ) +
    ggplot2::theme_light() +
    ggrepel::geom_text_repel(data = labels, ggplot2::aes(label = RSID),  min.segment.length = ggplot2::unit(0, 'lines')) +
    ggplot2::geom_hline(yintercept = -log10(5e-08), linetype = "dashed", alpha = 0.3)


}


#' Display genes based on position in the genome
#'
#' @param region_chr Chromosome, UCSC format: "chr11" not 11
#' @param region_start start of region
#' @param region_end end of region
#' @param plot_start_limit start of figure, in genomic location
#' @param plot_end_limit end of figure, in genomic location
#' @param build GRCh38 or GRCh37 to map genes to regions?
#' @param add_region input a tibble in the same format to add a custom region to draw
#' @param gene_type_in_name append the gene_type to the y-axis label on the gene plot
#' @param ... pass additional arguments to get_genes_in_region, for example gene_type = "protein_coding"
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples \dontrun{
#' create_gene_figure(chr = "chr12", 10, 100000, build = 37)
#' }
create_gene_figure <- function(region_chr, region_start, region_end, plot_start_limit, plot_end_limit, build, add_region, gene_type_in_name=TRUE, ...) {
  stopifnot("input start and end limits as integers" = all(is.integer(region_start), is.integer(region_end)))
  stopifnot("chr in UCSC format, 'chr11', not 11" = is.character(region_chr))

  if(missing(plot_start_limit) & missing(plot_end_limit)) {
    plot_start_limit <- region_start
    plot_end_limit <- region_end
  }

  # Find all genes that span the genomic region
  gene_data <- get_genes_in_region(
    region_chr = region_chr,
    region_start = region_start,
    region_end = region_end,
    build = build,
    ...
  )

  if(!missing(add_region)) {
    gene_data <-   dplyr::bind_rows(gene_data,add_region)
  }
  if(gene_type_in_name){
    gene_data <- gene_data |>
      dplyr::mutate(gene_name = paste(gene_name, gene_type, sep = " - "))
  }



  gene_data <- gene_data |>
    dplyr::mutate(
      # if the gene boundaries are outside the loci, set them to the loci boundaries
      start  = dplyr::if_else(start <= plot_start_limit, plot_start_limit, start),
      end    = dplyr::if_else(end >= plot_end_limit, plot_end_limit, end),
      # use strand information to determine direction of gene by flipping
      # end and start
      line_start = dplyr::if_else(strand == "+", start, end),
      line_end = dplyr::if_else(strand == "+", end, start)
    )

  gene_data |>
    ggplot2::ggplot(ggplot2::aes(x = line_start, xend = line_end,  y = gene_name, yend = gene_name)) +
    ggplot2::geom_segment(arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm"))) +
    ggplot2::scale_x_continuous(limits = c(plot_start_limit, plot_end_limit)) +
    ggplot2::labs(
      y = "",
      x = "Position"
    ) +
    ggplot2::theme_light()

}


