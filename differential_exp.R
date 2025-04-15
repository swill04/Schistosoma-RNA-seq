library(tidyverse)
library(DESeq2)
library(broom)

star_counts_df <- read_tsv("/data/classes/2025/spring/biol443/course_files/rnaseq_data/counts.tsv", comment = "#") |>
             mutate(across(where(is.numeric), as.integer))

star_counts_summary <- star_counts_df |>
    select(Geneid, contains('dedup/star')) |>
    rename_with(~str_remove(., "dedup/star/"), everything()) |>
    rename_with(~str_remove(., "_S[0-9]{2}_L005.bam:.*"), everything()) |>
    rowwise() |>
    mutate(total_counts = sum(c_across(where(is.numeric)), na.rm = T)) |>
    filter(total_counts >= 10)

star_sample_summary <- star_counts_df |>
    select(Geneid, contains('dedup/star')) |>
    rename_with(~str_remove(., "dedup/star/"), everything()) |>
    rename_with(~str_remove(., "_S[0-9]{2}_L005.bam:.*"), everything()) |>
    pivot_longer(-Geneid, names_to = 'sample', values_to = 'count') |>
    filter(count > 0) |>
    group_by(Geneid) |>
    tally() |>
    filter(n <= 3)

genes_to_remove = star_sample_summary$Geneid

star_counts_filt <- star_counts_summary |>
    filter(!Geneid %in% genes_to_remove) |>
    arrange(Geneid) |>
    select(-total_counts)

star_counts_m <- star_counts_filt |>
    select(-Geneid) |>
    as.matrix()
rownames(star_counts_m) <- star_counts_filt$Geneid

dists <- dist(t(star_counts_m))

star_dists_df <- as.matrix(dists) |>
    as_tibble(rownames = 'sample')

star_dist_plot <- star_dists_df |>
    pivot_longer(-sample, names_to = 'comp', values_to = 'dist') |>
    ggplot(aes(x = sample, y = comp, fill = dist)) +
    geom_tile() +
    scale_fill_viridis_c() +
    coord_equal() +
    NULL

star_pca_fit <- t(log10(star_counts_m + 1)) |> 
  prcomp(scale = TRUE)

star_pca_fit |>
    augment(t(star_counts_m)) |>
    dplyr::rename(sample = .rownames) |>
    mutate(sample = str_remove(sample, '-[0-9]+')) |>
    ggplot(aes(.fittedPC1, .fittedPC2, color = sample)) + 
    geom_point(size = 5)

star_metadata <- data.frame(sample_id = colnames(star_counts_m)) |>
    mutate(tissue = str_sub(sample_id, 1, 3),
           rep = str_sub(sample_id, 6))
rownames(star_metadata) <- star_metadata$sample_id
star_metadata <- select(star_metadata, -sample_id)
star_metadata

all(rownames(star_metadata) == colnames(star_counts_m))

dds <- DESeqDataSetFromMatrix(countData = star_counts_m,
                              colData = star_metadata,
                              design = ~ tissue)
dds <- DESeq(dds)

res <- results(dds)

volcano_data <- as_tibble(res, rownames = 'gene_id')

volcano_plot <- volcano_data |> 
    ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05)) +
    geom_vline(xintercept = 1) +
    geom_vline(xintercept = -1) +
    labs(
        x = 'log2FoldChange',
        y = '-log10(padj)'
    )
  
volcano_plot <- volcano_plot + theme_minimal() + theme(
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)
    )

ggsave("plots/volcano.png", volcano_plot)
# caption on plot: All points are located on the bottom of the plot indicating genes that aren't significantlly differentially expressed.
# caption ct'd: There are points that are located on the right/left - this indicates higher expression of that gene in the intestine vs liver.

vsd <- varianceStabilizingTransformation(dds)

post_deseq2 <- plotPCA(vsd, intgroup = c("tissue"))

ggsave("plots/post_deseq2.png")
# caption on plot: After performing deseq2 analysis, there is now separation of the points. 
# caption ct'd: There are pretty clear clusters of a few points, with only a few being slightly removed from the others



