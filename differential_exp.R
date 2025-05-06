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
    mutate(color = case_when(
        log2FoldChange > 1  & padj < 0.05 | log2FoldChange < -1 & padj < 0.05 ~ 'deg',
        TRUE ~ 'no deg'
    )) |> 
    ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = color),
        size = 5) +
    geom_hline(yintercept = -log10(0.05)) +
    geom_vline(xintercept = 1) +
    geom_vline(xintercept = -1) +
    labs(
        x = 'log2FoldChange',
        y = '-log10(padj)'
    ) + 
    scale_color_manual(values = c('#EBB6B3', '#334139')) +
    theme_minimal() + 
    theme(
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.position = 'none'
        )
volcano_plot

ggsave('plots/volcano.png', volcano_plot, width = 16, height = 10)
ggsave("plots/volcano.pdf", volcano_plot, width = 16, height = 10)
# caption on plot: All points are located on the bottom of the plot indicating genes that aren't significantlly differentially expressed.
# caption ct'd: There are points that are located on the right/left - this indicates higher expression of that gene in the intestine vs liver.

vsd <- varianceStabilizingTransformation(dds)

post_deseq2 <- plotPCA(vsd, intgroup = c("tissue"), returnData = TRUE) |>
    as_tibble() |>
    mutate(rep = str_extract(name, "[0-9]{2}")) |>
    ggplot() +
    geom_point(aes(x = PC1, y = PC2, color = group, shape = as.factor(rep)), 
                size = 15) +
    scale_shape_manual(values = c(15, 16, 17, 18)) +
    scale_color_manual(values = c('indianred', 'steelblue'),
                        labels = c('Intestine miracidia', 'Liver miracidia')) +
    labs(color = 'Sample', shape = "Replicate") +
    theme_minimal() +
    theme(
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.position = 'bottom'
        )
post_deseq2

ggsave('plots/post_deseq2.png', width = 16, height = 11)
ggsave("plots/post_deseq2.pdf", width = 16, height = 11)
# caption on plot: After performing deseq2 analysis, there is now separation of the points. 
# caption ct'd: There are pretty clear clusters of a few points, with only a few being slightly removed from the others

norm_counts <- read_csv('/data/classes/2025/spring/biol443/course_files/rnaseq_data/deseq_norm_counts.csv')

norm_counts_long <- norm_counts |>
    pivot_longer(-gene_id, names_to = 'sample', values_to = 'norm_count') |>
    mutate(sample = str_remove(sample, "_S[0-9]{2}_L005")) |>
    separate(sample, c('sample', 'replicate'), "-")

selected_genes <- norm_counts_long |>
    filter(gene_id == 'Smp_133770' | gene_id == 'Smp_210300' | gene_id == 'Smp_319430' | gene_id == 'Smp_319450')

DEG_plot <- selected_genes |>
    ggplot(aes(x = sample, y = norm_count, fill = sample)) +
        geom_boxplot() +
        geom_point(aes(shape = replicate), size = 5) +
    facet_wrap(facets = vars(gene_id), scales = "free") +
    scale_fill_manual(values = c('#CD5C5C', '#4682B4'), labels = c('Intestine miracidia', 'Liver miracidia')) +
    scale_shape_manual(values = c(15, 16, 17, 18)) +
    labs(
        title = 'Four Significantly Differentially Expressed Genes', 
        shape = 'Replicate',
        fill = 'Sample'
    ) +
    xlab('Sample Type') +
    ylab('Normalized Count') +
    theme_minimal() + 
    theme(
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 20),
        strip.text = element_text(size = 18),
        axis.text = element_text(size = 18)
)

DEG_plot

ggsave('plots/DEGs_plot.png', DEG_plot, width = 13, height = 9, bg = "white")
ggsave('plots/DEGs_plot.pdf', DEG_plot, width = 13, height = 9, bg = "white")
