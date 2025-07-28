library(tidyverse)
library(DESeq2)
library(broom)
library(ggrepel)

counts_df <- read_tsv("/data/users/willetse0745/Sm_Mira_IvT/pipeline/counts/star/counts.tsv", comment = "#") |>
             mutate(across(where(is.numeric), as.integer))

counts_summary <- counts_df |>
    select(Geneid, contains('bam')) |>
    rename_with(~str_remove(., "dedup/star/.*.bam:"), everything()) |>
    rowwise() |>
    mutate(total_counts = sum(c_across(where(is.numeric)), na.rm = T)) |>
    filter(total_counts >= 10)

sample_summary <- counts_df |>
    select(Geneid, contains('bam')) |>
    rename_with(~str_remove(., "dedup/star/.*.bam:"), everything()) |>
    pivot_longer(-Geneid, names_to = 'sample', values_to = 'count') |>
    filter(count > 0) |>
    group_by(Geneid) |>
    tally() |>
    filter(n <= 3)

genes_to_remove = sample_summary$Geneid

counts_filt <- counts_summary |>
    filter(!Geneid %in% genes_to_remove) |>
    arrange(Geneid) |>
    select(-total_counts)

counts_m <- counts_filt |>
    select(-Geneid) |>
    as.matrix()
rownames(counts_m) <- counts_filt$Geneid

metadata <- data.frame(sample_id = colnames(counts_m)) |>
    mutate(tissue = str_sub(sample_id, 1, 3),
           rep = str_sub(sample_id, 6, 6))
rownames(metadata) <- metadata$sample_id
metadata <- select(metadata, -sample_id)

dds <- DESeqDataSetFromMatrix(countData = counts_m,
                              colData = metadata,
                              design = ~ tissue)
dds <- DESeq(dds)

res <- results(dds)

vsd <- vst(dds)

pca_fit <- t(assay(vsd)) |> 
  prcomp(scale = TRUE)

(pca_plot <- pca_fit |>
  augment(t(assay(vsd))) |>
  dplyr::rename(sample = .rownames) |>
  separate(sample, into = c('tissue', 'rep'), sep = '-') |>
    mutate(tissue = case_when(
    tissue == 'Int' ~ 'Intestine',
    tissue == 'Liv' ~ 'Liver'
  )) |>
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, color = tissue, shape = tissue)) + 
  geom_point(size = 4) +
  scale_color_manual(values = c('indianred', 'steelblue')) +
  labs(x = "PC1 (29% of variance)", y = "PC2 (23% of variance)", color = "Tissue source", shape = "Tissue source") +
  theme_minimal() +
  NULL
)

ggsave('plots/pca.pdf', pca_plot, width = 4, height = 4)
ggsave('plots/pca.png', pca_plot, width = 4, height = 4)

volcano_data <- as_tibble(res, rownames = 'gene_id')

degs <- volcano_data |>
    filter(log2FoldChange > 2 | log2FoldChange < -2,
            -log10(padj) > -log10(0.05))

volcano_plot <- volcano_data |> 
    mutate(color = case_when(
        log2FoldChange > 1  & padj < 0.05 | log2FoldChange < -1 & padj < 0.05 ~ 'deg',
        TRUE ~ 'not deg'
    )) |> 
    ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = color),
        size = 5) +
    geom_label_repel(data = degs,
                     aes(x = log2FoldChange, y = -log10(padj), label = gene_id),
                     max.overlaps = 50, size = 1) +
    geom_hline(yintercept = -log10(0.05)) +
    geom_vline(xintercept = 1) +
    geom_vline(xintercept = -1) + 
    scale_color_manual(values = c('#EBB6B3', '#334139')) +
    theme_minimal()
volcano_plot

ggsave('plots/volcano.png', volcano_plot, bg = 'white', width = 6, height = 6)
ggsave("plots/volcano.pdf", volcano_plot, bg = 'white', width = 6, height = 6)
# caption on plot: All points are located on the bottom of the plot indicating genes that aren't significantlly differentially expressed.
# caption ct'd: There are points that are located on the right/left - this indicates higher expression of that gene in the intestine vs liver.

norm_counts <- counts(dds, normalized = TRUE) |>
    as_tibble(rownames = 'gene_id')
write_csv(norm_counts, '/data/users/willetse0745/Sm_Mira_IvT/pipeline/deseq_results/deseq_norm_counts.csv')



## from given normalized counts
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
