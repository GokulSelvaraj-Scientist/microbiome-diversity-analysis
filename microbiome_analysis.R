# ============================================================
# Microbiome Diversity Analysis: Real HMP IBD Study Data
# Author: Gokul Selvaraj
# GitHub: GokulSelvaraj-Scientist
# Description: Alpha and beta diversity analysis of gut microbiome
#              comparing healthy controls vs IBD patients using
#              real data from the Human Microbiome Project (HMP)
# Data: HMP_2019_ibdmdb via curatedMetagenomicData
#       426 healthy controls, 1201 IBD patients
#       585 microbial species, 1627 stool samples
# ============================================================

# --- Install packages if needed ---
# BiocManager::install("curatedMetagenomicData")
# install.packages(c("vegan", "ggplot2", "dplyr", "tidyr", "RColorBrewer"))

library(curatedMetagenomicData)
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

cat("============================================================\n")
cat("Microbiome Diversity Analysis\n")
cat("Human Microbiome Project — IBD Study (HMP_2019_ibdmdb)\n")
cat("Real data: 1627 stool samples, 585 microbial species\n")
cat("============================================================\n\n")

# ============================================================
# STEP 1: LOAD REAL HMP DATA
# ============================================================

cat("=== Loading real HMP microbiome data ===\n")

hmp_tse <- curatedMetagenomicData("HMP_2019_ibdmdb.relative_abundance",
                                   dryrun = FALSE, counts = FALSE)
hmp_se  <- hmp_tse[[1]]

# Extract abundance matrix and metadata
abund   <- as.matrix(assay(hmp_se))       # species x samples
meta_df <- as.data.frame(colData(hmp_se)) # sample metadata

cat("Samples:", ncol(abund), "\n")
cat("Species:", nrow(abund), "\n")
cat("Disease groups:\n")
print(table(meta_df$disease))

# Filter to healthy vs IBD, remove samples with all zeros
keep <- meta_df$disease %in% c("healthy", "IBD")
abund   <- abund[, keep]
meta_df <- meta_df[keep, ]

# Keep species present in at least 10% of samples
prev <- rowMeans(abund > 0)
abund <- abund[prev >= 0.10, ]
cat("\nAfter filtering - Species:", nrow(abund), "| Samples:", ncol(abund), "\n")
cat("Healthy:", sum(meta_df$disease == "healthy"), "\n")
cat("IBD:", sum(meta_df$disease == "IBD"), "\n")

# Transpose for vegan (samples x species)
abund_t <- t(abund)
groups  <- meta_df$disease

# ============================================================
# STEP 2: ALPHA DIVERSITY
# ============================================================

cat("\n=== Computing alpha diversity ===\n")

# Shannon diversity
shannon <- diversity(abund_t, index = "shannon")
# Simpson diversity
simpson <- diversity(abund_t, index = "simpson")
# Observed richness
richness <- specnumber(abund_t)

alpha_df <- data.frame(
  sample_id = rownames(abund_t),
  disease   = groups,
  shannon   = shannon,
  simpson   = simpson,
  richness  = richness
)

# Wilcoxon test
wt_shannon <- wilcox.test(shannon ~ groups, data = alpha_df)
wt_simpson <- wilcox.test(simpson ~ groups, data = alpha_df)
wt_rich    <- wilcox.test(richness ~ groups, data = alpha_df)

cat("Shannon diversity Wilcoxon p =", round(wt_shannon$p.value, 4), "\n")
cat("Simpson diversity Wilcoxon p =", round(wt_simpson$p.value, 4), "\n")
cat("Species richness Wilcoxon p  =", round(wt_rich$p.value, 4), "\n")

# Summary stats
cat("\nAlpha diversity summary:\n")
alpha_summary <- alpha_df %>%
  group_by(disease) %>%
  summarise(
    n           = n(),
    mean_shannon = round(mean(shannon), 3),
    sd_shannon   = round(sd(shannon), 3),
    mean_richness = round(mean(richness), 1)
  )
print(alpha_summary)

# ============================================================
# STEP 3: BETA DIVERSITY
# ============================================================

cat("\n=== Computing beta diversity ===\n")

# Bray-Curtis dissimilarity
bray <- vegdist(abund_t, method = "bray")

# PCoA
pcoa <- cmdscale(bray, k = 4, eig = TRUE)
pcoa_df <- data.frame(
  PC1     = pcoa$points[, 1],
  PC2     = pcoa$points[, 2],
  disease = groups
)

# Variance explained
var_exp <- pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]) * 100

# PERMANOVA
set.seed(42)
perm <- adonis2(bray ~ disease, data = data.frame(disease = groups),
                permutations = 999)
cat("\nPERMANOVA results:\n")
print(perm)

# ANOSIM
ano <- anosim(bray, groups, permutations = 999)
cat("\nANOSIM R statistic:", round(ano$statistic, 4),
    "| p-value:", ano$signif, "\n")

# ============================================================
# STEP 4: TOP DISCRIMINATING SPECIES
# ============================================================

cat("\n=== Identifying discriminating species ===\n")

# Wilcoxon test for each species
healthy_idx <- which(groups == "healthy")
ibd_idx     <- which(groups == "IBD")

species_tests <- apply(abund, 1, function(x) {
  wt <- wilcox.test(x[healthy_idx], x[ibd_idx])
  c(p_value    = wt$p.value,
    mean_healthy = mean(x[healthy_idx]),
    mean_ibd     = mean(x[ibd_idx]))
})

species_df <- as.data.frame(t(species_tests))
species_df$species    <- rownames(species_df)
species_df$p_adjusted <- p.adjust(species_df$p_value, method = "BH")
species_df$log2fc     <- log2((species_df$mean_ibd + 1e-6) /
                               (species_df$mean_healthy + 1e-6))

# Clean species names
species_df$species_clean <- gsub(".*s__", "", species_df$species)
species_df$species_clean <- gsub("_", " ", species_df$species_clean)

top_species <- species_df %>%
  filter(p_adjusted < 0.05) %>%
  arrange(p_adjusted) %>%
  head(20)

cat("Significant species (FDR < 0.05):", sum(species_df$p_adjusted < 0.05), "\n")
cat("\nTop 10 discriminating species:\n")
print(top_species[1:10, c("species_clean", "mean_healthy", "mean_ibd",
                            "log2fc", "p_adjusted")])

# ============================================================
# STEP 5: VISUALISATIONS
# ============================================================

cat("\n=== Generating visualisations ===\n")

GROUP_COLORS <- c("healthy" = "#2A9D8F", "IBD" = "#E76F51")

# --- Plot 1: Alpha diversity ---
alpha_long <- alpha_df %>%
  select(disease, shannon, simpson, richness) %>%
  pivot_longer(cols = c(shannon, simpson, richness),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
    shannon  = "Shannon Diversity",
    simpson  = "Simpson Diversity",
    richness = "Species Richness"))

p_vals <- data.frame(
  metric  = c("Shannon Diversity", "Simpson Diversity", "Species Richness"),
  p_label = c(
    paste0("p = ", formatC(wt_shannon$p.value, format="e", digits=2)),
    paste0("p = ", formatC(wt_simpson$p.value, format="e", digits=2)),
    paste0("p = ", formatC(wt_rich$p.value,    format="e", digits=2))
  )
)

p1 <- ggplot(alpha_long, aes(x = disease, y = value, fill = disease)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA) +
  geom_text(data = p_vals,
            aes(x = 1.5, y = Inf, label = p_label, fill = NULL),
            vjust = 1.5, size = 3.5, color = "#E63946", fontface = "bold") +
  scale_fill_manual(values = GROUP_COLORS) +
  facet_wrap(~metric, scales = "free_y") +
  labs(
    title    = "Alpha Diversity: Healthy vs IBD",
    subtitle = "Real HMP IBD Study Data (n=1,627 stool samples)",
    x        = "Disease Status",
    y        = "Diversity Value",
    fill     = "Group",
    caption  = "Source: Human Microbiome Project (HMP_2019_ibdmdb)"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title  = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text  = element_text(face = "bold", size = 11))

ggsave("01_alpha_diversity.png", p1, width = 12, height = 5, dpi = 300)
cat("Saved: 01_alpha_diversity.png\n")

# --- Plot 2: PCoA beta diversity ---
p2 <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = disease)) +
  geom_point(alpha = 0.5, size = 1.5) +
  stat_ellipse(level = 0.95, linewidth = 1.2) +
  scale_color_manual(values = GROUP_COLORS) +
  labs(
    title    = "Beta Diversity: PCoA of Bray-Curtis Dissimilarity",
    subtitle = paste0("PERMANOVA R² = ",
                      round(perm$R2[1], 3),
                      ", p = ", perm$`Pr(>F)`[1],
                      " | ANOSIM R = ", round(ano$statistic, 3)),
    x        = paste0("PC1 (", round(var_exp[1], 1), "% variance)"),
    y        = paste0("PC2 (", round(var_exp[2], 1), "% variance)"),
    color    = "Group",
    caption  = "Source: Human Microbiome Project (HMP_2019_ibdmdb)"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave("02_beta_diversity_pcoa.png", p2, width = 9, height = 7, dpi = 300)
cat("Saved: 02_beta_diversity_pcoa.png\n")

# --- Plot 3: Top discriminating species ---
top15 <- top_species %>%
  arrange(log2fc) %>%
  mutate(direction = ifelse(log2fc > 0, "Enriched in IBD", "Enriched in Healthy"))

p3 <- ggplot(top15, aes(x = reorder(species_clean, log2fc),
                         y = log2fc, fill = direction)) +
  geom_bar(stat = "identity", alpha = 0.85) +
  coord_flip() +
  scale_fill_manual(values = c("Enriched in IBD"     = "#E76F51",
                                "Enriched in Healthy" = "#2A9D8F")) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.5) +
  labs(
    title    = "Top Discriminating Microbial Species",
    subtitle = "Healthy vs IBD — Wilcoxon test, FDR < 0.05",
    x        = "Species",
    y        = "Log2 Fold Change (IBD / Healthy)",
    fill     = "",
    caption  = "Source: Human Microbiome Project (HMP_2019_ibdmdb)"
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 9),
        legend.position = "bottom")

ggsave("03_discriminating_species.png", p3, width = 10, height = 8, dpi = 300)
cat("Saved: 03_discriminating_species.png\n")

# --- Plot 4: Species abundance heatmap ---
top10_species <- top_species$species[1:10]

# Subsample for heatmap clarity
set.seed(42)
healthy_sub <- sample(healthy_idx, min(50, length(healthy_idx)))
ibd_sub     <- sample(ibd_idx,     min(50, length(ibd_idx)))
sub_idx     <- c(healthy_sub, ibd_sub)

heat_mat <- abund[top10_species, sub_idx]
rownames(heat_mat) <- gsub(".*s__", "", rownames(heat_mat))
rownames(heat_mat) <- gsub("_", " ", rownames(heat_mat))

heat_df <- as.data.frame(heat_mat) %>%
  mutate(species = rownames(.)) %>%
  pivot_longer(-species, names_to = "sample", values_to = "abundance") %>%
  mutate(
    log_abund = log10(abundance + 1e-5),
    disease   = ifelse(sample %in% colnames(abund)[healthy_sub],
                       "Healthy", "IBD")
  )

p4 <- ggplot(heat_df, aes(x = sample, y = species, fill = log_abund)) +
  geom_tile() +
  scale_fill_gradient2(low = "#457B9D", mid = "white", high = "#E63946",
                        midpoint = median(heat_df$log_abund),
                        name = "Log10\nAbundance") +
  facet_grid(. ~ disease, scales = "free_x", space = "free") +
  labs(
    title    = "Top 10 Discriminating Species — Abundance Heatmap",
    subtitle = "Subsampled (50 healthy + 50 IBD patients)",
    x        = "Sample",
    y        = "Species",
    caption  = "Source: Human Microbiome Project (HMP_2019_ibdmdb)"
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title   = element_text(face = "bold"),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text   = element_text(face = "bold", size = 11))

ggsave("04_species_heatmap.png", p4, width = 11, height = 7, dpi = 300)
cat("Saved: 04_species_heatmap.png\n")

# --- Plot 5: Phylum-level composition ---
# Extract phylum from species names
phylum_df <- data.frame(
  species = rownames(abund),
  phylum  = gsub(".*p__([^|]+).*", "\\1", rownames(abund))
)

# Sum abundances by phylum
phylum_abund <- abund %>%
  as.data.frame() %>%
  mutate(phylum = phylum_df$phylum) %>%
  group_by(phylum) %>%
  summarise(across(everything(), sum)) %>%
  as.data.frame()

rownames(phylum_abund) <- phylum_abund$phylum
phylum_mat <- as.matrix(phylum_abund[, -1])

# Get top 6 phyla
top_phyla <- names(sort(rowMeans(phylum_mat), decreasing = TRUE))[1:6]
phylum_plot <- phylum_mat[top_phyla, ] %>%
  as.data.frame() %>%
  mutate(phylum = rownames(.)) %>%
  pivot_longer(-phylum, names_to = "sample", values_to = "abundance") %>%
  mutate(disease = ifelse(sample %in% colnames(abund)[healthy_idx],
                          "Healthy", "IBD"))

phylum_summary <- phylum_plot %>%
  group_by(disease, phylum) %>%
  summarise(mean_abund = mean(abundance), .groups = "drop")

p5 <- ggplot(phylum_summary, aes(x = disease, y = mean_abund, fill = phylum)) +
  geom_bar(stat = "identity", position = "fill", alpha = 0.85) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title    = "Phylum-Level Microbiome Composition",
    subtitle = "Healthy vs IBD — mean relative abundance",
    x        = "Disease Status",
    y        = "Relative Abundance",
    fill     = "Phylum",
    caption  = "Source: Human Microbiome Project (HMP_2019_ibdmdb)"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave("05_phylum_composition.png", p5, width = 8, height = 6, dpi = 300)
cat("Saved: 05_phylum_composition.png\n")

# ============================================================
# SAVE RESULTS
# ============================================================

write.csv(alpha_summary, "alpha_diversity_summary.csv", row.names = FALSE)
write.csv(top_species[, c("species_clean","mean_healthy","mean_ibd",
                           "log2fc","p_adjusted")],
          "discriminating_species.csv", row.names = FALSE)

perm_results <- data.frame(
  Test      = c("PERMANOVA", "ANOSIM"),
  Statistic = c(round(perm$R2[1], 4), round(ano$statistic, 4)),
  P_value   = c(perm$`Pr(>F)`[1], ano$signif)
)
write.csv(perm_results, "beta_diversity_stats.csv", row.names = FALSE)

cat("\nAll results saved.\n")
cat("\n=== Analysis Complete ===\n")
cat("Real HMP samples analysed:", ncol(abund), "\n")
cat("Healthy:", sum(groups == "healthy"),
    "| IBD:", sum(groups == "IBD"), "\n")
cat("Significant species:", sum(species_df$p_adjusted < 0.05), "\n")
