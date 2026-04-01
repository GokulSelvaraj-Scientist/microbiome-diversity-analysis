# Microbiome Diversity Analysis: Real HMP IBD Study Data

## Overview
Alpha and beta diversity analysis of real gut microbiome data from the Human Microbiome Project (HMP) IBD study, comparing 426 healthy controls against 1,201 IBD patients across 1,627 stool samples and 585 microbial species. Demonstrates quantitative microbiome analysis methods directly applicable to microbiome drug development endpoints and clinical trial biomarker strategy.

## Data Source
- **Database:** Human Microbiome Project (HMP)
- **Study:** HMP_2019_ibdmdb (Inflammatory Bowel Disease Multi-site Database)
- **Access:** Public — via curatedMetagenomicData Bioconductor package
- **Samples:** 1,627 real stool samples
- **Groups:** 426 healthy controls vs 1,201 IBD patients
- **Species:** 585 microbial species (after prevalence filtering)

## Why This Matters
The gut microbiome is increasingly recognised as a key factor in disease and drug response. Microbiome analysis skills are directly relevant to:

- **IBD drug development** — biologics targeting gut inflammation (vedolizumab, ustekinumab) require microbiome biomarker strategies
- **Microbiome therapeutics** — FMT, live biotherapeutics, and microbiome modulators are an active drug development area
- **Clinical trial endpoints** — microbiome diversity metrics are used as secondary endpoints in IBD, metabolic disease, and oncology trials
- **Pharmacogenomics** — gut microbiome composition affects drug metabolism and response

## Analysis

### Alpha Diversity
- Shannon diversity, Simpson diversity, and species richness
- Wilcoxon rank-sum tests comparing healthy vs IBD
- All three metrics significantly reduced in IBD (p < 1e-8)

### Beta Diversity
- Bray-Curtis dissimilarity matrix
- PCoA ordination with 95% confidence ellipses
- PERMANOVA (999 permutations): R² = 0.012, p = 0.001
- ANOSIM: confirms significant community composition differences

### Discriminating Species
- Wilcoxon test for each of 585 species
- FDR correction (Benjamini-Hochberg)
- Most discriminating species depleted in IBD: Ruminococcus bicirculans, Roseburia hominis, Gemmiger formicilis, Akkermansia muciniphila
- Species enriched in IBD: Ruminococcus torques

### Phylum-Level Composition
- Firmicutes and Bacteroidetes dominate in both groups
- Subtle shifts in phylum composition between healthy and IBD

## Key Findings
- IBD is associated with significant reduction in gut microbiome diversity across all three alpha diversity metrics
- PERMANOVA confirms significant differences in overall community composition (p = 0.001)
- Key beneficial bacteria (Roseburia, Ruminococcus, Akkermansia) are depleted in IBD — consistent with published literature
- Ruminococcus torques is enriched in IBD — consistent with its role as a pro-inflammatory species

## Outputs
| File | Description |
|---|---|
| `01_alpha_diversity.png` | Shannon, Simpson, species richness violin plots |
| `02_beta_diversity_pcoa.png` | PCoA with PERMANOVA and ANOSIM statistics |
| `03_discriminating_species.png` | Top 20 discriminating species log2 fold change |
| `04_species_heatmap.png` | Top 10 species abundance heatmap |
| `05_phylum_composition.png` | Phylum-level composition stacked barplot |
| `alpha_diversity_summary.csv` | Alpha diversity summary statistics |
| `discriminating_species.csv` | Significant species with fold changes |
| `beta_diversity_stats.csv` | PERMANOVA and ANOSIM results |

## How to Run

### Step 1 — Install R packages
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("curatedMetagenomicData")
install.packages(c("vegan", "ggplot2", "dplyr", "tidyr", "RColorBrewer"))
```

### Step 2 — Run analysis
```r
setwd("path/to/microbiome-diversity-analysis")
source("microbiome_analysis.R")
```

Note: The first run downloads ~500MB of HMP data via ExperimentHub. Subsequent runs use the cached version.

## Technical Stack
- **Data access:** curatedMetagenomicData (Bioconductor)
- **Diversity analysis:** vegan (diversity, vegdist, adonis2, anosim)
- **Statistics:** Wilcoxon rank-sum test, FDR correction, PERMANOVA, ANOSIM
- **Visualisation:** ggplot2, RColorBrewer
- **Data:** Real HMP IBD study data (1,627 stool samples)

## Author
**Gokul Selvaraj, PhD**
GitHub: [GokulSelvaraj-Scientist](https://github.com/GokulSelvaraj-Scientist)
