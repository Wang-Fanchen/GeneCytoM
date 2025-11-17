# GeneCytoM: Multimodal Gene-to-Cytoband Profiling

[![R package version](https://img.shields.io/badge/Version-0.2.0-blue.svg)](https://github.com/Wang-Fanchen/GeneCytoM)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%3E%3D%203.5.0-blue.svg)](https://www.r-project.org/)

> **A comprehensive R package for comparing multimodal AI models (gene + image) against gene-only models to identify gene/chromosomal-regions with enhanced multimodal performance.**

---

## 📖 Overview

**GeneCytoM** (Gene-Cytoband-Multimodal) is designed for researchers working with **multimodal machine learning models** that integrate gene expression and imaging data (e.g., histopathology images). 

The package helps answer a critical question:

> *"Which genes and chromosomal regions show improved prediction importance when we add imaging data to gene-only models?"*

### What GeneCytoM Does

1. **Identifies multimodal-advantaged genes**: Genes where the multimodal model (gene + image) assigns higher importance than the gene-only model
2. **Maps genes to cytobands**: Automatically locates genes on chromosomal cytobands using built-in reference data
3. **Performs cytoband-level analysis**: Aggregates gene-level signals to identify chromosomal regions enriched for multimodal advantages
4. **Generates publication-ready visualizations**: Creates chromosome ideograms, comparative line plots, and lollipop charts

### Scientific Context

This package is particularly useful for:
- **Cancer genomics**: Comparing gene expression + tumor histology models vs gene-only models
- **Biomarker discovery**: Finding chromosomal hotspots where imaging adds predictive value
- **Multimodal AI evaluation**: Understanding which genomic regions benefit most from multimodal integration
- **Translational research**: Identifying cytobands for targeted experimental validation

---

## 🎯 Key Features

| Feature | Description |
|---------|-------------|
| 🧬 **Gene-Level Comparison** | Compare gene importance between multimodal (M) and gene-only (G) models |
| 📍 **Automatic Cytoband Mapping** | Built-in database of 206,757 human genes with chromosomal locations |
| 🔍 **Flexible Filtering** | Optional P-value thresholds for both model types |
| 📊 **Multiple Visualization Modes** | 7+ publication-ready plots with 4 different sorting strategies |
| ⚡ **One-Line Analysis** | Complete pipeline from raw data to visualizations |
| 🎨 **Comprehensive Reports** | CSV summaries + PDF figures for immediate interpretation |

---

## 🚀 Quick Start

### Installation

```r
# Method 1: Install from GitHub
devtools::install_github("Wang-Fanchen/GeneCytoM")

# Method 2: Install from local source
devtools::install("path/to/GeneCytoM")
```

**Required dependencies:**
```r
# CRAN packages
install.packages(c("dplyr", "ggplot2", "forcats", "magrittr"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("karyoploteR", "GenomicRanges", "regioneR"))
```

### Basic Usage

```r
library(GeneCytoM)

# Run complete analysis (one line!)
results <- run_cytoband_analysis(
  input_file = "your_gene_data.csv",
  output_dir = "results"
)

# Run complete analysis (detailed version)
run_cytoband_analysis(
  input_file = "your_gene_data.csv",
  output_dir = "results"
  comparison_filter = "all",
  pvalue_threshold_g = NULL,
  pvalue_threshold_m = 0.05,
  top_n = 20,
  save_intermediate = TRUE,
  create_plots = TRUE,
  verbose = TRUE
)


# View key results
cat("Multimodal-advantaged genes:", nrow(results$candidate_genes), "\n")
cat("Enriched cytobands:", nrow(results$cytoband_summary), "\n")
head(results$cytoband_summary)
```

---

## 📋 Input Data Format

### Understanding the Column Names

GeneCytoM uses a clear naming convention:
- **G suffix** = Gene-only model metrics
- **M suffix** = Multimodal model metrics (gene + image)

### Required Columns

Your input CSV must include these 5 columns:

| Column | Description | Example |
|--------|-------------|---------|
| `Gene` | Gene symbol (HGNC) | TP53, BRCA1, MYC |
| `Local.Index.G` | Gene importance rank in gene-only model | 0.523 |
| `Local.Index.M` | Gene importance rank in multimodal model | 0.687 |
| `Attribution.G` | Attribution score in gene-only model | -0.234 |
| `Attribution.M` | Attribution score in multimodal model | 0.456 |

* Note that contribution values cannot be negative. The negative sign here indicates high/low risk as determined by the model, not the actual magnitude of the contribution. Thus, all contribution values are automatically converted to their absolute values when included in calculations.

**Minimal example (minimal_input.csv):**
```csv
Gene,Local.Index.G,Local.Index.M,Attribution.G,Attribution.M
TP53,0.523,0.687,-0.234,0.456
BRCA1,0.445,0.512,0.123,0.234
MYC,0.678,0.891,-0.345,0.567
EGFR,0.234,0.456,0.089,0.178
```

### Optional Columns (Recommended)

Add P-value columns for statistical filtering:

| Column | Description |
|--------|-------------|
| `P.Value.G` | Statistical significance in gene-only model |
| `P.Value.M` | Statistical significance in multimodal model |

**Full example (with_pvalues.csv):**
```csv
Gene,Local.Index.G,Local.Index.M,Attribution.G,Attribution.M,P.Value.G,P.Value.M
TP53,0.523,0.687,-0.234,0.456,0.002,0.001
BRCA1,0.445,0.512,0.123,0.234,0.034,0.023
MYC,0.678,0.891,-0.345,0.567,0.0005,0.0001
```

---

## 💡 Usage Examples

### Example 1: Basic Analysis (No P-value Filtering)

```r
library(GeneCytoM)

# Identify genes where multimodal > gene-only
results <- run_cytoband_analysis(
  input_file = "mydata.csv",
  output_dir = "multimodal_analysis",
  comparison_filter = "up"  # "up" means M > G in Local_index
)
```

### Example 2: With P-value Filtering

```r
# Strict filtering: significant in both models
results <- run_cytoband_analysis(
  input_file = "mydata.csv",
  output_dir = "strict_analysis",
  comparison_filter = "up",
  pvalue_threshold_g = 0.05,  # Filter gene-only model
  pvalue_threshold_m = 0.01,  # Stricter filter for multimodal
  top_n = 30                  # Show top 30 cytobands
)
```

### Example 3: Exploring Different Comparisons

```r
# Find genes where multimodal is WORSE than gene-only
downregulated <- run_cytoband_analysis(
  input_file = "mydata.csv",
  output_dir = "downregulated",
  comparison_filter = "down"
)

# Analyze ALL genes regardless of direction
all_genes <- run_cytoband_analysis(
  input_file = "mydata.csv",
  output_dir = "all_analysis",
  comparison_filter = "all"
)
```

### Example 4: Accessing Results

```r
results <- run_cytoband_analysis("mydata.csv", "output")

# Candidate genes (multimodal-advantaged)
candidates <- results$candidate_genes
write.csv(candidates, "my_candidate_genes.csv")

# Cytoband summary with statistics
cytoband_stats <- results$cytoband_summary
top_cytobands <- head(cytoband_stats, 10)

# Find cytobands with highest multimodal advantage
advantaged <- cytoband_stats[cytoband_stats$Difference > 0, ]
```

---

## 📊 Output Files

GeneCytoM generates a comprehensive set of outputs organized in your specified directory:

### 📁 Directory Structure

```
output_dir/
├── 01_processed_data.csv                    # All genes with location info
├── 02_candidate_genes.csv                   # Filtered multimodal-advantaged genes
├── Candidate_cytoband_attribution_summary.csv  # Cytoband-level statistics
├── Cytoband_summary_5cols.csv               # Comprehensive 5-column summary
└── plots/
    ├── chromosome_ideogram_regions.pdf      # Genes shown as chromosomal regions
    ├── chromosome_ideogram_points.pdf       # Genes shown as points
    ├── Top20_Cytoband_lineplot.pdf          # Ranked by multimodal attribution
    ├── Top20_Cytoband_by_Difference.pdf     # Ranked by (M - G) difference ⭐
    ├── Top20_Cytoband_by_Ratio.pdf          # Ranked by (M / G) ratio ⭐
    ├── Cytoband_Difference_lollipop.pdf     # Difference visualization
    └── Cytoband_Ratio_lollipop.pdf          # Ratio visualization
```

### 📈 Understanding the Output Files

#### 1. **Candidate Genes CSV** (`02_candidate_genes.csv`)
Lists all genes where multimodal shows advantage:
- Gene symbols and chromosomal locations
- Attribution scores for both models
- Statistical significance (if P-values provided)
- Comparison status (up/down/mid)

#### 2. **Cytoband Summary CSV** (`Candidate_cytoband_attribution_summary.csv`)
Cytoband-level aggregated statistics:
- `Attribution.G_sum`: Total gene-only attribution
- `Attribution.M_sum`: Total multimodal attribution
- `Difference`: M - G (positive = multimodal advantage)
- `Ratio`: M / G (>1 = multimodal advantage)
- `n_genes`: Number of genes in the cytoband

#### 3. **Visualization PDFs**
Seven publication-ready plots showing:
- **Ideograms**: Chromosome-wide distribution of candidate genes
- **Line plots** (3 versions): Different ranking strategies to highlight different aspects
- **Lollipop charts**: Easy-to-interpret difference and ratio visualizations

---

## 🔧 Function Reference

### Core Pipeline Function

```r
run_cytoband_analysis(
  input_file,                    # Path to input CSV
  output_dir = "results",        # Output directory
  comparison_filter = "up",      # "up", "down", "mid", or "all"
  pvalue_threshold_g = NULL,     # P-value cutoff for gene-only model
  pvalue_threshold_m = NULL,     # P-value cutoff for multimodal model
  top_n = 20,                    # Number of top cytobands to visualize
  create_plots = TRUE,           # Generate visualizations?
  save_intermediate = TRUE       # Save intermediate results?
)
```

### Individual Component Functions

| Function | Purpose |
|----------|---------|
| `process_gene_data()` | Match genes to chromosomal locations |
| `screen_candidate_genes()` | Filter genes by comparison type and P-values |
| `analyze_cytoband_location()` | Aggregate to cytoband level |
| `generate_cytoband_summary()` | Create comprehensive 5-column summary |
| `plot_ideogram()` | Generate chromosome ideograms |
| `plot_top_cytobands()` | Create line plots (by multimodal attribution) |
| `plot_top_cytobands_by_difference()` | Create line plots (by M-G difference) |
| `plot_top_cytobands_by_ratio()` | Create line plots (by M/G ratio) |
| `plot_lollipop_charts()` | Generate lollipop visualizations |

---

## 📚 Real-World Use Case Example

### Scenario: Breast Cancer Multimodal Prediction

**Research question**: Does adding histopathology images to gene expression improve breast cancer subtype prediction? If so, which genomic regions benefit most?

**Workflow**:

1. **Train two models**:
   - Model G: Gene expression only
   - Model M: Gene expression + tumor histopathology images

2. **Extract feature importance**:
   ```r
   # After model training, extract gene importance scores
   # Format data with columns: Gene, Local.Index.G, Local.Index.M, etc.
   write.csv(gene_importance, "BRCA_comparison.csv")
   ```

3. **Run GeneCytoM**:
   ```r
   library(GeneCytoM)
   
   results <- run_cytoband_analysis(
     input_file = "BRCA_comparison.csv",
     output_dir = "BRCA_multimodal_analysis",
     comparison_filter = "up",
     pvalue_threshold_m = 0.05,
     top_n = 20
   )
   ```

4. **Interpret results**:
   ```r
   # How many genes benefit from images?
   cat("Multimodal-advantaged genes:", nrow(results$candidate_genes))
   
   # Which cytobands are enriched?
   top_cytobands <- head(results$cytoband_summary, 10)
   print(top_cytobands)
   
   # Example output:
   #   Cytoband  Attribution.G  Attribution.M  Difference  Ratio  n_genes
   #   17q12     105.3          143.2          37.9        1.36   7
   #   8q21.3    191.4          242.0          51.0        1.27   17
   ```

5. **Biological interpretation**:
   - Cytoband 17q12 (contains ERBB2/HER2) shows 36% increase with imaging
   - Suggests histology captures HER2 protein expression patterns
   - Validates known biology: HER2+ tumors have distinct morphology

---

## 🎨 Visualization Gallery

### Four Sorting Strategies

GeneCytoM provides **4 different ways** to rank and visualize cytobands:

| Sorting Method | When to Use | Highlights |
|----------------|-------------|------------|
| **By Multimodal Attribution** | Overall importance | Cytobands most important in multimodal model |
| **By Difference (M - G)** | Absolute improvement | Regions with largest increase from adding images |
| **By Ratio (M / G)** | Relative improvement | Regions with proportionally greatest benefit |
| **Lollipop Charts** | Quick comparison | Visual difference and ratio in one view |

### Example Visualization Outputs

**Chromosome Ideogram**: Shows spatial distribution of candidate genes
- Quickly identify chromosomes with many multimodal-advantaged genes
- See clustering patterns (e.g., multiple genes on chr17q)

**Line Plots**: Compare gene-only (blue) vs multimodal (red) attribution
- Red line above blue = multimodal advantage
- Larger gap = greater improvement from imaging

**Lollipop Charts**: Emphasize differences and ratios
- Longer sticks = larger differences
- Points above reference line (ratio=1) = multimodal benefit

---

## 📊 Interpreting Results

### Understanding the Metrics

**Difference (M - G)**:
- **Positive**: Multimodal model assigns more importance (good for imaging utility)
- **Negative**: Gene-only model performs better (imaging may add noise)
- **Large positive values**: Strong candidates for imaging-based biomarkers

**Ratio (M / G)**:
- **> 1**: Multimodal improvement
- **< 1**: Gene-only is better
- **>> 1**: Dramatic multimodal advantage (e.g., 2.0 = 100% increase)

### Statistical Significance

When P-values are provided:
- Genes must be significant in the specified model(s)
- Helps control false positives
- Recommended thresholds: 0.01-0.05 depending on study goals

### Biological Context

Consider these questions:
1. **Are enriched cytobands in known cancer-related regions?**
   - e.g., 17q12 (HER2), 17p13 (TP53), 8q24 (MYC)

2. **Do multimodal-advantaged genes have visual phenotypes?**
   - e.g., proteins visible in H&E staining

3. **Are results consistent with clinical knowledge?**
   - e.g., ER+ tumors have distinct morphology → ER genes benefit from imaging

---

## ⚠️ Important Notes

**About attribution score**: Note that contribution values cannot be negative. The negative sign here indicates high/low risk as determined by the model, not the actual magnitude of the contribution. Thus, all contribution values are automatically converted to their absolute values when included in calculations.

**Action required**: Update your CSV files to use new column names.

### Data Requirements

- **Minimum**: 5 required columns (Gene, Local.Index.G, Local.Index.M, Attribution.G, Attribution.M)
- **Recommended**: Include P-value columns for robust filtering
- **Gene names**: Use official HGNC symbols for best matching (TP53, not p53)
- **CSV format**: UTF-8 encoding, standard comma delimiter

---

## 🤝 Contributing

We welcome contributions! Ways to contribute:

1. **Report bugs**: [GitHub Issues](https://github.com/Wang-Fanchen/GeneCytoM/issues)
2. **Request features**: Describe your use case in an issue
3. **Submit pull requests**: 
   - Fork the repository
   - Create a feature branch
   - Make your changes
   - Submit a PR with clear description

### Development Setup

```r
# Clone repository
git clone https://github.com/Wang-Fanchen/GeneCytoM.git
cd GeneCytoM

# Install development dependencies
devtools::install_dev_deps()

# Run checks
devtools::check()

# Run tests (when available)
devtools::test()
```

---

## 📄 Citation

If you use GeneCytoM in your research, please cite:

```bibtex
@software{genecytom2025,
  title = {GeneCytoM: Multimodal Gene-to-Cytoband Profiling},
  author = {Wang, Fanchen},
  year = {2025},
  url = {https://github.com/Wang-Fanchen/GeneCytoM},
  version = {0.2.0},
  note = {R package for comparing multimodal and gene-only AI models}
}
```

---

## 📞 Support & Contact

- **GitHub Issues**: [Report bugs or request features](https://github.com/Wang-Fanchen/GeneCytoM/issues)
- **Documentation**: Run `?GeneCytoM` or `?run_cytoband_analysis` in R
- **Email**: wang.fanchen@outlook.com

---

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

**Key points**:
- ✅ Free for academic and commercial use
- ✅ Modify and distribute
- ✅ Private use
- ⚠️ No warranty provided
- ℹ️ Attribution required

---

## 🙏 Acknowledgments

**Data sources**:
- Gene location data: [NCBI Gene Database](https://www.ncbi.nlm.nih.gov/gene)
- Genome assembly: GRCh38/hg38

**Built with**:
- [R](https://www.r-project.org/) - Statistical computing language
- [dplyr](https://dplyr.tidyverse.org/) - Data manipulation
- [ggplot2](https://ggplot2.tidyverse.org/) - Data visualization
- [karyoploteR](https://bioconductor.org/packages/karyoploteR/) - Chromosome ideograms
- [Bioconductor](https://www.bioconductor.org/) - Genomics infrastructure

**Inspired by**:
- Multimodal AI research in cancer genomics
- Integration of imaging and molecular data
- Translational bioinformatics needs

---

## 🗺️ Roadmap

### Planned Features (v0.3.0)

- [ ] Batch analysis mode for multiple comparisons
- [ ] Gene set enrichment analysis integration

### Under Consideration

- [ ] Support for non-human organisms
- [ ] Integration with pathway databases (KEGG, Reactome)
- [ ] Statistical significance testing for cytobands
- [ ] Web-based interface (Shiny app)

**Have suggestions?** [Open an issue](https://github.com/Wang-Fanchen/GeneCytoM/issues) with your ideas!

---

<div align="center">

**⭐ If you find GeneCytoM useful, please star the repository! ⭐**

Made with ❤️ for the multimodal AI and cancer genomics community

</div>