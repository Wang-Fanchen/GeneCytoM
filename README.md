# GeneCytoM: Multimodal Gene-to-Cytoband Profiling

[![R package version](https://img.shields.io/badge/Version-0.2.0-blue.svg)](https://github.com/Wang-Fanchen/GeneCytoM)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%3E%3D%203.5.0-blue.svg)](https://www.r-project.org/)

GeneCytoM is a comprehensive R package for comparing multimodal (gene + image) models
against gene-only models at both the gene and cytoband levels. It identifies genes whose
importance is enhanced in multimodal settings, maps them to chromosomal cytobands, and
provides rich visualizations to highlight regions of multimodal advantage.

## 🎯 Key Features

- **📍 Gene Location Integration**: Automatically matches genes to chromosomal locations using built-in reference data
- **🔬 Candidate Gene Screening**: Filter genes based on expression patterns and statistical significance
- **📊 Rich Visualizations**: Generate chromosome ideograms, line plots, and lollipop charts
- **📈 Comprehensive Analysis**: Cytoband-level attribution analysis with multiple summary statistics
- **🚀 Easy to Use**: One-line execution of complete analysis pipeline
- **📦 Built-in Data**: Includes comprehensive gene location reference database

## 📦 Installation

### Quick Installation

```r
# Clone or download the package, then run:
source("install_GeneCytoM.R")
```

### Manual Installation

```r
# Install dependencies
install.packages(c("devtools", "roxygen2", "dplyr", "ggplot2", "forcats", "magrittr"))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("karyoploteR", "GenomicRanges"))

# Install the package
devtools::install("path/to/GeneCytoM")
```

## 🚀 Quick Start

### Basic Usage

```r
library(GeneCytoM)

# Run complete analysis with your data
results <- run_cytoband_analysis(
  input_file = "your_gene_data.csv",
  output_dir = "analysis_results"
)

# View results
cat("Found", nrow(results$candidate_genes), "candidate genes\n")
head(results$cytoband_summary)
```

### Input Data Format

#### Required Columns
Your CSV file must contain these columns:

```csv
Gene,Local.Index.G,Local.Index.M,Attribution.G,Attribution.M
TP53,0.523,0.687,-0.234,0.456
BRCA1,0.445,0.512,0.123,0.234
MYC,0.678,0.891,-0.345,0.567
```

#### Optional Columns
You can include P-value columns for filtering:

```csv
Gene,Local.Index.G,Local.Index.M,Attribution.G,Attribution.M,P.Value.G,P.Value.M
TP53,0.523,0.687,-0.234,0.456,0.002,0.001
BRCA1,0.445,0.512,0.123,0.234,0.034,0.023
MYC,0.678,0.891,-0.345,0.567,0.0005,0.0001
```

### Example Analysis

```r
# Run with P-value filtering
results <- run_cytoband_analysis(
  input_file = "gene_data.csv",
  output_dir = "results",
  comparison_filter = "up",
  pvalue_threshold_g = 0.05,     # Filter by P.Value.G
  pvalue_threshold_m = 0.05,     # Filter by P.Value.M
  top_n = 30
)

# Run without P-value filtering
results <- run_cytoband_analysis(
  input_file = "gene_data.csv",
  output_dir = "results",
  comparison_filter = "up"
)

# Access specific results
candidates <- results$candidate_genes
cytoband_summary <- results$cytoband_summary
plots <- results$plots
```

## 📊 Output Files

The analysis generates:

### Data Files (CSV)
- `01_processed_data.csv` - All genes with location information
- `02_candidate_genes.csv` - Filtered candidate genes
- `Candidate_cytoband_attribution_summary.csv` - Cytoband-level attribution
- `Cytoband_summary_5cols.csv` - Comprehensive 5-column summary

### Visualization Files (PDF)
- `chromosome_ideogram_regions.pdf` - Chromosome ideogram with regions
- `chromosome_ideogram_points.pdf` - Chromosome ideogram with points
- `Top20_Cytoband_lineplot.pdf` - Top cytobands line plot
- `Cytoband_Difference_lollipop.pdf` - Difference lollipop chart
- `Cytoband_Ratio_lollipop.pdf` - Ratio lollipop chart

## 🔧 Main Functions

### Core Analysis Functions

| Function | Description |
|----------|-------------|
| `run_cytoband_analysis()` | 🌟 Complete analysis pipeline |
| `process_gene_data()` | Data processing and location matching |
| `screen_candidate_genes()` | Candidate gene filtering |
| `analyze_cytoband_location()` | Cytoband attribution analysis |
| `generate_cytoband_summary()` | Comprehensive summary table |
| `run_example_analysis()` | 🎯 Quick example analysis |

### Visualization Functions

| Function | Description |
|----------|-------------|
| `plot_ideogram()` | Chromosome ideogram visualization |
| `plot_top_cytobands()` | Top cytobands line plot |
| `plot_lollipop_charts()` | Difference and ratio lollipop charts |

## 📚 Documentation

- **Complete User Guide**: See `USER_GUIDE.md` for detailed documentation
- **Installation Instructions**: See `INSTALL_GUIDE.md` for troubleshooting
- **Input Format**: See `INPUT_FORMAT.md` for data requirements

## 🎯 Use Cases

- **Cancer Genomics Research**: Identify chromosomal regions with significant gene expression changes
- **Biomarker Discovery**: Find cytoband regions enriched for candidate genes
- **Pathway Analysis**: Understand chromosomal clustering of regulated genes
- **Publication Figures**: Generate high-quality chromosome visualizations

## ⚠️ Breaking Changes in v0.2.0

Column names have been updated for better clarity:
- `Local.Index.1` → `Local.Index.G`
- `Local.Index.2` → `Local.Index.M`
- `Attribution.1` → `Attribution.G`
- `Attribution.2` → `Attribution.M`
- `P.Value.2` → `P.Value.M`
- New: `P.Value.G` (optional)

**Action required:** Update your input CSV files to use the new column names.

## 📋 Requirements

- R (>= 3.5.0)
- Required packages: dplyr, ggplot2, forcats, magrittr, karyoploteR, GenomicRanges

## 🤝 Contributing

We welcome contributions! Please:

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## 📄 Citation

If you use this package in your research, please cite:

```bibtex
@software{genecytom,
  title={GeneCytoM: Multimodal Gene-to-Cytoband Profiling},
  author={Package Authors},
  year={2025},
  url={https://github.com/Wang-Fanchen/GeneCytoM}
}
```

## 📞 Support

- **Issues**: Report bugs or request features on [GitHub Issues](https://github.com/Wang-Fanchen/GeneCytoM/issues)
- **Documentation**: Check the included guides and R help (`?GeneCytoM`)
- **Email**: your.email@example.com

## 📄 License

This package is licensed under the MIT License. See [LICENSE](LICENSE) for details.

## 🎉 Acknowledgments

- Built with R and Bioconductor
- Gene location data from NCBI Gene database
- Visualization powered by ggplot2 and karyoploteR
- Pipe operations powered by magrittr

---

**Ready to start your analysis?** See the [Quick Start Guide](USER_GUIDE.md) for step-by-step instructions!