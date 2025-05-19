# 🧪 Metadoon

***Metadoon*** is a user-friendly graphical interface and pipeline designed for processing and analyzing amplicon-based metagenomic data using tools like VSEARCH and R (with Phyloseq). It automates the workflow from FASTQ preprocessing to statistical visualization in R.

---

## *📦 Major Dependencies*

| *Dependency*                                              | *Version (Suggested)* | *Description*                                  |
| --------------------------------------------------------- | --------------------- | ---------------------------------------------- |
| *[Python](https://www.python.org/downloads/)*             | *3.12+*               | *Main interface (Tkinter GUI, logic control)*  |
| *[R](https://cran.r-project.org/)*                        | *4.4.3*               | *Statistical analysis and plotting*            |
| *[VSEARCH](https://github.com/torognes/vsearch/releases)* | *≥ 2.21.1*            | *FASTQ processing (dereplication, clustering)* |

---

## *🐍 Python Dependencies (minor)*

*These packages are included in the Conda environment:*

| *Package*  | *Purpose*                       |
| ---------- | ------------------------------- |
| *`tk`*     | *GUI interface with Tkinter*    |
| *`pillow`* | *Icon/image handling in Python* |

> \*Standard libraries used: **`os`**, **`sys`**, **`json`**, **`threading`**, \**`subprocess`*

---

## *📊 R Dependencies (minor)*

All packages listed below are automatically installed by the Conda environment.

### *📦 CRAN Packages*

* tidyverse, reshape2, igraph, foreach, lme4
* ggplot2, ggpubr, cowplot, dplyr, pheatmap, viridis
* ape, rprojroot

### *🧪 Bioconductor Packages*

* phyloseq, DESeq2, scater

### *🔧 GitHub Packages (installed post-env)*

```r
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("vlubitch/pairwiseAdonis")
devtools::install_github("microbiome/microbiome")
```

---

## *🛠️ System-Level Dependencies*

*These are installed via Conda or available on Unix-based systems:*

| *Tool/Library*                                           | *Description*                       |
| -------------------------------------------------------- | ----------------------------------- |
| `bash`, `wget`                                           | *Script automation and downloading* |
| `conda`, `mamba`                                         | *Environment management*            |
| `Rscript`                                                | *Execute R scripts via CLI*         |
| `libcurl`, `libxml2`                                     | *R package compilation*             |
| `openssl`, `zlib`, `gcc`, `make`, `libuv`, `gmp`, `mpfr` | *System/compiler libraries*         |

> *On macOS, **************XQuartz************** may be required for full R graphical support.*

---

## *🔗 Download Links*

* [Python](https://www.python.org/downloads/)
* [R](https://cran.r-project.org/)
* [VSEARCH](https://github.com/torognes/vsearch/releases)
* [Conda (recommended)](https://docs.conda.io/en/latest/)

---

## *🚀 Installation & Usage*

1. *Clone this repository:*

   ```bash
   git clone https://github.com/rdo-adan/metadoon.git
   cd metadoon
   ```

2. *Create the Conda environment:*

   ```bash
   bash setup.sh
   ```

3. Activate the conda enviroment:

   ```
   conda activate metadoon
   ```

4. *Run the GUI:*

   ```bash
   python metadoon.py
   ```
5. *Using Interface:*

<p align="justify">Upon launching the application, users are greeted with a clean dashboard divided into logical sections for data input, parameter configuration, and pipeline execution. Users begin by selecting the input FASTQ files(Upper Menu > Files > Upload Files) (paired-end) through a file browser. The interface provides a configuration screen (Configura Execution), and it shows fields for setting key parameters such as minimum read length, maximum expected errors (maxEE), clustering identity thresholds, and selection of chimera detection and taxonomic reference databases. Each field includes tooltips to assist users in understanding the function of each parameter.</p>

<p align="justify">Once all inputs and configurations are set, the user can initiate the pipeline by clicking the “Run Pipeline” button. The interface provides real-time status updates and alerts in a console panel, displaying each processing stage as it is executed—from merging and filtering reads to OTU clustering and diversity analysis. Upon completion, output folders are automatically created and populated with organized results, including quality plots, taxonomic bar charts, core microbiome heatmaps, and diversity metrics. These features make Metadoon accessible to researchers with limited command-line experience while preserving analytical transparency and control.</p>

---

## *📁 Project Structure*

```
Metadoon/
│
├── metadoon.py             # Main GUI
├── Analise.R               # R script for data analysis
├── setup.sh                # Environment setup script
├── metadoon_env.yaml       # Conda environment file
├── *.png, *.ico, *.icns    # Icons and GUI assets
└── Output/, Metadata/, OTUs/, Taxonomy/, Tree File/
```

---

## *📬 Contact*

*For issues or questions, please open an issue or contact the maintainer.*

[rdo.adan@gmail.com](mailto:rdo.adan@gmail.com)
