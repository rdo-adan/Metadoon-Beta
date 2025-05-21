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

> \*Standard libraries used: **`os`**, **`sys`**, **`json`**, **`threading`**, **`subprocess`**

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

> *On macOS, **XQuartz** may be required for full R graphical support.*

---

## *🔗 Download Links*

* [Python](https://www.python.org/downloads/)
* [R](https://cran.r-project.org/)
* [VSEARCH](https://github.com/torognes/vsearch/releases)
* [Conda (recommended)](https://docs.conda.io/en/latest/)

---

## *🚀 Installation & Usage*

> ⚠️ *Note: All required folders such as `Output/`, `Metadata/`, `OTUs/`, `Taxonomy/`, and `Tree File/` are automatically created during the pipeline execution if they do not exist.*

1. *Clone this repository:*

   ```bash
   git clone https://github.com/rdo-adan/metadoon.git
   cd metadoon
   ```

2. *Create the Conda environment:*

   ```bash
   bash setup.sh
   ```

3. *Activate the Conda environment:*

   ```bash
   conda activate metadoon
   ```

4. *Run the GUI:*

   ```bash
   python metadoon.py
   ```

---

## *📁 Project Structure*

```
Metadoon/
│
├── metadoon.py               # Main GUI script
├── Analise.R                  # R script for data analysis
├── generate_report.R          # Script to generate the final report
├── Metadoon_Report.Rmd        # RMarkdown template for the report
├── Metadoon-Beta.Rproj        # RStudio project file
├── metadoon_env.yaml          # Conda environment file
├── setup.sh                   # Environment setup script
├── LICENSE                    # License file
├── Readme.md                  # Project documentation
├── *.png, *.ico, *.icns       # Icons and GUI assets
│
├── Metadata/                  # Folder for metadata files
├── OTUs/                      # Folder for OTU tables
├── Taxonomy/                  # Folder for taxonomy files
├── Tree File/                 # Folder for phylogenetic tree files
├── Output/                    # Folder for generated results, plots, reports, tables, and the final report (HTML)
```

---

## *🗒️ How to Generate the Final Report*

- Inside the Metadoon interface, go to the **"Tools"** menu.
- Click **"Generate Final Report"**.
- This will run the R script that creates a complete report with all plots, alpha and beta diversity results, PERMANOVA, DESeq2 outputs, and summary.
- The report will be saved in the **Output/** folder as an HTML file.

## *💾 How to Save All Results to a Separate Folder*

- Inside the interface, go to **"Tools"**.
- Click **"Save and Clean Results"**.
- You will be prompted to select a folder where you want to save the results.
- Metadoon will move:
  - The **Output/** folder (containing all plots, tables, reports)
  - The parameter file **pipeline_params.json**
  - The **Rplots.pdf** file if generated
  - The **final report in HTML format**
- Once copied, the Output folder inside the project will be cleared.

---

## *📬 Contact*

*For issues or questions, please open an issue or contact the maintainer.*

[rdo.adan@gmail.com](mailto:rdo.adan@gmail.com)
