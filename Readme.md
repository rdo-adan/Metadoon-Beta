# 🧪 Metadoon

***Metadoon***\* is a user-friendly graphical interface and pipeline designed for processing and analyzing amplicon-based metagenomic data using tools like VSEARCH and R (with Phyloseq). It provides automation from FASTQ preprocessing to statistical visualization using R.\*

---

## *📦 Major Dependencies*

| *Dependency*                                              | *Version (Suggested)* | *Description*                                  |
| --------------------------------------------------------- | --------------------- | ---------------------------------------------- |
| [*Python*](https://www.python.org/downloads/)             | *3.12+*               | *Main interface (Tkinter GUI, logic control)*  |
| [*R*](https://cran.r-project.org/)                        | *≥ 4.2*               | *Statistical analysis and plotting*            |
| [*VSEARCH*](https://github.com/torognes/vsearch/releases) | *≥ 2.21.1*            | *FASTQ processing (dereplication, clustering)* |

---

## *🐍 Python Dependencies (minor)*

*These packages can be installed via **`pip`** or managed via Conda:*

```bash
pip install pillow tk
```

| *Package*                                                            | *Purpose*                       |
| -------------------------------------------------------------------- | ------------------------------- |
| *`tk`*                                                               | *GUI interface with Tkinter*    |
| *`pillow`*                                                           | *Icon/image handling in Python* |
| \*`os`\*\*, **`sys`**, **`json`**, **`threading`**, \**`subprocess`* | *Python standard libraries*     |

---

## *📊 R Dependencies (minor)*

### *📦 CRAN Packages*

```r
install.packages(c("tidyverse", "reshape2", "igraph", "foreach", "lme4",
                   "ggplot2", "ggpubr", "cowplot", "dplyr", "pheatmap", "viridis", "ape", "rprojroot"))
```

### *🧪 Bioconductor Packages*

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("phyloseq", "DESeq2", "scater"))
```

### *🔧 GitHub Packages*

```r
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("vlubitch/pairwiseAdonis")
devtools::install_github("microbiome/microbiome")
```

---

## *🛠️ System-Level Dependencies*

*Ensure these tools and libraries are available in your OS:*

| *Tool/Library*                                               | *Description*                            |
| ------------------------------------------------------------ | ---------------------------------------- |
| \*`bash`\*\*, \**`wget`*                                     | *Required for scripting and downloading* |
| *`conda`*                                                    | *Recommended for environment management* |
| *`Rscript`*                                                  | *Required to run R analysis from CLI*    |
| \*`libcurl`\*\*, **`zlib`**, **`libxml2`**, \**`libssl-dev`* | *Required for compiling R packages*      |
| \*`gcc`\*\*, \**`make`*                                      | *Compilation tools*                      |

> *On macOS, **`XQuartz`** may be required for full graphical R support.*

---

## *🔗 Download Links*

* ***Python***\*: \*[*https://www.python.org/downloads/*](https://www.python.org/downloads/)
* ***R***\*: \*[*https://cran.r-project.org/*](https://cran.r-project.org/)
* ***VSEARCH***\*: \*[*https://github.com/torognes/vsearch/releases*](https://github.com/torognes/vsearch/releases)
* ***Miniconda***\* (recommended): \*[*https://docs.conda.io/en/latest/miniconda.html*](https://docs.conda.io/en/latest/miniconda.html)

---

## *🚀 Usage*

1. *Clone this repository:*

   ```bash
   git clone https://github.com/seuusuario/metadoon.git
   cd metadoon
   ```

2. *Create environment:*

   ```bash
   conda env create -f environment.yaml
   conda activate metadoon
   ```

3. *Run the GUI:*

   ```bash
   python metadoon.py
   ```

---

## *📁 Project Structure*

```
Metadoon/
│
├── metadoon.py             # Main GUI
├── Analise.R               # R script for data analysis
├── install_r_packages.R    # (optional) R package installer
├── environment.yaml        # Conda environment file
├── metadoon.yaml           # Conda base config
├── *.png, *.ico, *.icns    # Icons and GUI assets
└── Output/, Metadata/, OTUs/, Taxonomy/, Tree File/
```

---

## *📬 Contact*

*For issues or questions, please open an issue or contact the maintainer.*
