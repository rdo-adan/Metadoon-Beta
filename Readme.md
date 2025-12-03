# 🧪 Metadoon
# 🧪 Metadoon

<div align="center">
  <img src="OP.png" alt="Metadoon Interface" width="50%">
</div>


***Metadoon*** is a user-friendly graphical interface and pipeline designed for processing and analyzing amplicon-based metagenomic data using tools like VSEARCH and R (with Phyloseq). It automates the workflow from FASTQ preprocessing to statistical visualization in R.

---

## *📦 Major Dependencies*

| *Dependency* | *Version (Suggested)* | *Description* |
| --------------------------------------------------------- | --------------------- | ---------------------------------------------- |
| *[Python](https://www.python.org/downloads/)* | *3.12+* | *Main interface (Tkinter GUI, logic control)* |
| *[R](https://cran.r-project.org/)* | *4.4.3* | *Statistical analysis and plotting* |
| *[VSEARCH](https://github.com/torognes/vsearch/releases)* | *≥ 2.21.1* | *FASTQ processing (dereplication, clustering)* |

---

## *🐍 Python Dependencies (minor)*

*These packages are included in the Conda environment:*

| *Package* | *Purpose* |
| ---------- | ------------------------------------------------------ |
| *`tkinter`*| *GUI interface (Standard Python Library)* |
| *`Pillow`* | *Icon/image handling in Python* |

> \*Standard libraries used: **`os`**, **`sys`**, **`json`**, **`threading`**, **`subprocess`**, **`shutil`**, **`glob`**, **`datetime`**

---

## *📊 R Dependencies (minor)*

All packages listed below are automatically installed by the Conda environment.

### *📦 CRAN Packages*

* tidyverse, reshape2, igraph, foreach, lme4
* ggplot2, ggpubr, cowplot, dplyr, pheatmap, viridis
* ape, rprojroot, wesanderson, RColorBrewer

### *🧪 Bioconductor Packages*

* phyloseq, DESeq2, scater

### *🔧 GitHub Packages (installed post-env)*
# R
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("vlubitch/pairwiseAdonis")
devtools::install_github("microbiome/microbiome")
## *🛠️ System-Level Dependencies*

*These are installed via Conda or available on Unix-based systems:*

| *Tool/Library*                                           | *Description*                       |
| -------------------------------------------------------- | ----------------------------------- |
| `bash`, `wget`                                           | *Script automation and downloading* |
| `conda`, `mamba`                                         | *Environment management*            |
| `Rscript`                                                | *Execute R scripts via CLI*         |
| `libcurl`, `libxml2`                                     | *R package compilation*             |
| pandoc                                                   | Report generation (HTML)            |
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

# Option 1: Using Easy Launchers (Recommended for Non-Technical Users)
We provide launcher scripts to simplify installation and execution without typing commands manually.

> Download and Extract:

- Download the repository as a ZIP file or clone it.

- Extract the folder to a location of your choice.

> Install Dependencies:

- Open the Installers/ folder.

- Windows: Double-click Windows_Install.bat.

- macOS: Double-click MacOS_Install.command.

- Linux: Open a terminal and run bash Linux_Install.sh.

- Wait for the installation to complete (this may take a few minutes).

> Run Metadoon:

- Open the Run/ folder.

- Windows: Double-click Windows_Run.bat.

- macOS: Double-click MacOS_Run.command.

- Linux: Run bash Linux_Run.sh.

-- The graphical interface will open automatically.

# Option 2: Using Command Line (Advanced)
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
## ⚙️ Pipeline Workflow
Metadoon executes a standard amplicon analysis workflow:

- Merge Pairs: Merges R1 and R2 FASTQ files using VSEARCH.

- Quality Filter: Filters reads based on maximum expected error (fastq_maxee).

- Dereplication: Identifies unique sequences to reduce computational load.

-- Clustering: Clusters sequences into OTUs (default 97% identity)/ ASV (ZOTU): Performs denoising (unoising) to resolve exact biological sequences.

- Chimera Removal: Removes chimeric sequences using both de novo and Reference-based detection.

- Taxonomy Assignment: Assigns taxonomy using the SINTAX algorithm.

- Statistical Analysis (R):

- Rarefaction curves.

- Alpha & Beta Diversity metrics.

- Core Microbiome analysis.

- Differential Abundance (DESeq2).

- Analysis of Compositions of Microbiomes with Bias Correction (ANCOM-BC).

## *📁 Project Structure*
Before run:
```
Metadoon/
│
├── metadoon.py              # Main GUI script (Python)
├── Analise.R                # Statistical analysis script (R)
├── generate_report.R        # Report generation script
├── Metadoon_Report.Rmd      # RMarkdown template
├── pipeline_params.json     # Configuration file
├── metadoon_env.yaml        # Conda environment definition
├── setup.sh                 # Installation script
├── LICENSE                  # License file
├── Readme.md                # Project documentation
├── *.png, *.ico, *.icns     # Icons and GUI assets
│
├── Installers/              # Scripts to install dependencies easily
│   ├── Windows_Install.bat
│   ├── MacOS_Install.command
│   └── Linux_Install.sh
│
├── Run/                     # Scripts to launch the tool easily
│   ├── Windows_Run.bat
│   ├── MacOS_Run.command
│   └── Linux_Run.sh
│
└── Example_Data.txt           # Links to Download Small dataset for testing the tool
```
After Run
```
Metadoon/
│
├── DB/                      # Downloaded reference databases (RDP, Silva, etc.)
├── Metadata File/           # Stores the uploaded metadata file
├── Tree File/               # Stores the phylogenetic tree (if provided)
│
├── Merged/                  # Paired-end reads merged by VSEARCH
├── FullFiles/               # Concatenated merged reads
├── Filtered/                # Quality filtered sequences
├── Dereplicated/            # Unique sequences (dereplication)
│
├── OTUs/                    # Clustering results
│   ├── centroids.fasta      # Representative sequences
│   ├── otus.fasta           # Final OTUs/ASVs (non-chimeric)
│   └── otutab.txt           # Abundance table
│
├── Taxonomy/                # Taxonomic classification results
│   ├── taxonomy_raw.txt     # Raw output from SINTAX
│   └── taxonomy.txt         # Cleaned taxonomy table for R
│
└── Output/                  # FINAL RESULTS
    ├── Plots (Alpha/Beta diversity, Heatmaps, Rarefaction)
    ├── Statistical Tables (DESeq2, ANCOM-BC, PERMANOVA)
    └── Metadoon_Report.html # Complete HTML Summary
│
├── metadoon.py              # Main GUI script (Python)
├── Analise.R                # Statistical analysis script (R)
├── generate_report.R        # Report generation script
├── Metadoon_Report.Rmd      # RMarkdown template
├── pipeline_params.json     # Configuration file
├── metadoon_env.yaml        # Conda environment definition
├── setup.sh                 # Installation script
├── LICENSE                  # License file
├── Readme.md                # Project documentation
├── *.png, *.ico, *.icns     # Icons and GUI assets
│
├── Installers/              # Scripts to install dependencies easily
│   ├── Windows_Install.bat
│   ├── MacOS_Install.command
│   └── Linux_Install.sh
│
├── Run/                     # Scripts to launch the tool easily
│   ├── Windows_Run.bat
│   ├── MacOS_Run.command
│   └── Linux_Run.sh
│
└── Example_Data.txt           # Links to Download Small dataset for testing the tool
```
---
## *🗒️ How to Run the Pipeline
## *🏃 How to Run the Pipeline*

### *⚠️ Input Data Requirements (Read Carefully)*

Before running Metadoon, ensure your files meet the following criteria to avoid errors:

1.  **Sequencing Platform:** Currently, Metadoon only supports **Illumina Paired-End** sequencing data (Forward and Reverse reads).
2.  **File Format:** Files must be in **`.fastq`** format.
3.  **Naming Convention:**
    * Files **must** contain `_R1_` (for forward) and `_R2_` (for reverse) identifiers to be recognized by the merger script.
    * **Sample Naming:** It is highly recommended to rename your files using a `Sample-Replica` structure before the R1/R2 tag.
    * *Example:* `M1-S1_R1.fastq` and `M1-S1_R2.fastq`.
4.  **⛔ File Naming Caution:**
    * **Avoid using extra hyphens (`-`)** or special characters in your sample names (e.g., avoid `Sample-Group-A-1_R1.fastq`).
    * Excessive hyphens can interfere with how R and some plotting libraries interpret group names.
    * **Recommended:** Use underscores (`_`) or alphanumeric characters for complex names (e.g., `SampleGroupA_1_R1.fastq`).

---

### *🖥️ Step-by-Step Execution*

1.  **Load FASTQ Files:**
    * Click on **"1. Load FASTQ Files"**.
    * Navigate to your data folder.
    * Select **all** your `.fastq` files (both R1 and R2) at once and click "Open".
    * The files will appear in the "Loaded Files" list on the right.

2.  **Configure Parameters:**
    * Click on **"2. Configure Parameters"**.
    * Adjust the processing settings according to your data needs (see the *Configuration Parameters Explained* section below for details).
    * **Clustering Method:** Choose between **OTU** (classic) or **ASV** (modern/denoising).
    * Click **"Save Parameters"** and close the window.

3.  **Run the Pipeline:**
    * Click the blue button **"3. RUN PIPELINE"**.
    * A processing window will appear. **Do not close the application.**
    * You can monitor the progress in the "Pipeline Log" terminal at the bottom right.
    * Wait for the message: *"Pipeline Finished"*.

4.  **Analyze & Report:**
    * Once the pipeline finishes, click **"4. Generate Report"** to view the HTML summary.

---

### *⚙️ Configuration Parameters Explained*

In the **"Configure Parameters"** window, you can fine-tune how Metadoon processes your data.

#### **1. General VSEARCH Parameters**
* **Threads:** The number of CPU cores to use. Higher values speed up processing.
* **Max Diffs (Merge):** Maximum number of mismatched bases allowed in the overlap region when merging R1 and R2 reads. Default is `30`.
* **Max EE (Filter):** "Maximum Expected Error". A quality filtering threshold. Reads with a cumulative error probability higher than this value are discarded. Lower values (e.g., `0.5` or `1.0`) are stricter and produce higher quality data.
* **Min Unique Size:** Minimum abundance for a sequence to be kept during dereplication. Default is `2` (removes singletons, which are often sequencing errors).
* **Analysis Type (Clustering):**
    * **OTU (Cluster 97%):** Traditional method grouping sequences with 97% similarity. Good for general diversity trends.
    * **ASV (Denoising/Unoising):** Uses a denoising algorithm to resolve single-nucleotide differences (ZOTUs). Provides higher resolution and accuracy.
* **Identity %:** The similarity threshold for OTU clustering (ignored if ASV is selected). Default is `0.97`.
* **SINTAX Cutoff:** Confidence threshold for taxonomic assignment (0.0 to 1.0). Default `0.8` means only classifications with ≥80% confidence are kept.
* **Strand:** Direction of reads to consider for taxonomy (`plus` or `both`). `both` is safer if orientation is unknown.

#### **2. Reference Databases**
* **Chimera DB:** Database used to detect and remove chimeric sequences (PCR artifacts). Default is `RDP Gold`.
* **16S Database:** Reference database for taxonomy assignment.
    * **RDP:** General purpose 16S database.
    * **Greengenes2:** Updated database aligned with genomic trees.
    * **Custom:** Allows you to upload your own `.fasta` database (e.g., Silva).

#### **3. R Analysis Parameters**
* **Color Palette:** Selects the color scheme for all plots (e.g., `viridis`, `plasma`, `magma`). These are color-blind friendly.
* **Rarefaction Settings:**
    * **Enable Rarefaction:** If checked, subsamples all libraries to the same depth (normalizes uneven sequencing effort).
    * **Depth:** The number of reads to subsample per sample (e.g., `1000`). Samples below this count are discarded.
* **Plot Limits (Top N):**
    * **Abundance Top N:** Number of most abundant taxa to show in bar plots (e.g., `15`).
    * **Core Top N:** Number of taxa to display in Core Microbiome heatmaps.

## *🗒️ How to Generate the Final Report
- Inside the Metadoon interface; Click "Generate Report".

- This will run the R script that creates a complete report with all plots, alpha and beta diversity results, PERMANOVA, DESeq2 outputs, and summary.

- The report will be saved in the root folder as Metadoon_Report.html.

## *💾 How to Save All Results
- Inside the interface, Click "Save Results".

- You will be prompted to select a destination folder.

- Metadoon will create a new folder named Metadoon_Results_YYYY-MM-DD_HH-MM-SS inside your selected directory.

- It copies all critical outputs (Output/, OTUs/, Taxonomy/, Reports) to this safe location.

- Optional Cleanup: After saving, the tool will ask if you want to delete the generated workspace folders (Merged, Filtered, DB, etc.) to free up disk space.

---

## *📬 Contact*

*For issues or questions, please open an issue or contact the maintainer.*

[rdo.adan@gmail.com](mailto:rdo.adan@gmail.com)
