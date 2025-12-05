# 🧪 Metadoon - Docker Edition

<div align="center">
  <img src="OP.png" alt="Metadoon Interface" width="600">
</div>

***Metadoon*** is a user-friendly graphical interface and pipeline designed for processing and analyzing amplicon-based metagenomic data using tools like **VSEARCH** and **R** (with Phyloseq). It automates the workflow from FASTQ preprocessing to robust statistical visualization in R.

### 🐳 **DOCKER VERSION**
This version runs in a fully containerized environment with all dependencies pre-installed. No need to install Python, R, VSEARCH, or any libraries manually. Compatible with **Linux** and **Windows** (via WSL2).

---

## *📦 What's Included in the Docker Image*

All major and minor dependencies are pre-installed:

| Component | Version | Purpose |
| :--- | :--- | :--- |
| **Python** | 3.10 | GUI interface (Tkinter) and pipeline logic |
| **R** | Latest | Statistical analysis and plotting |
| **VSEARCH** | ≥ 2.21.1 | FASTQ processing (merge, filter, cluster) |
| **Pillow** | Latest | Image handling in Python |

### *📊 R Packages (Pre-installed)*
* **CRAN:** `tidyverse`, `ggplot2`, `ggpubr`, `pheatmap`, `viridis`, `ape`, `RColorBrewer`, and more.
* **Bioconductor:** `phyloseq`, `DESeq2`, `scater`, `ANCOMBC`, `microbiome`.
* **GitHub:** `pairwiseAdonis`.

---

## *🚀 Installation & Usage*

You can run Metadoon using **Docker (Recommended)** or via **Native Installation (Conda)**.

### *Option 1: Docker Container (Recommended & Cross-Platform)*
This method runs Metadoon in an isolated container. It works on **Linux**, **Windows (via WSL2)**, and **macOS**. No need to install R, Python, or VSEARCH manually.

#### **Prerequisites**
1.  **[Docker Desktop](https://www.docker.com/products/docker-desktop/)** installed and running.
2.  **Windows Users:** Ensure **[WSL 2](https://learn.microsoft.com/en-us/windows/wsl/install)** is installed and Docker is configured to use it.
3.  **Graphics Support:**
    * *Linux:* Usually works out of the box.
    * *Windows:* WSL 2 usually handles graphics (WSLg). If not, use an XServer like **VcXsrv**.
    * *macOS:* Install **[XQuartz](https://www.xquartz.org/)** and allow network connections.

#### **How to Run**
1.  Download this repository.
2.  Open the `Run/` folder.
3.  Double-click (or run in terminal) the script for your OS:
    * 🪟 **Windows:** `Windows_Run.bat`
    * 🐧 **Linux:** `Linux_Run.sh`
    * 🍎 **macOS:** `MacOS_Run.command`

> **ℹ️ Under the Hood:** These scripts automatically build the image (if missing) and execute the following command to map your files and display the GUI:
> ```bash
> xhost +local:docker
>
> docker run --rm -it \
>   --user $(id -u):$(id -g) \
>   -e DISPLAY=$DISPLAY \
>   -v /tmp/.X11-unix:/tmp/.X11-unix \
>   -v "$(pwd)":/workspace:rw \
>   -v "$HOME":/host_home:ro \
>   metadoon
> ```

---

### *Option 2: Native Installation (Linux Only)*
This is the "classic" method using Conda. It is recommended primarily for **Linux** users who prefer running tools natively without Docker.

#### **Prerequisites**
* **Conda** (Anaconda or Miniconda) installed.

#### **Step-by-Step**
1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/rdo-adan/metadoon.git](https://github.com/rdo-adan/metadoon.git)
    cd metadoon
    ```

2.  **Run the Setup Script:**
    This script creates the `metadoon` environment and installs all dependencies (Python, R, VSEARCH, Libraries).
    ```bash
    bash setup.sh
    ```

3.  **Activate and Run:**
    ```bash
    conda activate metadoon
    python metadoon.py
    ```

---

## *⚙️ Pipeline Workflow*

Metadoon executes a standard amplicon analysis workflow:

1.  **Merge Pairs:** Merges R1 and R2 FASTQ files using VSEARCH.
2.  **Quality Filter:** Filters reads based on maximum expected error (`fastq_maxee`).
3.  **Dereplication:** Identifies unique sequences to reduce computational load.
4.  **Clustering / Denoising:**
    * **OTU:** Clusters sequences at 97% identity.
    * **ASV (ZOTU):** Performs denoising (unoising) to resolve exact biological sequences.
5.  **Chimera Removal:** Removes chimeric sequences using both *de novo* and Reference-based detection.
6.  **Taxonomy Assignment:** Assigns taxonomy using the SINTAX algorithm.
7.  **Statistical Analysis (R):**
    * **Rarefaction:** Curves to assess sequencing depth.
    * **Alpha Diversity:** Richness/Evenness metrics with Kruskal-Wallis/ANOVA tests.
    * **Beta Diversity:** NMDS/PCoA plots with PERMANOVA statistics.
    * **Differential Abundance:** Significant taxa identified by **DESeq2** and **ANCOM-BC** (robust to compositional bias).

---

## *📁 Project Structure*

Metadoon automatically manages file organization. Below is the structure before and after execution:

### *Core Files (Before Run)*
```text
Metadoon/
│
├── metadoon.py              # Main GUI script
├── Analise.R                # Statistical analysis script (R)
├── generate_report.R        # Report generation script
├── Metadoon_Report.Rmd      # RMarkdown template
├── Dockerfile               # Docker configuration
├── pipeline_params.json     # Configuration file
├── metadoon_env.yaml        # Conda environment definition
├── setup.sh                 # Native installation script (Linux)
├── LICENSE                  # License file
├── Readme.md                # Project documentation
│
├── Run/                     # Launcher scripts for Docker (All OS)
│   ├── Windows_Run.bat
│   ├── MacOS_Run.command
│   └── Linux_Run.sh
│
└── Example_Data.txt            # Small dataset for testing

```

### *Generated Directories (After Run)*
Once the pipeline runs, Metadoon creates specific folders to organize the workflow:

```text
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
```

---

## *📁 Input Data Requirements*

### *⚠️ Critical Requirements*
* **Platform:** Illumina Paired-End sequencing only.
* **Format:** `.fastq` files.
* **Naming Convention:**
    * Must contain `_R1_` (forward) and `_R2_` (reverse).
    * **Recommended:** `Sample-Replica_R1.fastq` format.
    * *Example:* `M1-S1_R1.fastq` and `M1-S1_R2.fastq`.
* **⛔ Avoid:**
    * Extra hyphens in sample names.
    * Special characters.
    * Use underscores (`_`) for complex names.

---

## *🗒️ How to Generate the Final Report*

1.  **Run Pipeline:** Ensure you have loaded FASTQ files and clicked "RUN PIPELINE". Wait for the "Pipeline Finished" message.
2.  **Generate:** Go to the **"Tools"** menu > Click **"Generate Final Report"**.
3.  **View:** The `Metadoon_Report.html` will be created in the root folder and opened automatically.

## *💾 How to Save All Results*

1.  **Save:** Go to **"Tools"** > Click **"Save and Clean Results"**.
2.  **Select Folder:** Choose a destination directory.
3.  **Automatic Organization:** Metadoon creates a timestamped folder (e.g., `Metadoon_Results_2025-10-20_14-30`) containing all critical outputs.
4.  **Cleanup:** You can choose to automatically delete the working directories (`Merged`, `Filtered`, etc.) to free up space.

---

## *⚙️ Configuration Parameters Explained*

In the **"Configure Parameters"** window, you can fine-tune how Metadoon processes your data.

* **Threads:** The number of CPU cores to use.
* **Max Diffs (Merge):** Maximum mismatches allowed in overlap region (Default: `30`).
* **Max EE (Filter):** "Maximum Expected Error". Lower values (e.g., `1.0`) are stricter.
* **Min Unique Size:** Minimum abundance to keep a sequence (Default: `2`).
* **Analysis Type:**
    * **OTU:** 97% clustering.
    * **ASV:** Denoising/Unoising (High resolution).
* **Identity %:** Similarity threshold for OTU (ignored for ASV).
* **SINTAX Cutoff:** Taxonomy confidence (Default `0.8`).
* **Databases:** Supports RDP, Greengenes2, and Custom databases.

---

## *📬 Contact*

*For issues or questions, please open an issue or contact the maintainer.*

[rdo.adan@gmail.com](mailto:rdo.adan@gmail.com)
