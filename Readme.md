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

Metadoon operates entirely locally using **Docker**. This ensures that all dependencies (R, Python, VSEARCH) work perfectly on your machine without complex manual installation.

### 📋 1. Prerequisites

Before running Metadoon, you must have **Docker** installed and running:

* **Windows Users:** Install [Docker Desktop for Windows](https://www.docker.com/products/docker-desktop/).
    * *Note:* Ensure "Use WSL 2 based engine" is checked in Docker settings.
* **macOS Users:**
    1.  Install [Docker Desktop for Mac](https://www.docker.com/products/docker-desktop/).
    2.  **Crucial:** Install [XQuartz](https://www.xquartz.org/) to allow the graphical interface to appear.
        * *After installing XQuartz, go to Preferences > Security and check "Allow connections from network clients".*
* **Linux Users:** Install [Docker Engine](https://docs.docker.com/engine/install/).

---

### 🏃 2. How to Run (One-Click Launchers)

You do **not** need to open a terminal or type complex Docker commands manually. We provide automatic launchers for every system.

#### 🪟 Windows
1.  Download and unzip this repository.
2.  Navigate to the `Run/` folder.
3.  Double-click **`Windows_Run.bat`**.
4.  A terminal window will open to start the engine, and the Metadoon GUI will appear shortly.

#### 🍎 macOS
1.  Open **XQuartz** first.
2.  Navigate to the `Run/` folder.
3.  Right-click **`MacOS_Run.command`** and select *Open*.
    * *Note: If macOS prevents execution due to security, go to System Settings > Privacy & Security and allow the script.*

#### 🐧 Linux
1.  Open a terminal in the `Run/` folder.
2.  Make the script executable (only needed once):
    ```bash
    chmod +x Linux_Run.sh
    ```
3.  Run the script:
    ```bash
    ./Linux_Run.sh
    ```

---

### 🧪 3. Testing with Example Data

To ensure everything is working correctly, you can use our lightweight synthetic dataset.

1.  [**📥 Download Example Data (.zip)**](https://mega.nz/file/7ywinBYB#uaISNfE7d-9veK9earSEaI2vjR50CSByBKiHwgcToSU)
2.  Unzip the file.
3.  Open Metadoon.
4.  Click **"Load FASTQ Files"** and select the `.fastq` files from the extracted folder.
5.  When prompted for Metadata, select `metadata.tsv`.
6.  Click **"Run Pipeline"**.

---

### 📂 4. Where are my files? (Docker Mapping)

When Metadoon launches via these scripts, it automatically maps your local folders into the container for easy access:

* **/workspace**: This is the folder where you ran the script (the project root). **All results (Output folder) are saved here.**
* **/app/YOUR_DATA**: This maps to your computer's **User Profile** (Documents, Downloads, Desktop) for easy file finding.
* **/app/C_Drive** *(Windows Only)*: Direct access to your entire C: drive.

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
├── Readme.md                # Project documentation                    # Launcher scripts for Docker (All OS)
├── Windows_Run.bat
├── MacOS_Run.command
├── Linux_Run.sh
└── Example_Data.txt         # Links to Download a Small dataset for testing
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
    * Try aways to Use your Metadata file un tsv format.
* **⛔ Avoid:**
    * Extra hyphens in sample names.
    * Special characters.
    * Use underscores (`_`) for complex names.

---

---
### *📂 Windows Users: How to Find Your Files (Docker)*

When running via Docker on Windows, the file selection window works slightly differently. We have created a shortcut to make it easy to find your Windows files.

**Steps to Load Files:**

1.  Click **"1. Load FASTQ Files"**.
2.  In the window that opens, you will see a folder named **`YOUR_DATA`**.
3.  **Double-click `YOUR_DATA`**.
    * This automatically opens your Windows User folder (`C:\Users\YourName\`).
    * From there, you can navigate to your **Desktop**, **Documents**, or **Downloads** to select your FASTQ files.

> 💡 **Tip:** If you placed your FASTQ files inside the `Metadoon` folder itself, they will appear right away in the initial list, next to the `Run` folder.
---


## *🗒️ How to Generate the Final Report*

1.  **Run Pipeline:** Ensure you have loaded FASTQ files, Configured your parameters and clicked "RUN PIPELINE". Wait for the "Pipeline Finished" message.
2.  **Generate:** Go to the **"Tools"** menu > Click **"Generate Final Report"**.
3.  **View:** The `Metadoon_Report.html` will be created in the root folder and opened automatically.

## *💾 How to Save All Results*

1.  **Save:** Go to **"Tools"** > Click **"Save and Clean Results"**.
2.  **Select Folder:** Choose a destination directory.
3.  **Automatic Organization:** Metadoon creates a timestamped folder (e.g., `Metadoon_Results_2025-10-20_14-30`) containing all critical outputs.
4.  **Cleanup:** You can choose to automatically delete the working directories (`Merged`, `Filtered`, etc.) to free up space after saving.

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
