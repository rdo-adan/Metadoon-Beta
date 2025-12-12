# 🧪 Metadoon - Docker Edition

<div align="center">
  <img src="OP.png" alt="Metadoon Interface" width="600">
  
  <br><br>

  ![Docker](https://img.shields.io/badge/Docker-Containerized-blue?style=flat&logo=docker)
  ![Python](https://img.shields.io/badge/Python-3.10-yellow?style=flat&logo=python)
  ![R](https://img.shields.io/badge/R-Statistics-blue?style=flat&logo=r)
  ![VSEARCH](https://img.shields.io/badge/VSEARCH-%E2%89%A52.21.1-green)
  ![License](https://img.shields.io/badge/license-MIT-lightgrey)

  <p>
    <b>User-friendly graphical interface and pipeline for amplicon-based metagenomic data analysis.</b>
  </p>
</div>

---

**Metadoon** automates the workflow from FASTQ preprocessing to robust statistical visualization in R, utilizing tools like **VSEARCH** and **Phyloseq**.

### 🐳 **Why Docker Edition?**
This version runs in a fully containerized environment with **all dependencies pre-installed**. 
* ❌ No need to install Python, R, or VSEARCH manually.
* ✅ Compatible with **Windows** (via WSL2), **macOS**, and **Linux**.

---

## 📦 What's Included

All major and minor dependencies come pre-installed in the image:

| Component | Version | Purpose |
| :--- | :--- | :--- |
| **Python** | 3.10 | GUI interface (Tkinter) and pipeline logic |
| **R** | Latest | Statistical analysis and plotting |
| **VSEARCH** | ≥ 2.21.1 | FASTQ processing (merge, filter, cluster) |
| **Pillow** | Latest | Image handling in Python |

### 📊 R Packages (Pre-installed)
* **CRAN:** `tidyverse`, `ggplot2`, `ggpubr`, `pheatmap`, `viridis`, `ape`, `RColorBrewer`.
* **Bioconductor:** `phyloseq`, `DESeq2`, `scater`, `ANCOMBC`, `microbiome`.
* **GitHub:** `pairwiseAdonis`.

---

## 🚀 Quick Start (One-Click Launchers)

Metadoon provides shortcut scripts for a seamless experience. You don't need to type commands manually.

### 1. Prerequisites

* **Windows Users:** Install [Docker Desktop](https://www.docker.com/products/docker-desktop/).
    * *⚠️ Important:* Ensure **"Use WSL 2 based engine"** is checked in Docker settings.
* **macOS Users:**
    1.  Install [Docker Desktop](https://www.docker.com/products/docker-desktop/).
    2.  **Crucial:** Install [XQuartz](https://www.xquartz.org/) (required for the GUI).
    3.  *Setup:* Open XQuartz > Preferences > Security > Check **"Allow connections from network clients"**.
* **Linux Users:** Install [Docker Engine](https://docs.docker.com/engine/install/).

### 2. How to Run

First, ensure execute permissions (Unix systems only). In the root directory, run:
```bash
chmod +x *
```

#### 🪟 Windows
1.  Download and unzip this repository.
2.  Double-click **`Windows_Run.bat`**.
3.  A terminal will start the engine, and the GUI will appear shortly.

#### 🍎 macOS
1.  Double-click **`MacOS_Run.command`**.
    * *Note:* If macOS blocks execution, go to **System Settings > Privacy & Security** and allow the script.

#### 🐧 Linux
1.  Double-click or run via terminal:
    ```bash
    ./Linux_Run.sh
    ```

---

## 📂 Handling Files & Docker Mapping

When Metadoon launches via Docker, it maps your local folders to the container:

* `/workspace` ⮕ Maps to the **Metadoon folder** (project root). **All results are saved here.**
* `/app/YOUR_DATA` ⮕ Maps to your computer's **User Profile** (Documents, Downloads, Desktop).
* `/app/C_Drive` *(Windows Only)* ⮕ Maps to your C: drive.

> **💡 WINDOWS TIP: How to find your files**
> 1. Click **"Load FASTQ Files"** in the GUI.
> 2. Double-click the folder named **`YOUR_DATA`**.
> 3. This opens your user folder (`C:\Users\You`). Navigate to Desktop/Downloads from there.

---

## 🧪 Testing with Example Data

Validate your installation with our synthetic dataset.

1.  [**📥 Download Example Data (.zip)**](https://mega.nz/file/7ywinBYB#uaISNfE7d-9veK9earSEaI2vjR50CSByBKiHwgcToSU)
2.  Unzip the file.
3.  Open Metadoon.
4.  Click **"Load FASTQ Files"** and select the `.fastq` files.
5.  Select `metadata.tsv` when prompted.
6.  Click **"Run Pipeline"**.

---

## ⚙️ Pipeline Workflow

Metadoon executes a standard amplicon analysis workflow:

1.  **Merge Pairs:** Merges R1 and R2 FASTQ files using VSEARCH.
2.  **Quality Filter:** Filters reads based on Max Expected Error (`fastq_maxee`).
3.  **Dereplication:** Identifies unique sequences.
4.  **Clustering / Denoising:**
    * *OTU:* 97% clustering.
    * *ASV (ZOTU):* Denoising/Unoising (High resolution).
5.  **Chimera Removal:** Removes chimeric sequences (De novo + Reference-based).
6.  **Taxonomy Assignment:** Uses SINTAX algorithm.
7.  **Statistical Analysis (R):**
    * Rarefaction Curves.
    * **Alpha Diversity:** Richness/Evenness (Kruskal-Wallis/ANOVA).
    * **Beta Diversity:** NMDS/PCoA (PERMANOVA).
    * **Differential Abundance:** DESeq2 and ANCOM-BC.

---

## 📁 Project Structure

Metadoon organizes your workspace automatically.

```text
Metadoon/
├── DB/                  # Reference databases (RDP, Silva)
├── Metadata File/       # Your metadata
├── Tree File/           # Phylogenetic tree (optional)
│
├── Merged/              # Merged reads
├── Filtered/            # Quality filtered sequences
├── Dereplicated/        # Unique sequences
│
├── OTUs/                # Clustering results
│   ├── otus.fasta       # Final Sequences (OTUs/ASVs)
│   └── otutab.txt       # Abundance table
│
├── Taxonomy/            # Classification results
│   └── taxonomy.txt     # Cleaned taxonomy for R
│
├── Output/              # 🏆 FINAL RESULTS
│   ├── Plots/           # Heatmaps, PCoA, Alpha Div, Rarefaction
│   └── Stats/           # DESeq2, ANCOM-BC tables
│
└── Metadoon_Report.html # Complete HTML Summary
```

---

## ⚠️ Input Data Requirements

To avoid errors, ensure your data meets these criteria:

* **Platform:** Illumina Paired-End sequencing.
* **Format:** `.fastq` files.
* **Naming Convention:**
    * Must contain `_R1_` and `_R2_`.
    * *Recommended:* `Sample-Name_R1.fastq`.
    * ❌ **Avoid:** Spaces, special characters, or extra hyphens/underscores in sample names.

---

## 🛠️ Advanced: Native Installation (Linux Only)

If you prefer running without Docker using Conda, click below.

<details>
<summary><b>🔻 Click to expand Native Installation Instructions</b></summary>

This method is recommended for Linux users who want to run tools natively.

**Prerequisites:** Conda (Anaconda or Miniconda).

**1. Clone the repository:**
```bash
git clone https://github.com/rdo-adan/Metadoon.git
cd Metadoon/
```

**2. Install Dependencies:**
This script creates the `metadoon` environment and installs Python, R, and VSEARCH.
```bash
bash setup.sh
```

**3. Run:**
```bash
conda activate metadoon
python metadoon.py
```
</details>

---

## 🗒️ Tools & Reports

### Generating the Report
1.  Run the pipeline until finished.
2.  Go to **Tools > Generate Final Report**.
3.  The `Metadoon_Report.html` will open automatically.

### Saving Results
1.  Go to **Tools > Save and Clean Results**.
2.  Select a destination folder.
3.  Metadoon creates a timestamped archive (e.g., `Metadoon_Results_2025-10-20`) and can optionally delete temporary files (`Merged`, `Filtered`) to save space.

---

## 📬 Contact

For issues, bugs, or questions, please open an Issue on GitHub or contact the maintainer.

📧 **Maintainer:** [rdo.adan@gmail.com](mailto:rdo.adan@gmail.com)
