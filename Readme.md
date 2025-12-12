# 🧪 Metadoon

<div align="center">
  <img src="OP.png" alt="Metadoon Interface" width="600">
  
  <br><br>

  ![Docker](https://img.shields.io/badge/Docker-Supported-blue?style=flat&logo=docker)
  ![Conda](https://img.shields.io/badge/Conda-Supported-green?style=flat&logo=anaconda)
  ![Python](https://img.shields.io/badge/Python-3.10-yellow?style=flat&logo=python)
  ![R](https://img.shields.io/badge/R-Statistics-blue?style=flat&logo=r)
  ![License](https://img.shields.io/badge/license-MIT-lightgrey)

  <p>
    <b>User-friendly graphical interface and pipeline for amplicon-based metagenomic data analysis.</b>
  </p>
</div>

---

**Metadoon** automates the workflow from FASTQ preprocessing to robust statistical visualization in R, utilizing tools like **VSEARCH** and **Phyloseq**. It features a streamlined **5-step interface** and runs easily via **Docker** or **Natively** via Conda.

---

## 📦 What's Included

The environment includes:

| Component | Purpose |
| :--- | :--- |
| **Python 3.10** | GUI interface (Tkinter) and pipeline logic |
| **R (Latest)** | Statistical analysis and plotting |
| **VSEARCH** | FASTQ processing (merge, filter, cluster) |
| **Libraries** | `phyloseq`, `DESeq2`, `ggplot2`, `vegan`, etc. |

---

## 🚀 Option 1: One-Click Launchers
*Easy start scripts for all platforms.*

### ⚠️ First-Time Setup (Permissions)
**For macOS (`.command`) and Linux (`.sh`) users only:**
Before running the scripts for the first time, you must grant execution permissions via terminal.
1. Open a terminal inside the Metadoon folder.
2. Run the command:
   ```bash
   chmod +x *
   ```
*Note: **Windows users (`.bat`) DO NOT need this step.** You can run the file directly.*

### 1. Prerequisites by OS
* **Windows & Linux:** [Docker](https://www.docker.com/) installed (Enable WSL 2 for Windows).
* **macOS:** [Conda](https://www.anaconda.com/download) installed.
    * *The macOS `.command` launcher runs the **Native Conda** version, not Docker.*

### 2. How to Run
Just double-click the launcher for your OS:

* 🪟 **Windows:** Double-click `Windows_Run.bat` (Runs Docker).
* 🍎 **macOS:** Double-click `MacOS_Run.command` (Runs Conda/Native).
* 🐧 **Linux:** Run `./Linux_Run.sh` (Runs Docker).

---

## 🐍 Option 2: Manual Installation (Terminal)
*Recommended for Linux users or advanced users who prefer manual control.*

Follow these steps to run Metadoon directly on your system without the one-click scripts.

### 1. Prerequisites
* **Conda** (Anaconda or Miniconda) must be installed.

### 2. Installation & Execution
Open your terminal and run the following commands in order:

**Step 1: Clone the repository**
```bash
git clone https://github.com/rdo-adan/Metadoon.git
```

**Step 2: Enter the directory**
```bash
cd Metadoon/
```

**Step 3: Grant execution permissions**
Essential to ensure all scripts can run.
```bash
chmod +x *
```

**Step 4: Install dependencies**
This script creates the `metadoon` environment and installs R, Python, and VSEARCH.
```bash
bash setup.sh
```

**Step 5: Activate environment & Run**
```bash
conda activate metadoon
python metadoon.py
```

---

## 🖥️ Interface & Workflow

The new interface guides you through 5 simple steps:

1.  **Load FASTQ Files:** Select your raw data (must contain `_R1_` and `_R2_`).
2.  **Configure Parameters:** Adjust threads, max errors, and databases (optional).
3.  **RUN PIPELINE:** Starts the analysis (Merge -> Filter -> Cluster -> Taxonomy -> Stats).
4.  **Generate Report:** Creates the final HTML summary after the run finishes.
5.  **Save Results:** Exports all tables, plots, and reports to a clean folder.

---

## 📂 Handling Files (Docker Users)

If using Docker (Windows/Linux script), Metadoon maps your local folders:
* `/workspace` ⮕ **Metadoon folder** (Results saved here).
* `/app/YOUR_DATA` ⮕ **User Profile** (Documents, Downloads).
* `/app/C_Drive` ⮕ **C: Drive** (Windows only).

> **💡 Native/macOS Users:** You have direct access to your entire file system.

---

## ⚙️ Pipeline Details

1.  **Merge Pairs:** Merges R1 and R2 using VSEARCH.
2.  **Quality Filter:** Filters reads based on MaxEE.
3.  **Dereplication:** Identifies unique sequences.
4.  **Clustering:** OTU (97%) or ASV (Denoising).
5.  **Chimera Removal:** De novo + Reference-based.
6.  **Taxonomy:** SINTAX algorithm.
7.  **Statistics (R):** Alpha/Beta Diversity, Rarefaction, DESeq2, ANCOM-BC.

---

## 📁 Project Structure

Metadoon automatically manages file organization.

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
├── Windows_Run.bat          # Launcher scripts for Docker (All OS)
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

## ⚠️ Input Data Requirements

* **Format:** Illumina Paired-End `.fastq`.
* **Naming:** Must contain `_R1_` and `_R2_`.
* **No Special Characters:** Avoid spaces or extra hyphens in sample names.

---

## 📬 Contact

For issues or questions:
📧 [rdo.adan@gmail.com](mailto:rdo.adan@gmail.com)

