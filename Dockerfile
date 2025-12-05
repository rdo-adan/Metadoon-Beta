# Use base Miniconda image
FROM continuumio/miniconda3

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libgmp3-dev \
    libmpfr-dev \
    pandoc \
    libx11-6 \
    libxext6 \
    libxrender1 \
    libxtst6 \
    tk \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy environment file
COPY metadoon_env.yaml .

# Create environment
RUN conda env create -f metadoon_env.yaml

# Add environment to PATH
ENV PATH /opt/conda/envs/metadoon/bin:$PATH

# Install extra R packages
RUN Rscript -e "if (!requireNamespace('devtools', quietly = TRUE)) install.packages('devtools', repos='https://cloud.r-project.org')" && \
    Rscript -e "devtools::install_github('pmartinezarbizu/pairwiseAdonis/pairwiseAdonis')"

# Copy source code
COPY . .

# Default command
CMD ["python", "metadoon.py"]