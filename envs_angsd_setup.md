# ANGSD v0.935 — Environment Setup

This document describes how to install ANGSD v0.935 and supporting tools used in Script 01 (SNP calling) of the *Rumex hastatus* phylogenomics pipeline.

---

## Option A — Conda installation (recommended)

ANGSD, BWA, and SAMtools are all available through the [bioconda](https://bioconda.github.io/) channel.

### 1. Set up bioconda channels (one-time setup)

If you haven't configured bioconda before, run:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

### 2. Create a dedicated environment

```bash
conda create -n rumex-snp -c bioconda \
    angsd=0.935 \
    bwa=0.7.17 \
    samtools=1.17 \
    -y
```

### 3. Activate the environment

```bash
conda activate rumex-snp
```

### 4. Verify installation

```bash
angsd --version       # should report 0.935
bwa 2>&1 | head -3    # BWA usage info
samtools --version | head -1
```

---

## Option B — Manual installation from source

If conda is not available, compile ANGSD v0.935 from source:

```bash
# Install dependencies (Ubuntu/Debian)
sudo apt install build-essential zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev

# Download and compile ANGSD v0.935
wget https://github.com/ANGSD/angsd/archive/refs/tags/0.935.tar.gz
tar -xzf 0.935.tar.gz
cd angsd-0.935
make

# Add to PATH (replace /path/to with your install location)
export PATH=/path/to/angsd-0.935:$PATH
```

---

## Option C — HPC module system

If you are on an HPC cluster, ANGSD may already be installed as a module:

```bash
module avail angsd
module load angsd/0.935
```

---

## Notes

- All other pipeline tools (SMC++, TreeMix, ADMIXTURE, Circuitscape) have their own conda environments — see `envs/conda_env.yml` for the full setup.
- The `rumex-snp` environment is isolated so you can manage SNP-calling tools without affecting other analyses.
