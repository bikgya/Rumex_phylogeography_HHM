# Script 01 — SNP Calling

This folder contains the scripts for nuclear SNP discovery in *Rumex hastatus*.

## Pipeline overview

```
Cleaned reads  ──► BWA-MEM mapping  ──► Sorted BAM files  ──► ANGSD SNP calling
                  (01a_bam_generation)   (per sample)         (01b_angsd_snp_calling)
```

## Scripts

| Script | Description | Input | Output |
|--------|-------------|-------|--------|
| `01a_bam_generation.sh` | Map paired-end cleaned reads to the reference genome with BWA-MEM, sort and index with SAMtools | Cleaned FASTQ files + reference | Sorted, indexed BAM files |
| `01b_angsd_snp_calling.sh` | Identify nuclear SNPs across all samples using ANGSD (GATK genotype likelihood model) | BAM list + reference | MAF, posterior genotypes, depth summaries |

## Before running

1. **Install tools** — see `envs/setup_angsd_env.md`:
   ```bash
   conda activate rumex-snp
   ```

2. **Index your reference** (one-time):
   ```bash
   bwa index data/raw/reference/genome.fasta
   samtools faidx data/raw/reference/genome.fasta
   ```

3. **Prepare the sample name list** at `data/sample_info/sample_names.txt`:
   ```
   sample001
   sample002
   sample003
   ...
   ```

## Running

```bash
# Step 1: Map reads → sorted BAMs (this is the slow step; consider HPC/SLURM)
bash scripts/01_snp_calling/01a_bam_generation.sh

# Step 2: Generate the BAM list for ANGSD
ls results/snp_calling/bam/*.sorted.bam > data/sample_info/bam_list.txt

# Step 3: Call SNPs with ANGSD
bash scripts/01_snp_calling/01b_angsd_snp_calling.sh
```

## Parameters used in the published analysis

These are documented at the top of each script under `USER CONFIGURATION`. Key settings:

- **BWA-MEM:** `-k 19 -B 3 -t 50`
- **ANGSD:** `-GL 2 -minMapQ 25 -minQ 20 -minInd 0.5 -minMaf 0.01 -SNP_pval 1e-6 -maxDepth 1000`

Adjust these for your own dataset by editing the USER CONFIGURATION block at the top of each script.
