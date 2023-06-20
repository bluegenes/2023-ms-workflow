# 2023-ms-workflow

## RNA-Seq workflow for de novo transcriptome

### Quickstart

1. First, `git clone` the software and `cd` into the folder:

```
git clone https://github.com/bluegenes/2023-ms-workflow.git
cd 2023-ms-workflow
```

2. Install the `conda` environment that contains essential software for running (Snakemake, python pandas library).

```
mamba env create -f environment.yml
```
> You may need to install the conda/mamba package manager first. Recommend installing via [mambaforge](https://github.com/conda-forge/miniforge#mambaforge). Root not needed.

3. Activate the conda environment 
```
conda activate 2023-ms-workflow
```

4. Make sure the reads are in the right spot (the 'inputs' folder). Two options:
- Copy or link reads to the 'inputs' folder
- Use the 'download_reads' workflow target to download:
```
snakemake download_reads -j1 -n
```
> - `-j1` number of simultaneous jobs
> - `-n` specifies a 'dryrun'; remove to actually execute

Do a full workflow 'dryrun' to see what steps will be run:
```
snakemake -n -j 1
```
> To actually execute the steps, remove the `-n` (dryrun) and specify the number of jobs you'd like to run at once with `-j `


### Workflow Steps:
  - Adapter trim (and lightly quality trim) reads
    - `fastp`
  - Assess and visualize raw and trimmed read quality
    - `fastqc`
    - `multiqc`
  - Assemble reads
    - `Trinity`
  - Quantify reads
    - `Salmon`

To be added:
  - Annotation
    - `Trinotate`
    - `dammit` (v2)
  - Differential Expression
    - `DESeq2`


## Resource requirements

[Trinity Resource Requirements](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Computing-Requirements)
