# 2023-ms-workflow

## RNA-Seq workflow for de novo transcriptome

### Quickstart

First, `git clone` the software and `cd` into the folder:

```
git clone https://github.com/bluegenes/2023-ms-workflow.git
cd 2023-ms-workflow
```

Install the `conda` environment that contains essential software for running (Snakemake, python pandas library).

```
mamba env create -f environment.yml
```
> You may need to install the conda/mamba package manager first. Recommend installing via [mambaforge](https://github.com/conda-forge/miniforge#mambaforge). Root not needed.

To run, we activate the conda environment and use `snakemake -n` to show the steps the will be run:

```
conda activate 2023-ms-workflow
snakemake -n
```


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

