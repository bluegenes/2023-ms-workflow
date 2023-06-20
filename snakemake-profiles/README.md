# Snakemake Profiles

To use:

```
mkdir -p ~/.config/snakemake
```
# copy profile folders to `~/.config/snakemake`
```
cd snakemake-profiles
cp -r * ~/.config/snakemake
```
> each profile folder has a `config.yaml` file that contains configuration options that are passed to snakemake

Use the profiles when running snakemake
```
snakemake --profile default # default interactive session
snakemake --profile debug # debug/testing (do not rerun jobs if they fail)
snakemake --profile slurm # submit jobs via sbatch
```
