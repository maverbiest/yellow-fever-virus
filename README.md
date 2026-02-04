# Nextclade build for Yellow Fever Virus (full genome)

## Instructions
To run this build yourself, you will first need to install some dependencies.
We provide a `.yaml` file you can use to set up a conda/mamba environment, e.g.:
```bash
micromamba env create -f environment.yaml
micromamba activate yellow-fever-virus
```

Next, create the directory structure expected by the workflow:
```bash
mkdir -p data/translations results/ out-dataset/
```

Following this, you can run the workflow using `snakemake --cores <n> all`, where you substitute `<n>` with the number of cores you want to allocate the run.

The resulting file `auspice/yellow-fever.json` can be viewed using auspice. The individual files making up the Nextclade dataset will be stored in `out-dataset`.
