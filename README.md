# Private repo for the analysis of Fks hotspots DMS

## gyōza

The first step was to analyze the DMS dataset with [gyōza](https://github.com/durr1602/gyoza) (at this time `--tag 37652ad`)

Single mutants of Hotspots 1, 2 and 3 were analyzed using [this config](./config/config.yaml).

For simplicity, Fks homologous hotspots were analyzed **separately** using the `provided` mode of gyōza. This mode required that expected sequences be provided in a specific format. To generate the proper input (silent mutants from the first dataset + orthologous hotspot sequences), we ran [this custom script](./pre/scripts/generate_gyoza_input_orthologs.py).

For all gyōza analyses, read count threshold at T0 was set at `10`, which means only variants (amino acid sequences) which had more than `10` reads in **all** T0 replicates were used to calculate an average selection coefficient across replicates.

## Requirements

The following steps need to be run in a Python 3.13 virtual environment. With `uv`:

```
uv venv fks_gyoza --python 3.13
source fks_gyoza/bin/activate
uv pip install -r post/requirements.txt
```

## Analysis

### How to reproduce analyses
Simply run [this notebook](./post/notebooks/driver.ipynb)

### Classification of scores to get labels
A Gaussian Mixture Model (GMM) was [trained](./post/notebooks/train_GMM.ipynb) and used to [classify fitness scores](./post/notebooks/classify_gyoza_data.ipynb) calculated by gyōza to obtain labels reflecting mutational effects (deleterious, WT-like, intermediary, resistant, etc).

Thresholds were set to resolve overlaps between the different gaussians predicted.

Classification was collapsed (deleterious/non-deleterious in the control condition, sensitive/resistant upon selection).

### Classification of scores to get labels
DMS data were [compared to growth data](./post/notebooks/20240129_validations_test3.ipynb) from reconstructed mutants to perform a linear regression and infer a score for a few single mutants in FKS1-HS1 missing from the DMS dataset.

### Plot heatmaps
[This notebook](./post/notebooks/heatmaps.ipynb) was used to plot heatmaps of fitness scores.

### Generate defattr files for ChimeraX
[This notebook](./post/notebooks/defattr.ipynb) was used to generate the files used in ChimeraX to plot onto the protein structures the proportion of resistant mutants per position for each hotspot.