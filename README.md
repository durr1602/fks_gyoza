# Private repo for the analysis of Fks hotspots DMS

## gyōza

The first step was to analyze the DMS dataset with [gyōza](https://github.com/durr1602/gyoza)

This was done by `snakeploy`ing gyōza with tag [`89c40d6`](https://github.com/durr1602/gyoza/commit/89c40d6e0c46e47cd2b936f219e1819f372869ac). Config is in `config/`.

Read count threshold at T0 was set at `10`, which means only variants (amino acid sequences) which had more than `10` reads in **all** T0 replicates were used to calculate an average selection coefficient across replicates.

`snakemake==8.24.0` was used for the analysis.

## Requirements

The following steps need to be run in a Python 3.13 virtual environment. With `uv`:

```
uv venv fks_gyoza --python 3.13
source fks_gyoza/bin/activate
uv pip install -r post/requirements.txt
```

## Analysis

### Classification of scores to get labels
A Gaussian Mixture Model (GMM) was used to [classify fitness scores](./post/notebooks/classify_gyoza_data.ipynb) calculated by gyoza to obtain labels reflecting mutational effects (deleterious, WT-like, intermediary, resistant, etc).

Thresholds were set to resolve overlaps between the different gaussians predicted.

Classification was collapsed (deleterious/non-deleterious in the control condition, sensitive/resistant upon selection).

### Classification of scores to get labels
DMS data were [compared to growth data](./post/notebooks/20240129_validations_test3.ipynb) from reconstructed mutants to perform a linear regression and infer a score for a few mutants missing from the DMS dataset.