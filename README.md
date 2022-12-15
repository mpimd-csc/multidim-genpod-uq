Galerkin POD for PDEs with Uncertainties
---

[![DOI](https://zenodo.org/badge/291024430.svg)](https://zenodo.org/badge/latestdoi/291024430)

This is the code of the numerical experiments in our paper

> Benner, Heiland (2020): *Space and Chaos-expansion Galerkin POD Low-order
> Discretization of PDEs for Uncertainty Quantification*

in the second version from July 2022.

## Installation

Install `dolfin` and `gmesh`.

Then clone this repo and install the package with dependencies via

```
pip install -e .  # make sure you use Python 3
```

if the installation of `multim-galerkin-pod` fails because of `scikit-sparse`
try `pip install --no-deps multidim-galerkin-pod==1.1.0` instead. 

The source are in `gen_pod_uq` and the files for the simulations are in `scripts`.

## Rerun the simulations

**NOTE**: For reproduction of the results, use version `1.1.1` of the package to be installed like

```sh
pip install gen-pod-uq==1.1.1
```

from the [`pypi repo`](https://pypi.org/project/gen-pod-uq/)

### Generate the mesh
```
cd mesh
mkdir 3D-mshs
source maketheme-3D.sh
```

### Results of the PCE and POD approximations

To reproduce the results of the manuscript

```
cd scripts
source runitall.sh
```

You may want to comment out some parts.

The raw data of our simulations is provided in the folder `rawdata`. In order to postprocess copy it to the `scripts/cached-data` folder

```sh
# ## caution: this may replace computed data
# cp rawdata/*json scripts/cached-data/
# ## caution: this may replace computed data
```

### Post Processing

```
cd scripts
source postprocess.sh
```

