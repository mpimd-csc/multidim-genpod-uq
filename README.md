Galerkin POD for PDEs with Uncertainties
---

This is the code of the numerical experiments in our paper

> Benner, Heiland (2020): *Space and Chaos-expansion Galerkin POD Low-order
> Discretization of PDEs for Uncertainty Quantification*

## Installation

Install `dolfin` and `gmesh`.

Then install the python dependencies via

```
pip install dolfin-navier-scipy==1.0.0
pip install sadptprj-riclyap-adi==1.0.0
pip install multidim-galerkin-pod==1.0.2
```
if the installation of `multim-galerkin-pod` fails because of `scikit-sparse`
use `pip install --no-deps multidim-galerkin-pod==1.0.2` instead. 

Then clone/download/unzip this repository

The source are in `gen_pod_uq` and the files for the simulations in `scripts`.

## Rerun the simulations

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

### Meshtest

```
cd scripts
source addpypath.sh
python3 test_the_meshs.py
```

### The small examples of Section 6

```
cd scripts
source addpypath.sh
python3 test_pce_mc.py
```

The file `test_pce_mc.py` has appended the commands that can be processed in
`Mathematica` to compute the reference values.
