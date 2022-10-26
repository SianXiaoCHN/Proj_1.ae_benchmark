# Proj_1.ae_benchmark
Autoencoder benchmark


## About openmm
The openmm is installed with cuda 10.0.
```bash
conda install -c conda-forge openmm cudatoolkit=10.0
```

## About environment
After change the environment, env.txt must be updated.

Export the conda env:
```bash
conda list -e > <environment-name>.txt
```

Install from req.txt:
```bash
conda create -n <environment-name> --file <environment-name>.txt
```