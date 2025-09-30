# SSIM

[Second School of Integrated Metabolomics](https://escoladeproteomica2.brprot.com.br/)

## Local installation

Install conda

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh 
```

Create a dedicated conda environment (`env` folder for each notebook) and activate

```
conda env create -f env/environment_pyopenms.yml
conda activate pyopenms
pip install jupyter
jupyter notebook
```

Alternatively, create dedicated environment

```
conda create -n pyopenms python=3.12
conda activate pyopenms
pip install pyopenms==3.2
pip install ipykernel
# Use multiple kernels on your local jupyter install
python -m ipykernel install --user --name pyopenms --display-name pyopenms
# update env file after installing new libraries
conda env export | grep -v "^prefix: " > environment.yml 
```

## Index

- [GNPS2 data inspection](./notebooks/gnps2_data_inspection.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](http://colab.research.google.com/github/computational-chemical-biology/SSIM/blob/master/notebooks/gnps2_data_inspection.ipynb)
- [pyOpenMS](./notebooks/pyopenms-api.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](http://colab.research.google.com/github/computational-chemical-biology/SSIM/blob/master/notebooks/pyopenms-api.ipynb)
- [Compare pyOpenMS and MZmine](./notebooks/comp_mzmine_pyopenms.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](http://colab.research.google.com/github/computational-chemical-biology/SSIM/blob/master/notebooks/comp_mzmine_pyopenms.ipynb)
- [Format database for ChemWalker](./notebooks/formatdb.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](http://colab.research.google.com/github/computational-chemical-biology/SSIM/blob/master/notebooks/formatdb.ipynb)
- [MassQL](./notebooks/massql.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](http://colab.research.google.com/github/computational-chemical-biology/SSIM/blob/master/notebooks/massql.ipynb)
- [GNPS2 API](./notebooks/gnps-api.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](http://colab.research.google.com/github/computational-chemical-biology/SSIM/blob/master/notebooks/gnps-api.ipynb)

