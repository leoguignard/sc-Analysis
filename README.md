# sc-Analysis

This is just a notebook to read and analyse single-cell "omics" data based on the very good [scanpy](https://scanpy.readthedocs.io/en/latest/installation.html) library.

## Description of the repository

- utilities: folder containing some useful functions
- setup.py: Installation script
- README.md: This file
- Data-Analysis.ipynb: an jupyter notebook with which you might be able to analyse single cell data

## Basic usage

You want to run the notebook in a version of python >3.6 and then follow the instruction within the notebook. One way to run the notebook is to run the following command, from a terminal when position within the folder

```shell
jupyter notebook
```

## Dependencies

Some dependecies are requiered:

- general python dependecies:
  - numpy, matplotlib
- single cell analysis dependencies:
  - [scanpy](https://scanpy.readthedocs.io/en/latest/installation.html)
  - [anndata](https://anndata.readthedocs.io/en/stable/index.html)
  - [leidenalg](https://github.com/vtraag/leidenalg)
  - [igraph](https://igraph.org/python/doc/igraph-module.html)

## Quick install

To quickly install the necessary dependencies you can run the following command from a terminal:

```shell
pip install .
```

It is possible that you lack the rights to install the necessary libraries. In that case I would advise to use virtual environements (using conda for example) or appending the option `--user` at the end of either of the command:

```shell
pip install . --user
```


