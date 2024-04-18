[![PyPI version](https://badge.fury.io/py/dysregnet.svg)](https://badge.fury.io/py/dysregnet)

# DysRegNet package
DysRegNet, is a  method for inferring patient-specific regulatory alterations (dysregulations) from gene expression profiles. DysRegNet uses linear models to account for confounders and residual-derived z-scores to assess significance.
## Installation
To install the package from PyPI please run:
```bash
pip install dysregnet
```

or you can install it from git:
```bash
git clone https://github.com/biomedbigdata/DysRegNet_package.git  && cd DysRegNet_package
python setup.py install
```

## Data input
The inputs of the  package are the following Pandas DataFrame objects:

- expression_data  - Gene expression matrix in the format: patients as rows (first column - patients/samples ids), and genes as columns.
- GRN - Gene Regulatory Network (GRN) with two columns in the following order ['TF', 'target'].
- meta -  Metadata with the first column containing patients/samples ids and other columns for the condition and the covariates.

The patients id or samples ids must be the same in the "expression_data" and  "meta". Additionally, gene names or ids must match the ones in the "GRN" DataFrame. 

In the condition column of the meta DataFrame, the control samples should be encoded as 0 and case samples as 1.

The gene regulatory network should be provided by the user. You can either use an experimental validated GRN or learn it from control samples. We recommend using software like [arboreto](https://github.com/aertslab/arboreto) since you can use its output directly to DysRegNet.

## Parameters 
Additionally, you can provide the following parameters:

- conCol: Column name for the condition in the meta DataFrame.

- CatCov: List of categorical variable names. They should match the name of their columns in the meta Dataframe.

- ConCov: List of continuous covariates. They should match the name of their columns in the meta Dataframe.

- zscoring: If True, DysRegNet will scale the expression of each gene and all continuous confounders based on their mean and standard deviation in the control samples.

- bonferroni_alpha: P-value threshold for multiple testing correction

- normaltest: If True, DysRegNet runs a normality test for residuals "scipy.stats.normaltest". If residuals are not normal, the edge will not be considered in the analysis. 

- normaltest_alpha: P-value threshold for normaltest (if True).

- R2_threshold: R-squared (R2) threshold from 0 to 1 (optional).  If the fit is weaker, the edge will not be considered in the analysis. 

- direction_condition:  If True, DysRegNet will only consider case samples with positive residuals (target gene overexpressed) for models with a negative TF coefficient as potentially dysregulated. Similarly, for positive TF coefficients, only case samples with negative residuals are considered. Please check the paper for more details.

The parameters are also annotated with dockstrings for more details.

## Get Started
Import the package and pandas:
```python
import dysregnet
import pandas as pd
```

Define the confounding variables or the design matrix 
```python
# define condition column (0 indicated control, 1 indicates case)
conCol='condition'

# define categorical confounder columns in meta dataframe 
CatCov=['race','gender']  

# define continuous confounder columns in meta dataframe.
ConCov=['birth_days_to']
```

Run DysRegNet
```python
data=dysregnet.run(expression_data=expr,
                   meta=meta, 
                   GRN=grn,
                   conCol=conCol
                   CatCov=CatCov,
                   ConCov=ConCov,
                   direction_condition=True,
                   normaltest=True,
                   R2_threshold=.2)

# get the patient-specific dysregulate networks
data.get_results()

# or with binary edges
data.get_results_binary()

# get R2 values, coefficients, and coefficient p-values for all models/edges
data.get_model_stats()
```

## The output
The package outputs a data frame that represents patient-specific dysregulated edges. The columns represent edges, and the rows are patient IDs. 

In the result table, a value of 0 means that the edge is not significantly dysregulated (different from control samples). Otherwise, the z-score is reported. 

The method "get_results_binary()" outputs binarized dysregulations instead of z-scores. 

"get_model_stats()" outputs R2 values, coefficients, and coefficient p-values for all models/edges.

## Example

A simple example for running DysRegNet:
([Notebook](https://github.com/biomedbigdata/DysRegNet_package/blob/main/test.ipynb)/[Google Colab](https://colab.research.google.com/github/biomedbigdata/DysRegNet_package/blob/main/test.ipynb)).

You will need to download the demo dataset and extract the files into test dataset/

Link for the demo dataset: https://figshare.com/ndownloader/files/35142652

## Cite
"DysRegNet: Patient-specific and confounder-aware dysregulated network inference"
Johannes Kersting*, Olga Lazareva*, Zakaria Louadi*, David B. Blumenthal, Jan Baumbach, Markus List. bioRxiv 2022.04.29.490015; doi: https://doi.org/10.1101/2022.04.29.490015. * equal first-authors