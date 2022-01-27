# DysRegNet package



## Installation


To install the package from git:

`git clone https://github.com/biomedbigdata/DysRegNet_package.git  && cd DysRegNet_package`

`pip install .`



## Data input

The inputs of the  package are the following Pandas DataFrame object:


- expression_data  - Gene expression matrix with the format: patients as rows (first column - patients/samples ids), and genes as columns.
- GRN - Gene Regulatory Network (GRN) with two columns in the following order ['TF', 'target']
- meta -  Metadata with the first column containing patients/samples ids and other columns for covariates.


Note that: the patients' or samples ids must be the same in the "expression_data" and  "meta". Additionally, gene names or ids must match the ones in the "GRN" DataFrame. GRN network should be provided a prior, we recommend using the software like  arboreto ([arboreto](https://github.com/aertslab/arboreto), and feed use its output to DysRegNet.

Additionally, you can provide the following parameters:


            
- conCol: str, default=='condition'
        Column name for the condition in the metadata. Should be provided in case of desing=="two".


- CatCov: List of strings.
        List of categorical variable names. They should match the name of their columns in the meta Dataframe.

- ConCov: List of strings.
        List of continuous covariates. They should match the name of their columns in the meta Dataframe.


- zscoring: boolean, default: True 
     zscoring of expression data (if needed).

- bonferroni_alpha: Float
        P-value threshold for multiple testing correction

- normaltest: Bool
        If True. Run a normality test for residuals "scipy.stats.normaltest". If residuals are not normal, the edge will not be considered in the analysis. 

- normaltest_alpha: Float
     normaltest p value threshold.

- R2_threshold: float from 0 to 1 (optional)




## Get Started


Please note, that the functions are annotated with dockstrings for more details.

Import the package and pandas:


```python
import dysregnet
import pandas as pd
```



Define the confounding variables or the design matrix 

```python
# categorical variable columns in meta dataframe.
# these columns will be transformed to variables for regression 
CatCov=['race','gender']  

# continuous variable columns in meta dataframe.
ConCov=['birth_days_to']
```


Run DysRegNet

```python
data=dysregnet.run(expression_data=expr,
                   meta=meta, 
                   GRN=grn,
                   conCol='condition'
                   CatCov=CatCov,
                   ConCov=ConCov,
                   normaltest=True,
                   R2_threshold=.2,
                   conCol='sample')

# results table
data.get_results()

```


## Cite
