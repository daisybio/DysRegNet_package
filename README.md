# DysRegNet package



## Installation


To install the package from git:

`git clone https://github.com/biomedbigdata/DysRegNet_package.git  && cd DysRegNet_package`

`python setup.py install`



## Data input

The inputs of the  package are the following Pandas DataFrame object:


- expression_data  - Gene expression matrix with the format: patients as rows (first column - patients/samples ids), and genes as columns.
- GRN - Gene Regulatory Network (GRN) with two columns in the following order ['TF', 'target'].
- meta -  Metadata with the first column containing patients/samples ids and other columns for the condition and the covariates.


The patients id or samples ids must be the same in the "expression_data" and  "meta". Additionally, gene names or ids must match the ones in the "GRN" DataFrame. 

In the condition column of the meta DataFrame, the control samples should be encoded as 0 and case samples as 1.

GRN network should be provided a prior, You can either use an experimental validated GRN or learn it from control samples, we recommend using software like [arboreto](https://github.com/aertslab/arboreto), since you can use its output directly to DysRegNet.





## Parameters 


Additionally, you can provide the following parameters:


            
- conCol: Column name for the condition in the meta DataFrame.

- CatCov: List of categorical variable names. They should match the name of their columns in the meta Dataframe.

- ConCov: List of continuous covariates. They should match the name of their columns in the meta Dataframe.

- zscoring: Boolean, default: True. zscoring of expression data (if needed).

- bonferroni_alpha:P-value threshold for multiple testing correction

- normaltest: Boolean. If True, Run a normality test for residuals "scipy.stats.normaltest". If residuals are not normal, the edge will not be considered in the analysis. 

- normaltest_alpha: p-value threshold for normaltest (if True).

- R2_threshold: R-squared (R2) threshold from 0 to 1 (optional).  If the fit is weaker, the edge will not be considered in the analysis. 

- direction_condition: Boolean. If True: only include dysregulation that are relevalant for the interactions (down regulation of an activation or up regulation of a supressions). Please check the paper for more details.


## Get Started


Please note, that the functions are annotated with dockstrings for more details.

Import the package and pandas:


```python
import dysregnet
import pandas as pd
```



Define the confounding variables or the design matrix 

```python
# The condition column
conCol='condition'

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
                   conCol=conCol
                   CatCov=CatCov,
                   ConCov=ConCov,
                   direction_condition=True,
                   normaltest=True,
                   R2_threshold=.2 )

# results table
data.get_results()

# or a binary result

data.get_results_binary()

```
## The output

The package output a DataFrame that represents patient-specific dysregulated edges. The columns represent edges and the rows patient ids. 

In the result table, a value of 0 means that the edge is not significantly dysregulated (different from control samples). Otherwise, the z-score is reported, with a positive in case of activation and a negative sign in case of repression (different than the sign of the residual). 

The method "get_results_binary()", outputs binarized dysregulations instead of z-scores. 


## Example

A simple example for running DysRegNet:
([Notebook](https://github.com/biomedbigdata/DysRegNet_package/blob/main/test.ipynb)/[Google Colab](https://colab.research.google.com/github/biomedbigdata/DysRegNet_package/blob/main/test.ipynb)).




## Cite
