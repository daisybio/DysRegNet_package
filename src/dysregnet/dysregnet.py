
import pandas as pd

from . import functions 




# Main code

class run(object):

    
    

        def __init__(self, 
                     expression_data,
                     GRN,
                     meta, 
                     conCol='condition', 
                     CatCov=[],  
                     ConCov=[],
                     zscoring=False,
                     bonferroni_alpha= 1e-2,
                     R2_threshold=None,
                     normaltest=False,
                     normaltest_alpha=1e-3,
                     direction_condition=True):

                    """
                    Raw data processing for further analysis

                    expression_data: a pandas DataFrame (rows=samples, columns=genes)
                             Gene expression matrix with the format: patients as rows (first column - patients/samples ids), and genes as columns.
                             Patients/samples ids must match the ones in meta DataFrame.
                             Gene names or ids must match the ones in GRN DataFrame

                    GRN: a pandas DataFrame 
                          Gene Regulatory Network (GRN) with two columns in the following order ['TF', 'target'].

                    meta: a pandas DataFrame  
                            Meta data: First column should contain patients/samples ids and other column for covariates/condition. 
                            Please make sure to have condition column in the meta DataFrame with 0 as control and 1 as the condition.
                            Specify the condition Column name in "conCol". 
                            Optionally :
                                Specify categorical variable columns in the parameter CatCov.
                                Specify continuous variable columns in the parameter ConCov.

                    conCol: str, default=='condition'
                            Column name for the condition in the meta data. Should be provided in case of desing=="two".


                    CatCov: List of strings.
                            List of categorical variable names. They should match the name of their columns in meta Dataframe.

                    ConCov: List of strings.
                            List of continuous covariates. They should match the name of their columns in meta Dataframe.


                    zscoring: boolean, default: False 
                         zscoring of expression data (if needed).

                    bonferroni_alpha: Float
                            P value threshold for multiple testing correction

                    normaltest: Bool
                            If True. Run a normality test for residuals "scipy.stats.normaltest". If residuals are not normal, the edge will not be considered in the analysis. 

                    normaltest_alpha: Float
                         normaltest p value threshold.

                    R2_threshold: float from 0 to 1 (optional)

                        Coefficient of determination threshold for every edge in GRN. If the R2 is less that threshold, the edge will not be considered in the analysis. 

                    direction_condition: Bool
                         If True: only include dysregulation that are relevalant for the interactions: down regulation of an activation or up regulation of a supressions. Please check the paper for more details.

                          """



                    # init method 
                    self.conCol=conCol
                    self.CatCov=CatCov
                    self.ConCov=ConCov
                    self.zscoring=zscoring
                    self.bonferroni_alpha=bonferroni_alpha
                    self.R2_threshold=R2_threshold
                    self.normaltest=normaltest
                    self.normaltest_alpha=normaltest_alpha
                    self.direction_condition=direction_condition


                    # quality check of parameters

                    # set sample as indexes
                    meta = meta.set_index(meta.columns[0])
                    expression_data = expression_data.set_index(expression_data.columns[0])


                    # check sample ids
                    samples=[ s for s in  list(meta.index) if s in list(expression_data.index) ]
                    if not samples:
                          raise ValueError("Sample columns are not found or the ids don't match. Please make sure that the first column in 'expression_data' and 'meta' are both sample ids.")

                    self.meta=meta.loc[samples]
                    self.expression_data=expression_data.loc[samples]



                    #check condition column
                    if self.conCol not in self.meta.columns:
                            raise ValueError(" Invalid conCol value. Could not find the column '%s' in meta DataFrame" % self.conCol)

                    if set(self.meta[conCol].unique())!={0,1}:
                            raise ValueError(" Invalid values in '%s' column in meta DataFrame. Please make sure to have condition column in the meta DataFrame with 0 as control and 1 as the condition (int)." % self.conCol)

                    # split sample ids (cases and control)
                    self.control= list( self.meta[self.meta[conCol]==0].index )
                    self.case= list( self.meta[self.meta[conCol]==1].index )



                    # check GRN and gene ids

                    GRN_genes=list(set(GRN.iloc[:,0].values.tolist() + GRN.iloc[:,1].values.tolist()))
                    GRN_genes=[g for g in GRN_genes if g in expression_data.columns]

                    if not GRN_genes:
                         raise ValueError('Gene id or name in GRN DataFrame do not match the ones in expression_data DataFrame')

                    self.expression_data=self.expression_data[GRN_genes]
                    self.GRN=GRN[GRN.iloc[:,0].isin(GRN_genes) ]
                    self.GRN=self.GRN[ self.GRN.iloc[:,1].isin(GRN_genes) ].drop_duplicates()


                    self.cov_df,self.expr, self.control, self.case = functions.process_data(self)


                    self.results, self.model_stats = functions.dyregnet_model(self)
                

                
                
        def get_results(self):
            return  self.results

        
        
        
        def get_results_binary(self):
                res_binary=self.results.copy()
                res_binary = res_binary.where(res_binary==0, other=1)
                
                return res_binary

    
        def get_model_stats(self):
            return self.model_stats

        
        
        

            
            
            
            
