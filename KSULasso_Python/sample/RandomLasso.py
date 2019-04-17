#' Performs variable selection and regularization using Random Lasso.
#'
#'
#' independent: independent data.
#' dependent: Dependent data.
#' bootstraps: Number of times random sampling occurs.
#' suppress: Supresses all printing and time estimation.
#' cutoff: A cut off value for reducing insignificant values to zero.

import numpy as np
# import adaptive 
import pandas as pd
from sklearn import linear_model

class RandomLasso:
    
    def __init__ (self,algorithm1, independent = None, dependent = None, bootstraps = None, supress = None, cutoff = None):
        self.algorithm1 = algorithm1
        self.independent = pd.read_csv(independent) if independent is not None else None
        self.dependent = pd.read_csv(dependent) if dependent is not None else None
        self.sampleCount = self.independent.shape[0]
        self.featureCount = self.independent.shape[1]-1
        self.bootstraps = bootstraps if bootstraps is not None else int(self.featureCount/self.sampleCount)
        self.supress = self.supress if supress is not None else False
        self.cutoff = self.cutoff if cutoff is not None else None
        self.procedure1(self.bootstraps,algorithm1)


    def readData(self,pathx, pathy):
        self.independent = pd.read_csv(pathx)
        self.dependent = pd.read_csv(pathy)

    def printData(self):
        print(self.independent.join(self.dependent))

        
    def procedure1(self, bootstraps, algorithm1):
        # For every boot strap we select sample of random feature names, num of features = number of samples
        # Mix up the rows
        #replace = no repetition.
        # get row of random samples and the columns will be random features
        # Dependent 1 row selects the random samples. 
        # Averaging for normalization
        # Centering the independant  
        #STD
        # scale = Centering
        # RETURN THE COEFFECIENTS for 1 random sample test.
        #step 2 you just sample with the results of the last step
        s1 = []
        s2 = []
        s3 = []
        for i in range(bootstraps):
            s1.append(i)
            s2.append(i)
            s3.append(i)
        
        for num in s1:
            s1[num] = self.independent.sample(self.sampleCount,axis=1,replace=False)
            appended = s1[num].join(self.dependent)
            s1[num] = appended.sample(self.sampleCount,axis=0,replace=True)

            s2[num] = s1[num].loc[:, s1[num].columns != 'V1']
            s3[num] = s1[num].loc[:, s1[num].columns == 'V1']
        
        self.x = s2
        self.y = s3

        algorithm1.fit(self.x[0],self.y[0])

        self.proc1Coeff = (algorithm1.coef_)

    def printCoeff(self):
        print(self.proc1Coeff)


        

        