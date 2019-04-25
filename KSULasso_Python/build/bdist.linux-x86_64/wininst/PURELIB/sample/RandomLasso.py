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
from sklearn.preprocessing import scale
from sklearn import linear_model

class RandomLasso:
    
    def __init__ (self,independent = None, dependent = None, bootstraps = None, supress = None, cutoff = None):
        self.independent = independent if independent is not None else None
        self.dependent = dependent if dependent is not None else None
        self.sampleCount = self.independent.shape[0] if independent is not None else 1
        self.featureCount = self.independent.shape[1]-1 if independent is not None else 1
        self.bootstraps = bootstraps if bootstraps is not None else int((self.featureCount/self.sampleCount)*80)
        self.supress = self.supress if supress is not None else False
        self.cutoff = self.cutoff if cutoff is not None else None

        
        
        


    def readData(self,pathx, pathy):
        self.independent = pd.read_csv(pathx)
        self.dependent = pd.read_csv(pathy)
        self.sampleCount = self.independent.shape[0]
        self.featureCount = self.independent.shape[1]-1
        self.bootstraps = int((self.featureCount/self.sampleCount)*80)

        return self.independent, self.dependent

    def printData(self):
        print(self.independent.join(self.dependent))

        
    def procedure1(self, bootstraps, algorithm1):
        all_coeff = pd.DataFrame(np.zeros((bootstraps,self.featureCount)))
        for num in range(bootstraps):



            listFeatures = np.random.choice(range(0, self.featureCount),self.sampleCount,replace=False)
            listSamples = np.random.choice(range(0, self.sampleCount),self.sampleCount,replace=True)

            randomIndependent = self.independent.iloc[listSamples, listFeatures]
            randomDependent = self.dependent.iloc[listSamples]

            #math

            randomDependentmean = randomDependent.mean()
            randomDependentscaled = randomDependent - randomDependentmean

            independentmean = pd.DataFrame(randomIndependent.mean())
            independentScale = pd.DataFrame(scale(randomIndependent),index=randomIndependent.index, columns=randomIndependent.columns)
            independentStd = np.std(randomIndependent)

            independantScale2 = pd.DataFrame(scale(independentScale),index=independentScale.index, columns=independentScale.columns)

            algorithm1.fit(independantScale2,randomDependentscaled)

            std_coeff = algorithm1.coef_ / independentStd
            std_coeff = np.array(std_coeff,dtype=np.float32)
            all_coeff.iloc[num,listFeatures] = std_coeff

        
        summedCoeff = all_coeff.sum(axis=0,skipna=True)
        summedCoeff = np.abs(summedCoeff)
        summedCoeff = np.array(summedCoeff)
        summedCoeff /= summedCoeff.sum().astype(float)
        
        return summedCoeff


    def procedure2(self, bootstraps, algorithm2, Probabilities):
        

        
        all_coeff = pd.DataFrame(np.zeros((bootstraps,self.featureCount)))
        for num in range(bootstraps):



            listFeatures = np.random.choice(range(0, self.featureCount),self.sampleCount,replace=False,p=Probabilities)
            listSamples = np.random.choice(range(0, self.sampleCount),self.sampleCount,replace=True)

            randomIndependent = self.independent.iloc[listSamples, listFeatures]
            randomDependent = self.dependent.iloc[listSamples]

            #math

            randomDependentmean = randomDependent.mean()
            randomDependentscaled = randomDependent - randomDependentmean

            independentmean = pd.DataFrame(randomIndependent.mean())
            independentScale = pd.DataFrame(scale(randomIndependent),index=randomIndependent.index, columns=randomIndependent.columns)
            independentStd = np.std(randomIndependent)

            independantScale2 = pd.DataFrame(scale(independentScale),index=independentScale.index, columns=independentScale.columns)

            algorithm2.fit(independantScale2,randomDependentscaled)

            std_coeff = algorithm2.coef_ / independentStd
            std_coeff = np.array(std_coeff,dtype=np.float32)
            all_coeff.iloc[num,listFeatures] = std_coeff


        summedCoeff = all_coeff.sum(axis=0,skipna=True)
        summedCoeff = np.abs(summedCoeff)
        
        summedCoeff = summedCoeff/bootstraps

        return summedCoeff


    def fit(self,algorithm1,algorithm2,x=None,y=None,bootstraps=None):
        self.independent = x if x is not None else None
        self.dependent = y if y is not None else None
        bootstraps = bootstraps if bootstraps is not None else self.bootstraps
        
        probabilities = self.procedure1(bootstraps,algorithm1)

        values = self.procedure2(bootstraps,algorithm2,probabilities)

        return values




        











        

        