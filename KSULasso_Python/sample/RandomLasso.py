#' Performs variable selection and regularization using Random Lasso.
#'
#'
#' independent: independent data.
#' dependent: Dependent data.
#' bootstraps: Number of times random sampling occurs.
#' suppress: Supresses all printing and time estimation.
#' cutoff: A cut off value for reducing insignificant values to zero.

import numpy
import adaptive 
import pandas as pd

class RandomLasso:
    def __init__ (self, independent, dependent, bootstraps, supress, cutoff):
        self.independent = independent # x
        self.dependent = dependent # y
        self.bootstraps = bootstraps
        self.supress = supress
        self.cutoff = cutoff # optional
        
    def RandomLasso(self, independent, dependent, bootstraps, supress, cutoff):
        print("This is the main function")

    def procedure1(self, bootstraps, algorithm1, algorithm2):
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
        while(bootstraps!=0):

    