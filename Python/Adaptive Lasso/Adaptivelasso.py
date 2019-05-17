# Implemented Adaptive Lasso 

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.linear_model import Lasso
from sklearn.linear_model import base
from sklearn.linear_model import coordinate_descent
from sklearn.utils import check_array, check_X_y, deprecated
from sklearn.utils import check_array, check_X_y, column_or_1d, deprecated
from sklearn.utils.validation import check_random_state

class AdaptiveLasso:
    '''
    The multi-step adaptive Lasso iteratively solves Lasso estimates with 
    penalty weights applied to the regularization of the coefficients.
    The optimization objective for the adaptive Lasso is::
    
        1/n * ||y - X Beta||^2_2 + alpha * w |Beta|_1
    
    Where w is a weight vector calculated in the previous stage by::
    
        w_j = alpha/(|Beta_j|^gamma + eps)
    '''
    def __init__(self, n_lasso_iterations = 5, gamma = 1.0, alpha=1.0,
                 eps = None, fit_intercept=True, normalize=False,
                 precompute=False, copy_X=True, max_iter=1000,
                 tol=1e-4, positive=False,
                 random_state=None, selection='cyclic'):
        super(AdaptiveLasso, self).__init__(
            alpha=alpha, fit_intercept=fit_intercept,
            normalize=normalize, precompute=precompute, copy_X=copy_X,
            max_iter=max_iter, tol=tol, warm_start=False,
            positive=positive, random_state=random_state,
            selection=selection)
        self.n_lasso_iterations = n_lasso_iterations
        self.gamma = gamma
        self.eps = eps # optional parameter
        self.alpha = alpha
        self.copy_X = copy_X
        self.fit_intercept = fit_intercept
        
    def fit(self, X, y):
        if self.eps is None:
            eps_ = np.finfo(float).eps
        else:
            eps_ = self.eps
        if self.gamma <= 0:
            raise ValueError('gamma must be positive')
       # else:
       #     self.update_ = lambda w: np.power(np.abs(w) + eps_, -self.gamma)

        alphas = column_or_1d(np.atleast_1d(self.alpha))


        if alphas.shape[0] != 1 and alphas.shape[0] != self.n_lasso_iterations:
            raise ValueError("alpha must be a float or an array of length=%s" % repr(self.n_lasso_iterations))
        if alphas.shape[0] != self.n_lasso_iterations:
            alphas = column_or_1d(np.repeat(np.atleast_1d(self.alpha),  self.n_lasso_iterations))

        X, y = check_X_y(X, y, accept_sparse='csc', dtype=np.float64,
                 order='F', copy=self.copy_X and self.fit_intercept,
                 multi_output=True, y_numeric=True)

       # X, y, X_mean, y_mean, X_std, precompute, Xy = \
       #     _pre_fit(X, y, None, self.precompute, self.normalize,
       #              self.fit_intercept, copy=True)

        n_samples, n_features = X.shape
        weights = np.ones(n_features)

        for k in range(self.n_lasso_iterations):
            # X_w = X / weights[np.newaxis, :]
            X_w = np.divide(X,weights[np.newaxis, :])
            self.alpha = alphas[k]
            super(AdaptiveLasso, self.fit(X_w, y))
            self.coef_ /= weights
            #weights = self.update_(self.coef_)
            weights = np.power(np.abs(self.coef_) + eps_, -self.gamma) 

        return self
    