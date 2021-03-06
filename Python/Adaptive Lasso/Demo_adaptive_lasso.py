
"""Example of adaptive Lasso to produce event sparser solutions
Adaptive lasso consists in computing many Lasso with feature
reweighting. It's also known as iterated L1.
"""

import numpy as np
import pandas as pd
from sklearn.datasets import make_regression
from sklearn.linear_model import Lasso
from sklearn.preprocessing import LabelEncoder

#data = pd.read_csv("/Users/danielchoo/PythonProjects/Adaptive Lasso/test2_x.csv")
#data.drop(['Unnamed: 0'], axis=1, inplace = True)
#data.head(10000)

#def preProcessDataset():
#    data = pd.read_csv('test2_x.csv')
#    data.drop(['Unnamed: 0'], axis=1, inplace = True)
#    data.head(10000)


X, y, coef = make_regression(n_samples=306, n_features=8000, n_informative=50,
                             noise=0.1, shuffle=True, coef=True, random_state=42)

X /= np.sum(X ** 2, axis=0)  # scale features
alpha = 0.1

# g = lambda w: np.sqrt(np.abs(w))
def g(w): 
    return np.sqrt(np.abs(w))


def gprime(w): 
    return 1. / (2. * np.sqrt(np.abs(w)) + np.finfo(float).eps)

# Or another option:
# ll = 0.01
# g = lambda w: np.log(ll + np.abs(w))
# gprime = lambda w: 1. / (ll + np.abs(w))


n_samples, n_features = X.shape


def p_obj(w): 
    return 1. / (2 * n_samples) * np.sum((y - np.dot(X, w)) ** 2) \
    + alpha * np.sum(g(w))


weights = np.ones(n_features)
n_lasso_iterations = 5

for k in range(n_lasso_iterations):
    X_w = X / weights[np.newaxis, :]
    clf = Lasso(alpha=alpha, fit_intercept=False)
    clf.fit(X_w, y)
    coef_ = clf.coef_ / weights
    weights = gprime(coef_)
    print(p_obj(coef_))  # should go down

print(np.mean((clf.coef_ != 0.0) == (coef != 0.0)))

