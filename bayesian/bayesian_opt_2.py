#!/usr/bin/env python3

# Adapted from https://scikit-learn.org/stable/auto_examples/gaussian_process/plot_gpc.html#sphx-glr-auto-examples-gaussian-process-plot-gpc-py

import numpy as np

from matplotlib import pyplot as plt

from sklearn.metrics import accuracy_score, log_loss
#from sklearn.gaussian_process import GaussianProcessClassifier
import sklearn.gaussian_process as gp
from sklearn.gaussian_process.kernels import RBF

def surrogate(model, X):
	return model.predict(X)

def acquisition(X, Xsamples, model):
	yhat = surrogate(model, X)
	best = max(yhat)
	mu = surrigate(model, Xsamples)
	mu = mu[:,0]
	probs = norm.cdf((mu-best) / (std+1e-9))

# Generate data
ma_data = np.loadtxt('along_0p.dat')
X = ma_data[:,0][:,np.newaxis]
y = ma_data[:,1]
y = np.where(y < 0, 1, 0)

# Specify Gaussian Processes with fixed and optimized hyperparameters
gpc_model = gp.GaussianProcessClassifier(kernel=1.0 * RBF(length_scale=0.05))
gpc_model.fit(X[:], y[:])

print("Log Marginal Likelihood (optimized): %.3f"%gpc_model.log_marginal_likelihood(gpc_model.kernel_.theta))
print("Accuracy: %.3f "%accuracy_score(y[:], gpc_model.predict(X[:])))
print("Log-loss: %.3f "%log_loss(y[:], gpc_model.predict_proba(X[:])[:, 1]))
Xtest = np.array(0.02).reshape(-1,1)
print("predict: {}".format(gpc_model.predict(Xtest)))
print("predict proba: {}".format(gpc_model.predict_proba(Xtest)))

# Plot posteriors
plt.figure()
plt.scatter(X[:, 0], y[:], c='k', label="Train data", edgecolors=(0, 0, 0))
#X_ = np.linspace(0, 5, 100)
X_ = np.linspace(0, 0.05, 100)
plt.plot(X_, gpc_model.predict_proba(X_[:, np.newaxis])[:, 1], 'b', label="Optimized kernel: %s" % gpc_model.kernel_)
plt.xlabel("Feature")
plt.ylabel("Class 1 probability")
#plt.xlim(0, 5)
plt.xlim(0, 0.05)
plt.ylim(-0.25, 1.5)
plt.legend(loc="best")

plt.show()
