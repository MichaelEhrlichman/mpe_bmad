#!/usr/bin/env python3

# Adapted from https://scikit-learn.org/stable/auto_examples/gaussian_process/plot_gpc.html#sphx-glr-auto-examples-gaussian-process-plot-gpc-py

import numpy as np

from matplotlib import pyplot as plt

from sklearn.metrics import accuracy_score, log_loss
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF

# Generate data
ma_data = np.loadtxt('along_0p.dat')
X = ma_data[:,0][:,np.newaxis]
y = ma_data[:,1]
y = np.where(y < 0, 1, 0)

# Specify Gaussian Processes with fixed and optimized hyperparameters
gp_opt = GaussianProcessClassifier(kernel=1.0 * RBF(length_scale=0.05))
gp_opt.fit(X[:], y[:])

print("Log Marginal Likelihood (optimized): %.3f"%gp_opt.log_marginal_likelihood(gp_opt.kernel_.theta))

print("Accuracy: %.3f "%accuracy_score(y[:], gp_opt.predict(X[:])))
print("Log-loss: %.3f "%log_loss(y[:], gp_opt.predict_proba(X[:])[:, 1]))

# Plot posteriors
plt.figure()
plt.scatter(X[:, 0], y[:], c='k', label="Train data",
            edgecolors=(0, 0, 0))
#X_ = np.linspace(0, 5, 100)
X_ = np.linspace(0, 0.05, 100)
plt.plot(X_, gp_opt.predict_proba(X_[:, np.newaxis])[:, 1], 'b', label="Optimized kernel: %s" % gp_opt.kernel_)
plt.xlabel("Feature")
plt.ylabel("Class 1 probability")
#plt.xlim(0, 5)
plt.xlim(0, 0.05)
plt.ylim(-0.25, 1.5)
plt.legend(loc="best")

plt.show()
