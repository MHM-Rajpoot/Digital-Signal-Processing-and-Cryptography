import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Random vectors
X = pd.read_csv('random_vectors.csv').iloc[:,1:]
plt.figure()
plt.plot(X)
plt.title("Random Vectors")
plt.xlabel("Sample Index")
plt.ylabel("Amplitude")
plt.grid(True)
plt.show()

# Covariance matrix heatmap
cov = pd.read_csv('covariance_matrix.csv', index_col=0)
plt.figure()
plt.imshow(cov.values, cmap='viridis', interpolation='none')
plt.colorbar()
plt.title("Covariance Matrix")
plt.show()

# PCA transformed coefficients (first 10 vectors)
Y_pca = pd.read_csv('pca_transformed_vectors.csv').iloc[:,1:]
plt.figure()
plt.plot(Y_pca.iloc[0:10,:])
plt.title("PCA Transformed Vectors")
plt.xlabel("Principal Component Index")
plt.ylabel("Amplitude")
plt.grid(True)
plt.show()

# DCT transformed coefficients (first 10 vectors)
Y_dct = pd.read_csv('dct_transformed_vectors.csv').iloc[:,1:]
plt.figure()
plt.plot(Y_dct.iloc[0:10,:])
plt.title("DCT Transformed Vectors")
plt.xlabel("DCT Coefficient Index")
plt.ylabel("Amplitude")
plt.grid(True)
plt.show()