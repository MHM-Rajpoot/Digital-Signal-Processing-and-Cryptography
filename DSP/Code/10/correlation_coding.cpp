/*
 * correlation_coding_pca_dct.cpp
 *
 * Demonstration of correlation coding:
 *  - Random vectors
 *  - Covariance and correlation
 *  - PCA (Karhunen–Loève Transform)
 *  - Eigen decomposition
 *  - DCT as orthogonal transform coding
 *  - Comparison of PCA vs DCT energy compaction
 *
 * CSV export for MATLAB/Python/GNUplot.
 *
 * Inspired by Dr Markus Kuhn, Computer Laboratory, University of Cambridge.
 *
 * Compile with: g++ -std=c++11 -o correlation_coding_pca_dct correlation_coding_pca_dct.cpp -lm
 */

#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <algorithm>

using namespace std;

/* ===== CSV export ===== */
void export_csv(const string &filename, const vector<vector<double>> &matrix) {
    ofstream file(filename);
    int N = matrix.size();
    int M = matrix[0].size();
    file << "n";
    for(int j=0;j<M;j++) file << ",col" << j;
    file << "\n";
    for(int i=0;i<N;i++){
        file << i;
        for(int j=0;j<M;j++)
            file << "," << matrix[i][j];
        file << "\n";
    }
    file.close();
    cout << "Exported " << filename << endl;
}

/* ===== Generate correlated random vectors ===== */
vector<vector<double>> generate_random_vectors(int N, int M, double rho=0.8) {
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> dis(0,1);
    vector<vector<double>> X(N, vector<double>(M,0.0));
    for(int n=0;n<N;n++){
        double common = dis(gen);
        for(int m=0;m<M;m++){
            double independent = dis(gen);
            X[n][m] = rho*common + sqrt(1.0 - rho*rho)*independent;
        }
    }
    return X;
}

/* ===== Compute mean of columns ===== */
vector<double> column_mean(const vector<vector<double>> &X){
    int N = X.size();
    int M = X[0].size();
    vector<double> mean(M,0.0);
    for(int m=0;m<M;m++){
        for(int n=0;n<N;n++) mean[m] += X[n][m];
        mean[m] /= N;
    }
    return mean;
}

/* ===== Compute covariance matrix ===== */
vector<vector<double>> covariance_matrix(const vector<vector<double>> &X){
    int N = X.size();
    int M = X[0].size();
    vector<double> mean = column_mean(X);
    vector<vector<double>> cov(M, vector<double>(M,0.0));
    for(int i=0;i<M;i++){
        for(int j=0;j<M;j++){
            for(int n=0;n<N;n++)
                cov[i][j] += (X[n][i]-mean[i])*(X[n][j]-mean[j]);
            cov[i][j] /= (N-1);
        }
    }
    return cov;
}

/* ===== PCA via Eigen decomposition (Jacobi) ===== */
void jacobi_eigen(const vector<vector<double>> &A, vector<double> &eigenvalues, vector<vector<double>> &eigenvectors, int max_iter=100) {
    int N = A.size();
    vector<vector<double>> V(N, vector<double>(N,0.0));
    eigenvectors = V;
    eigenvalues.resize(N,0.0);
    vector<vector<double>> D = A;

    for(int i=0;i<N;i++) eigenvectors[i][i]=1.0;

    for(int iter=0;iter<max_iter;iter++){
        double max_off=0.0; int p=0,q=1;
        for(int i=0;i<N;i++)
            for(int j=i+1;j<N;j++)
                if(fabs(D[i][j])>max_off){ max_off=fabs(D[i][j]); p=i; q=j;}
        if(max_off<1e-10) break;

        double theta = 0.5*atan2(2*D[p][q], D[q][q]-D[p][p]);
        double c = cos(theta), s = sin(theta);

        double Dpp = c*c*D[p][p]-2*s*c*D[p][q]+s*s*D[q][q];
        double Dqq = s*s*D[p][p]+2*s*c*D[p][q]+c*c*D[q][q];
        double Dpq = 0.0;

        for(int j=0;j<N;j++){
            if(j!=p && j!=q){
                double Djp = c*D[j][p]-s*D[j][q];
                double Djq = s*D[j][p]+c*D[j][q];
                D[j][p]=D[p][j]=Djp;
                D[j][q]=D[q][j]=Djq;
            }
        }
        D[p][p]=Dpp; D[q][q]=Dqq; D[p][q]=D[q][p]=Dpq;

        for(int i=0;i<N;i++){
            double vip = c*eigenvectors[i][p]-s*eigenvectors[i][q];
            double viq = s*eigenvectors[i][p]+c*eigenvectors[i][q];
            eigenvectors[i][p]=vip;
            eigenvectors[i][q]=viq;
        }
    }

    for(int i=0;i<N;i++) eigenvalues[i]=D[i][i];
}

/* ===== Apply linear transform to data ===== */
vector<vector<double>> apply_transform(const vector<vector<double>> &X, const vector<vector<double>> &T){
    int N = X.size();
    int M = X[0].size();
    vector<vector<double>> Y(N, vector<double>(M,0.0));
    for(int n=0;n<N;n++)
        for(int k=0;k<M;k++)
            for(int m=0;m<M;m++)
                Y[n][k] += X[n][m]*T[m][k];
    return Y;
}

/* ===== DCT Type-II transform matrix ===== */
vector<vector<double>> dct_matrix(int M){
    vector<vector<double>> T(M, vector<double>(M,0.0));
    double factor = sqrt(2.0/M);
    for(int k=0;k<M;k++){
        for(int n=0;n<M;n++){
            T[n][k] = factor*cos(M_PI*(2*n+1)*k/(2.0*M));
        }
    }
    for(int n=0;n<M;n++) T[n][0] /= sqrt(2.0);
    return T;
}

int main(){
    cout << "=== Correlation Coding: PCA vs DCT ===" << endl;

    int N=128, M=8;
    auto X = generate_random_vectors(N,M,0.8);
    export_csv("random_vectors.csv", X);

    // Covariance matrix
    auto cov = covariance_matrix(X);
    export_csv("covariance_matrix.csv", cov);

    // PCA
    vector<double> eigenvalues;
    vector<vector<double>> eigenvectors;
    jacobi_eigen(cov, eigenvalues, eigenvectors);
    export_csv("pca_eigenvectors.csv", eigenvectors);

    auto Y_pca = apply_transform(X, eigenvectors);
    export_csv("pca_transformed_vectors.csv", Y_pca);

    // DCT
    auto T_dct = dct_matrix(M);
    export_csv("dct_matrix.csv", T_dct);

    auto Y_dct = apply_transform(X, T_dct);
    export_csv("dct_transformed_vectors.csv", Y_dct);

    cout << "PCA and DCT transform coding done. CSV exported." << endl;
    return 0;
}

/*
# Python plotting code

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
*/
