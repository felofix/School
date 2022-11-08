//
//  problem6.cpp
//  project2
//
//  Created by Felix Aarekol Forseth on 21/09/2022.
//

#include <stdio.h>
#include "utilities.hpp"

int main(int argc, const char * argv[]){
    // Problem 6 x = 10.
    
    // Constants.
    double n = 100;
    double h = 1/n;
    double a = -1/pow(h, 2);
    double d = 2/pow(h, 2);
    double error = 0.002;
    
    // Creaating eigenvalues, eigenvectors and armadillo vectors.
    int N = n - 1;
    double eps = 0.001;
    arma::vec xhat = arma::linspace(0, 1, N);
    arma::vec eigenvalues(N, fill::zeros);
    arma::mat eigenvectors(N, N, fill::zeros);
    arma::mat A = tridiagonal(a, d, a, N);
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, A);
    int maxiter = 10e3;
    int iterations = 0;
    bool converged = true;
    jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged, N);
    arma::uvec sort = arma::sort_index(eigenvalues, "ascend");

    // Writing to file.
    for (int i = 0; i < 3; i++){
        std::string direc = "vectors/eigenvector";
        direc += std::to_string(i+10) + ".txt";
        eigenvectors.col(sort(i)) = minusconst(eigenvectors.col(sort(i)), eigvec.col(i));
        writetofilefloat(xhat, eigenvectors.col(sort(i)), direc);
    }
    
    for (int i = 0; i < 3; i++){
        std::string direc = "vectors/eigenvectoranal";
        direc += std::to_string(i+10) + ".txt";
        writetofilefloat(xhat, eigvec.col(i), direc);
    }
    return 0;
}
