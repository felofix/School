//
//  problem4.cpp
//  project2
//
//  Created by Felix Aarekol Forseth on 20/09/2022.
//

#include <stdio.h>
#include "utilities.hpp"

int main(int argc, const char * argv[]){
    // Problem 4.
    double n = 7;
    double h = 1/n;
    double a = -1/pow(h, 2);
    double d = 2/pow(h, 2);
    double error = 0.02;
    
    // Creaating eigenvalues, eigenvectors and armadillo vectors. 
    int N = 6;
    double eps = 1e-10;
    arma::vec eigenvalues(N, fill::zeros);
    arma::mat eigenvectors(N, N, fill::zeros);
    arma::mat A = tridiagonal(a, d, a, N);
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, A);
    int maxiter = 10e3;
    double maxerror = 1e-10;
    int iterations = 0;
    bool converged = true;
    jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged, N);
    
    // Sorting:
    arma::uvec sort = arma::sort_index(eigenvalues, "ascend");
    
    // Checking if similar.
    bool sameval = true;
    bool samevec = true;
    
    for (int i = 0; i < N; i ++){ // wtf?
        if (sameval == false && samevec == false){
            std::cout << "They were not similar :(" << std::endl;
        }
        if (abs(eigenvalues(sort(i)) - eigval(i)) < maxerror){
            sameval = true;
        }
        else{
            sameval = false;
        }
        samevec = checkequaleig(arma::normalise(abs(eigenvectors.col(sort(i)))), abs(arma::normalise(eigvec.col(i))), 0.0002);
    }
    
    // Letting know if they are similar.
    if (sameval == samevec == true){
        std::cout << "Jacobi and analytical are equal and good friends :)" << std::endl;
    }
}
