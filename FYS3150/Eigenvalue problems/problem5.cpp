//
//  problem5.cpp
//  project2
//
//  Created by Felix Aarekol Forseth on 20/09/2022.
//

#include <stdio.h>
#include "utilities.hpp"

int main(int argc, const char * argv[]){
    // Problem 5.
    // Constants.
    std::string direction = "iterations.txt";
    double n = 7;
    double h = 1/n;
    double a = -1/pow(h, 2);
    double d = 2/pow(h, 2);
    double eps = 0.001;
    int iterations = 0;
    std::vector<int> iterationslist;
    std::vector<int> N;
    int maxiter = 10e4;
    bool converged = true;
    
    // Iterating over different values of N.
    for (int Ni = 1; Ni <= 120; Ni++){
        arma::vec eigenvalues(Ni, fill::zeros);
        arma::mat eigenvectors(Ni, Ni, fill::zeros);
        arma::mat A = tridiagonal(a, d, a, Ni);
        jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged, Ni);
        N.push_back(Ni);
        iterationslist.push_back(iterations);
        iterations = 0;
    }
    
    // Writing results to file. 
    writetofile(N, iterationslist, direction);
    return 0;
}

