
//  problem2.cpp
//  project2
//
//  Created by Felix Aarekol Forseth on 20/09/2022.
//

#include <stdio.h>
#include "utilities.hpp"

int main(int argc, const char * argv[]){
    // Problem 2.
    // Constants.
    double n = 10;
    double h = 1/n;
    double a = -1/pow(h, 2);
    double d = 2/pow(h, 2);
    int N = 6;
    double error = 0.002;
    
    // Armadillo eigenvalues and eigenvectors.
    arma::vec eigval;
    arma::mat eigvec;
    arma::mat A = tridiagonal(a, d, a, N);
    arma::eig_sym(eigval, eigvec, A);
    
    // Analytical eigenvalues and eigenvectors.
    arma::vec aeigvalues = analeigenval(a, d, N);
    arma::mat aeigvectors(N, N, fill::zeros);
    for (int i = 0; i < N; i++){        // Creaating eigenvectors. 
        aeigvectors.col(i) = analeigenvec(N, i + 1);
    }
    
    // Checking if they have similar eigenvalues and eigenvectors.
    bool sameval = checkequaleig(eigval, aeigvalues, error);
    bool samevec = checkequalmatrix(abs(arma::normalise(eigvec)), abs(arma::normalise(aeigvectors)), error);
    if (samevec == true and sameval == true){
        std::cout << "Armadillo and analytical are equal and good friends :)" << std::endl;
    }
    return 0;
}

