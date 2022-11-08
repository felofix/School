//
//  problem3.cpp
//  project2
//
//  Created by Felix Aarekol Forseth on 20/09/2022.
//

#include <stdio.h>
#include "utilities.hpp"

int main(int argc, const char * argv[]){
    // Problem 3.
    // Constants.
    int N = 4;
    int k = 0;
    int l = 0;
    
    // Testing maxoffvalue. 
    arma::mat A(N, N, fill::zeros);
    arma::mat R(N, N, fill::zeros);
    A(0, 0) = A(1, 1) = A(2, 2) = A(3, 3) = 1;
    A(3, 0) = A(0, 3) = 0.5;
    A(1, 2) = A(2, 1) = -0.7;
    double maxval = maxoffvalue(A, k, l, N);
    std::cout << maxval << " " << k << " " << l << std::endl;
    return 0;
}
