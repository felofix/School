//
//  utilities.hpp
//  project2
//
//  Created by Felix Aarekol Forseth on 14/09/2022.
//

#ifndef utilities_hpp
#define utilities_hpp

#include <stdio.h>
#include <vector>
#include <armadillo>
#include <math.h>
#include <typeinfo>
using namespace arma;

arma::mat tridiagonal(double a, double b, double c, int N); // Creating a tridaigonal matrix.
arma::vec analeigenval(double a, double d, int N); // Creating analytical eigenvalues.
arma::vec analeigenvec(int N, int i); // Creating analytical eigenvectors.
bool checkequaleig(arma::vec A, arma::vec B, double error); // Checking if vectors are equal.
bool checkequalmatrix(arma::mat A, arma::mat B, double error); // Checking if vectors are equal.
double maxoffvalue(const arma::mat& A, int& k, int& l, int N); // Finding largest off-diagonal value.
double jacobi_rotate(arma::mat& A, arma::mat& R, int& k, int& l, int N); // One iteration of Jacobi-rotation.
void jacobi_eigensolver(arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, // Solving Jacobi. 
                                                   const int maxiter, int& iterations, bool& converged, int N);
void writetofile(std::vector<float> xvalues, std::vector<float> ux, std::string direc);
void writetofilefloat(arma::vec xvalues, arma::vec ux, std::string direc);
arma::vec minusconst(arma::vec a, arma::vec b); // Checking if there is a factor -1 mistake. 
#endif /* utilities_hpp */
