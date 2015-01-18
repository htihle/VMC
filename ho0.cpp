#include "ho0.h"
#include <iostream>
#include <armadillo>

using namespace std;

double HO0::eval(arma::vec x,double a)
{
    return exp(-arma::dot(x,x)/(2*a*a));
}

//double HO0::gradient(arma::vec x, double a)
//{
//    return -arma::norm(x)/(a*a)*exp(-arma::dot(x,x)/(2*a*a));
//}

//double HO0::laplacian(arma::vec x, double a)
//{
//    return (arma::dot(x,x)/(a*a*a*a)- 1/(a*a))*exp(-arma::dot(x,x)/(2*a*a));
//}
