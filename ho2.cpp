#include "ho2.h"
#include <iostream>
#include <armadillo>

using namespace std;

double HO2::eval(arma::vec x,double a)
{
    return (2*arma::dot(x,x) -1)*exp(-arma::dot(x,x)/(2*a*a));
}

//double HO2::laplacian(arma::vec x, double a)
//{
//    double h = 0.0001;
//    return (this->eval(x+h,a) + this->eval(x-h,a) -2*this->eval(x,a))/(h*h);
//}
