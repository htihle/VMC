#include "ho1.h"
#include <iostream>
#include <armadillo>

using namespace std;

double HO1::eval(arma::vec x,double a,double r)
{
    return r*exp(-r*r/(2*a*a));
}

double HO1::laplacian(arma::vec x, double a,double r)
{
    double rr = r*r;
//    double h = 0.0001;
//    cout << (this->eval(x+h,a) + this->eval(x-h,a) -2*this->eval(x,a))/(h*h)- (x*x/(a*a*a*a)- 3.0/(a*a))*x*exp(-x*x/(2*a*a)) << endl;
    return (rr/(a*a*a*a)- 3.0/(a*a))*r*exp(-rr/(2*a*a));
}
