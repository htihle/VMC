#include "ho2.h"
#include <iostream>
#include <armadillo>

using namespace std;

double HO2::eval(double x,double a)
{
    return (2*x*x -1)*exp(-x*x/(2*a*a));
}

double HO2::laplacian(double x, double a)
{
    double h = 0.0001;
    return (this->eval(x+h,a) + this->eval(x-h,a) -2*this->eval(x,a))/(h*h);
}
