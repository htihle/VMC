#include "ho1.h"
#include <iostream>
#include <armadillo>

using namespace std;

double HO1::eval(double x,double a)
{
    return x*exp(-x*x/(2*a*a));
}

double HO1::laplacian(double x, double a)
{

//    double h = 0.0001;
//    cout << (this->eval(x+h,a) + this->eval(x-h,a) -2*this->eval(x,a))/(h*h)- (x*x/(a*a*a*a)- 3.0/(a*a))*x*exp(-x*x/(2*a*a)) << endl;
    return (x*x/(a*a*a*a)- 3.0/(a*a))*x*exp(-x*x/(2*a*a));
}
