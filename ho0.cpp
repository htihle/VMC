#include "ho0.h"
#include <iostream>
#include <armadillo>

using namespace std;

double HO0::eval(double x,double a)
{
    return exp(-x*x/(2*a*a));
}

double HO0::gradient(double x, double a)
{
    return -x/(a*a)*exp(-x*x/(2*a*a));
}

double HO0::laplacian(double x, double a)
{
    return (x*x/(a*a*a*a)- 1/(a*a))*exp(-x*x/(2*a*a));
}
