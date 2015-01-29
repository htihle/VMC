#pragma once
#include <orbital.h>
#include <armadillo>

class H2s : public Orbital
{
public:
    H2s()  : Orbital() {}
    double eval(arma::vec x, double a);
};

