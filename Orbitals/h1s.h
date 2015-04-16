#pragma once
#include <Orbitals/orbital.h>
#include <armadillo>

class H1s : public Orbital
{
public:
    H1s()  : Orbital() {}
    double eval(arma::vec x, double a, double r);
};
