#pragma once
#include <orbital.h>
#include <armadillo>

class H2p : public Orbital
{
public:
    H2p(int state);
    double eval(arma::vec x, double a);
    int state;
};

