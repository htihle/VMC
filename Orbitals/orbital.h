#ifndef ORBITAL_H
#define ORBITAL_H

#include <armadillo>

class Orbital {
public:
    Orbital();
    virtual double eval(arma::vec x,double a, double r) = 0;
    virtual arma::vec gradient(arma::vec x, double a, double r);
    virtual double laplacian(arma::vec x, double a, double r);
};

#endif // ORBITAL_H
