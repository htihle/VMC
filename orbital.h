#ifndef ORBITAL_H
#define ORBITAL_H

#include <armadillo>

class Orbital {
public:
    Orbital();
    virtual double eval(arma::vec x,double a) = 0;
    virtual double gradient(arma::vec x, double a);
    virtual double laplacian(arma::vec x, double a);
};

#endif // ORBITAL_H
