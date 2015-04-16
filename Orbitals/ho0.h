#ifndef HO0_H
#define HO0_H
#include <Orbitals/orbital.h>

class HO0 : public Orbital
{
public:
    HO0()  : Orbital() {}
    double eval(arma::vec x, double a,double r);
//    double gradient(arma::vec x, double a);
//    double laplacian(arma::vec x,double a);
};

#endif // HO0_H
