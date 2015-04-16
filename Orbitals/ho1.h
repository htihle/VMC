#ifndef HO1_H
#define HO1_H
#include <Orbitals/orbital.h>

class HO1 : public Orbital
{
public:
    HO1()  : Orbital() {}
    double eval(arma::vec x,double a,double r);
    double laplacian(arma::vec x,double a,double r); //
};

#endif // HO1_H
