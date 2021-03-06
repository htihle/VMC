#ifndef HO2_H
#define HO2_H

#include <Orbitals/orbital.h>

class HO2 : public Orbital
{
public:
    HO2()  : Orbital() {}
    double eval(arma::vec x,double a,double r);
    //double laplacian(arma::vec x,double a);
};
#endif // HO2_H
