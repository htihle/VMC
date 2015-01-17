#ifndef HO2_H
#define HO2_H

#include <orbital.h>

class HO2 : public Orbital
{
public:
    HO2()  : Orbital() {}
    double eval(double x,double a);
    double laplacian(double x,double a);
};
#endif // HO2_H
