#ifndef HO0_H
#define HO0_H
#include <orbital.h>

class HO0 : public Orbital
{
public:
    HO0()  : Orbital() {}
    double eval(double x, double a);
    double gradient(double x, double a);
    double laplacian(double x,double a);
};

#endif // HO0_H
