#ifndef HO1_H
#define HO1_H
#include <orbital.h>

class HO1 : public Orbital
{
public:
    HO1()  : Orbital() {}
    double eval(double x,double a);
    double laplacian(double x,double a);
};

#endif // HO1_H
