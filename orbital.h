#ifndef ORBITAL_H
#define ORBITAL_H

class Orbital
{
public:
    Orbital();
    virtual double eval(double x,double a);
    virtual double gradient(double x, double a);
    virtual double laplacian(double x, double a);
};

#endif // ORBITAL_H
