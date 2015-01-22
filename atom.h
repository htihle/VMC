#ifndef ATOM_H
#define ATOM_H
#include <hamiltonian.h>
#include <armadillo>

class Atom : public Hamiltonian
{
public:
    double Z;
    Atom(double Z);
    double calculatePotential(arma::mat x);
};

#endif // ATOM_H
