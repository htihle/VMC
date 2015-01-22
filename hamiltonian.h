#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H
#include <armadillo>
class Hamiltonian
{
public:
    Hamiltonian();
    virtual double calculatePotential(arma::mat x)=0;
};

#endif // HAMILTONIAN_H
