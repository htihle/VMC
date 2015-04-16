#pragma once
#include <armadillo>
#include <Hamiltonians/hamiltonian.h>

class HarmonicOscillator : public Hamiltonian
{
public:
    int nParticles;
    double omega;

    HarmonicOscillator(int nParticles, double omega);

    double calculatePotential(arma::mat x);
};

