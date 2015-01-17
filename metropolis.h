#ifndef METROPOLIS_H
#define METROPOLIS_H
#include <armadillo>
#include <expectationvalues.h>
#include <wavefunction.h>

class Metropolis
{
public:
    ExpectationValues *expect;
    WaveFunction *wave;
    int WhichParticle;
    Metropolis(ExpectationValues *expect, WaveFunction *wave);
    void Run(int n);
};

#endif // METROPOLIS_H
