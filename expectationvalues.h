#ifndef EXPECTATIONVALUES_H
#define EXPECTATIONVALUES_H
#include <armadillo>
#include <wavefunction.h>

class ExpectationValues
{
public:
    int numberofEVs;
    double LocalEnergy;
    arma::vec ev;
    WaveFunction *wave;
    ExpectationValues(int numberofev, WaveFunction *wave);
    ExpectationValues();
    void Sample(arma::mat &xnew, arma::mat &x, int WhichParticle);
    void ReSample(arma::mat &xnew,arma::mat &x);
    double calculateEnergy(arma::mat x);
};

#endif // EXPECTATIONVALUES_H
