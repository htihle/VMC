#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H
#include <armadillo>

class WaveFunction
{
public:
    double a,Rsd;
    int NumberOfParticles,NumberOfDimensions;
    double acceptanceCounter;
    WaveFunction(arma::vec a, int N, int Ndim);
    WaveFunction() {}
    virtual double laplacianLog(arma::mat x) = 0;
    virtual bool newStep(arma::mat &xnew, arma::mat x, int &WhichParticle) = 0;
    virtual void setUpForMetropolis(arma::mat &x) = 0;
};

#endif // WAVEFUNCTION_H
