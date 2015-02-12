#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H
#include <orbital.h>
#include <ho0.h>
#include <ho1.h>
#include <ho2.h>
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
//    virtual void getSlaterInverse(arma::vec x);
    virtual void updateSlaterInverse(arma::mat x, int i)=0;
    virtual void setUpForMetropolis(arma::mat &x) = 0;
};

#endif // WAVEFUNCTION_H
