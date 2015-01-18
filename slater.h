#pragma once
#include <wavefunction.h>
#include <orbital.h>
#include <ho0.h>
#include <ho1.h>
#include <ho2.h>
#include <armadillo>


class Slater : public WaveFunction {

public:
    int       NumberOfParticles; // Number of particles in each slater, i.e. numberOfTotalParticles/2.
    int       NumberOfDimensions;
    int       acceptanceCounter;
    int       whichSlater;
    int       splitSlater;
    double    a,Rsd;
    Orbital*  myorbital[5];
    arma::mat slaterInverse[2];

    Slater(double a, int N);

    bool      newStep                     (arma::mat &xnew, arma::mat x, int &whichParticle);
    double    laplacianLog                (arma::mat x);
    double    evaluateSlater              (arma::mat x);
    double    slaterNumericalLaplacianLog (arma::mat x);
    double    slaterAnalyticalLaplacianLog(arma::mat x);
    void      updateSlaterInverse         (arma::mat x, int i);
    void      setUpForMetropolis          (arma::mat &x);
    void      getSlaterInverse            (arma::mat x);
    arma::mat calculateSlater             (arma::mat x);
};
