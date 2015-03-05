#pragma once
#include <wavefunction.h>
#include <armadillo>


class HeliumWaveFunction : public WaveFunction {

public:
//    int       NumberOfParticles; // Number of particles in each slater, i.e. numberOfTotalParticles/2.
//    int       NumberOfDimensions;
    double    a,b,Rsd,Rc,stepSize,D,r1,r2,r12,wfold,wfnew; // a and b are Alpha and Beta

    bool      interacting;
    HeliumWaveFunction(arma::vec a, int N,int Ndim, bool interacting);
    bool      newStep                     (arma::mat &xnew, arma::mat x, int &whichParticle);
    double    laplacianLog                (arma::mat x);
    void      setUpForMetropolis          (arma::mat &x);
    void      computeR                    (arma::mat x);
    double    wf();
};
