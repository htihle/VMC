#pragma once
#include <wavefunction.h>
#include <orbital.h>
#include <ho0.h>
#include <ho1.h>
#include <ho2.h>
#include <h1s.h>
#include <h2s.h>
#include <h2p.h>
#include <armadillo>


class Slater : public WaveFunction {

public:
//    int       NumberOfParticles; // Number of particles in each slater, i.e. numberOfTotalParticles/2.
//    int       NumberOfDimensions;
    int       whichSlater;
    int       splitSlater;
    int       sizeOfSlater;
    double    a,b,Rsd,Rc,stepSize,D; // a and b are Alpha and Beta
    Orbital*  myorbital[5];
    arma::mat slaterInverse[2];
    arma::mat R,Rold,spins,correlationsMatNew,correlationsMatOld;
    arma::mat gradient, gradientOld, Jgradient, JgradientOld;
    arma::mat quantumForceNew, quantumForceOld;
    arma::mat jastrowLaplacian, jastrowGradient;
    bool      interacting;


    Slater(arma::vec a, int N,int Ndim, bool interacting);

    bool      newStep                     (arma::mat &xnew, arma::mat x, int &whichParticle);
    double    laplacianLog                (arma::mat x);
    double    evaluateSlater              (arma::mat x);
    double    slaterNumericalLaplacianLog (arma::mat x);
    double    slaterAnalyticalLaplacianLog(arma::mat x);
    double    computeGreensFunction       (arma::mat x, arma::mat xOld);
    double    computeJastrowEnergy        (arma::mat& jastrowGradient);
    void      updateSlaterInverse         (arma::mat x, int i);
    void      setUpForMetropolis          (arma::mat &x);
    void      getSlaterInverse            (arma::mat x);
    void      makeR                       (arma::mat x);
    void      updateR                     (arma::mat x, int i);
    void      fillSpinMatrix              ();
    void      fillCorrelationsMatrix      ();
    void      updateCorrelationsMatrix    (int i);
    void      updateSlaterGradient        (arma::mat x, int i);
    void      computeSlaterGradient       (arma::mat x);
    void      makeJgradient               (arma::mat x);
    void      computeJastrowGradient      (int particle);
    void      computeJastrowLaplacian     (int particle);
    void      computeRc                   (int i);
    void      setUpJastrow                ();
    arma::mat calculateSlater             (arma::mat x, int upordown);
};
