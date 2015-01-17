#ifndef ONEDIMENSIONALSLATER_H
#define ONEDIMENSIONALSLATER_H
#include <wavefunction.h>
#include <orbital.h>
#include <ho0.h>
#include <ho1.h>
#include <ho2.h>
#include <armadillo>

class OneDimensionalSlater : public WaveFunction
{
public:
    OneDimensionalSlater(double a, int N);  // : WaveFunction(a,N) {}
    double a,Rsd;
    int NumberOfParticles,NumberOfDimensions,acceptanceCounter;
    Orbital* myorbital[5];
    arma::mat SlaterInverse;
    double laplacianLog(arma::vec x);
    bool newStep(arma::mat &xnew, arma::mat x, int &WhichParticle);
    double evaluateSlater(arma::vec x);
    arma::mat calculateSlater(arma::vec x);
    double slaterNumericalLaplacianLog(arma::vec x);
    double slaterAnalyticalLaplacianLog(arma::vec x);
    void getSlaterInverse(arma::vec x);
    void updateSlaterInverse(arma::vec x, int i);
    void setUpForMetropolis(arma::mat &x);
};

#endif // ONEDIMENSIONALSLATER_H
