#include "heliumwavefunction.h"
#include <armadillo>

using namespace std;

using arma::mat;
using arma::vec;
using arma::randn;
using arma::randu;
using arma::zeros;

HeliumWaveFunction::HeliumWaveFunction(vec a, int N, int Ndim, bool interacting) : WaveFunction(a,N,Ndim)
{
    this->interacting = interacting;
    this->a = a(0);
    this->b = a(1);
    NumberOfParticles = N;
    NumberOfDimensions = Ndim;
    D = 0.5;
    stepSize = 0.0005;
}

void HeliumWaveFunction::setUpForMetropolis(mat &x) {
    acceptanceCounter = 0;
//    arma::arma_rng::set_seed(200); // not sure if this helps
    x = a*randn<mat>(this->NumberOfDimensions, this->NumberOfParticles);
    this->computeR(x);
    wfold = this->wf();
}
double HeliumWaveFunction::laplacianLog(mat x)
{
    if(this->interacting){
        return 0;
    } else {
        return  2*a*a - a*(2/r1 +2/r2);
    }
}

void HeliumWaveFunction::computeR(mat x)
{
    double sum = 0;
    r1 = norm(x.col(0));
    r2 = norm(x.col(1));
    r12 = norm(x.col(0)-x.col(1));
}
bool HeliumWaveFunction::newStep(mat &xnew, mat x,int &whichParticle) {
    double Rcoeff= 0;
    xnew = x;
    vec num  = randn<vec>(NumberOfDimensions);
    vec num2 = randu<vec>(2);
    whichParticle = round(num2(1)*(this->NumberOfParticles-1));
    //xnew.col(whichParticle) = x.col(whichParticle) + this->stepSize/NumberOfParticles/NumberOfDimensions*a*num; //newstep(x)
    xnew.col(whichParticle) = x.col(whichParticle)
            + num*sqrt(stepSize);
            //+ this->quantumForceOld.col(whichParticle) * stepSize * D; //newstep(x)
    this->computeR(xnew);
    wfnew = this->wf();
    Rsd = wfnew/wfold;
    bool accept = false;
    if(interacting){
//        this->computeRc(whichParticle);
        Rcoeff = Rsd*Rc;
    } else {
        Rcoeff = Rsd;
    }

    double coeff = Rcoeff*Rcoeff; // * this->computeGreensFunction(xnew, x);

    if(coeff > num2(0)){
        accept = true;
        wfold = wfnew;
//        gradientOld        = gradient;
//        quantumForceOld    = quantumForceNew;
        acceptanceCounter +=1;
    } else {
        wfnew = wfold;
//        gradient           = gradientOld;
//        quantumForceNew    = quantumForceOld;
    }
    return accept;
}
double HeliumWaveFunction::wf(){
    return exp(-a*(r1 +r2));
}
