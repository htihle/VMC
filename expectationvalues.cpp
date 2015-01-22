#include "expectationvalues.h"
#include <armadillo>

using namespace arma;

ExpectationValues::ExpectationValues(int numberofev, WaveFunction *wave, Hamiltonian *ham)
{
    numberofEVs = numberofev;
    ev = zeros<vec>(numberofEVs);
    this->wave =wave;
    this->ham  = ham;
    LocalEnergy = 0;
}
ExpectationValues::ExpectationValues()
{
}

void ExpectationValues::Sample(mat &xnew,mat &x, int WhichParticle)
{

    x = xnew;
    wave->updateSlaterInverse(x,WhichParticle);
    vec f=zeros<vec>(numberofEVs);
    LocalEnergy = calculateEnergy(x);
    f(0) = LocalEnergy;
    f(1) = LocalEnergy*LocalEnergy;
    f(2) = arma::norm(x.col(0));
    f(3) = arma::dot(x.col(0), x.col(0));
    this->ev +=f;
}

void ExpectationValues::ReSample(mat &xnew, mat &x)
{
    xnew = x;
    vec f=zeros<vec>(numberofEVs);
    f(0) = LocalEnergy;
    f(1) = LocalEnergy*LocalEnergy;
    f(2) = x(0);
    f(3) = x(0)*x(0);
    this->ev +=f;
}

double ExpectationValues::calculateEnergy(mat x) {
//    double sum = 0;
//    double Z = wave->NumberOfParticles*2;
//    for(int i = 0;i<wave->NumberOfParticles*2;i++)
//    {

//        sum -= Z/ (norm(x.col(i)));
//        for(int j= i+1;j<wave->NumberOfParticles*2;j++)
//        {
//            sum += 1/ (norm(x.col(i)-x.col(j)));
//        }

//        for(int j= 0;j<wave->NumberOfDimensions;j++)
//        {
//            sum += 1.0/2*x(j,i)*x(j,i);
//        }
//    }
    return -wave->laplacianLog(x)/2.0 + ham->calculatePotential(x);
}
