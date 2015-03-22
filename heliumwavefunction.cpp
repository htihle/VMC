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
    stepSize = 0.001;
}

void HeliumWaveFunction::setUpForMetropolis(mat &x) {
    acceptanceCounter = 0;
    arma::arma_rng::set_seed(200); // not sure if this helps
    x = a*randn<mat>(this->NumberOfDimensions, this->NumberOfParticles);
    this->computeR(x);
    wfold = this->wf();
}

double HeliumWaveFunction::laplacianLog(mat x) {

    double EL1 = 2*a*a - a*(2/r1 +2/r2);
    double EL2 = -2 * (1/(2*(1+b*r12)*(1+b*r12)) * (a*(r1+r2)/r12 * (1 - dot(x.col(0),x.col(1))/(r1*r2)) - 1/(2*(1+b*r12)*(1+b*r12)) -2/r12 + 2*b/(1+b*r12)));


    // numerical laplacian
//    double numE = 0;
//    double sum = 0;

//    double h = 0.0001;
//    double h2 = h*h;
//    for (int j = 0; j< NumberOfParticles;j++){
//        numE = 0;
//        for (int i = 0; i < NumberOfDimensions; i++) {
//            x(i,j) += h;
//            numE += wf(x);
//            x(i,j) -= 2*h;
//            numE += wf(x);
//            x(i,j) += h;
//        }
//        numE /=this->wf();
//        numE -= 6.0;
//        numE /= h2;
//        sum += numE;
//    }
//    cout << EL1-sum << endl;
    if(this->interacting) {
        return EL1 + EL2;
    } else {
        return EL1;
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

    xnew = x;
    vec num  = randn<vec>(NumberOfDimensions);
    vec num2 = randu<vec>(2);
    whichParticle = round(num2(1)*(this->NumberOfParticles-1));
    //xnew.col(whichParticle) = x.col(whichParticle) + this->stepSize/NumberOfParticles/NumberOfDimensions*a*num; //newstep(x)
    xnew.col(whichParticle) = x.col(whichParticle)
            + num*sqrt(stepSize)
            + this->QF(whichParticle,x) * stepSize * D; //newstep(x)
//    cout << this->QF(whichParticle,x) << " " << whichParticle << endl;
    this->computeR(xnew);
    wfnew = this->wf();
    R = wfnew/wfold;
    bool accept = false;


    if(true) { //R*R > num2(0)){
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

double HeliumWaveFunction::wf() {
    double waveFun = exp(-a*(r1 +r2));
    if (this->interacting) {
        return waveFun * exp(r12 / (2*(1+b*r12)));
    } else {
        return waveFun;
    }
}

double HeliumWaveFunction::wf(mat X) {
    double r1_ = norm(X.col(0));
    double r2_ = norm(X.col(1));
    double r12_ = norm(X.col(0)-X.col(1));

    double waveFun = exp(-a*(r1_ +r2_));
    if (this->interacting) {
        return waveFun * exp(r12_ / (2*(1+b*r12_)));
    } else {
        return waveFun;
    }
}

vec HeliumWaveFunction::QF(int whichParticle, mat x){

    if (true) { //this->interacting) {
        vec qf = zeros<vec>(NumberOfDimensions);

        double h = 0.00001;
        for (int i = 0; i < NumberOfDimensions; i++) {
            x(i,whichParticle) += h;
            qf(i) = wf(x);
            x(i,whichParticle) -= h;
            qf(i) -= wf(x);
        }
        qf = qf / (this->wf()*h);
        return 2*qf;
    } else  {
        vec qf2 = zeros<vec>(NumberOfDimensions);
        if(whichParticle == 0) {
            qf2 = -2*x.col(0)/r1;
            return qf2;
        } else {
            qf2 = -2*x.col(1)/r2;
            return qf2;
        }
    }



//    for (int i = 0; i < NumberOfDimensions; i++) {
//        cout << fabs(qf(i)-qf2(i)) << ", ";
//    }
//    cout << endl;
//    return qf2;
}

