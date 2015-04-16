//#include "onedimensionalslater.h"
//#include <armadillo>

//using namespace arma;

//OneDimensionalSlater::OneDimensionalSlater(double a, int N) : WaveFunction(a,N,1) {
//    this->a = a;
//    myorbital[0] = new HO0();
//    myorbital[1] = new HO1();
//    myorbital[2] = new HO2();
//    NumberOfParticles = N;
//    NumberOfDimensions = 1;
//}




//double OneDimensionalSlater::laplacianLog(mat x) {
//    double analytic = this->slaterAnalyticalLaplacianLog(x);
////    double num = this->slaterNumericalLaplacianLog(x);
//    return analytic;
//}

//bool OneDimensionalSlater::newStep(mat &xnew, mat x,int &WhichParticle) {
//    xnew = x;
//    vec num = randn<vec>(1);
//    vec num2 = randu<vec>(2);
//    WhichParticle = round(num2(1)*(this->NumberOfParticles-1));
//    xnew.row(WhichParticle) = x.row(WhichParticle) + 2.0/NumberOfParticles*a*num(0);    //newstep(x)

//    Rsd = 0;
//    for(int j = 0;j<NumberOfParticles;j++) {
//        Rsd += myorbital[j]->eval(xnew.row(WhichParticle),a, R(WhichParticle,WhichParticle))*SlaterInverse(j,WhichParticle);
//    }
//    bool accept = false;
//    if(Rsd*Rsd>num2(0)){
//        accept = true;
//        this->acceptanceCounter+=1;
//    }
//    return accept;
//}

//double OneDimensionalSlater::evaluateSlater(arma::mat x) {
//    return det(calculateSlater(x));
//}

//mat OneDimensionalSlater::calculateSlater(mat x) {
//    mat Slater = zeros<mat>(NumberOfParticles,NumberOfParticles);

//    for(int i = 0;i<NumberOfParticles;i++) {
//        for(int j = 0; j<NumberOfParticles;j++) {
//            Slater(i,j) = myorbital[j]->eval(x.row(i),a, R(i,i));
//        }
//    }
//    return Slater;
//}

//double OneDimensionalSlater::slaterNumericalLaplacianLog(mat x) {
//    double h= 0.0001;
//    mat mod = zeros<mat>(NumberOfParticles);
//    double sum = 0;
//    for(int i = 0;i<NumberOfParticles;i++) {
//        mod(i)+=h;
//        sum += this->evaluateSlater(x+mod) +this->evaluateSlater(x-mod);
//        mod(i)=0;
//    }
//    sum/=this->evaluateSlater(x);
//    sum -=NumberOfParticles*2;
//    sum/=h*h;
//    return sum;
//}

//void OneDimensionalSlater::getSlaterInverse(arma::mat x) {
//    SlaterInverse = this->calculateSlater(x);
//    SlaterInverse = inv(SlaterInverse);
//}

//double OneDimensionalSlater::slaterAnalyticalLaplacianLog(arma::mat x) {
//    // we should keep the laplacian wrt all the different particles separately and update only
//    // the ones we have moved (or is that correct?)
//    double sum = 0;
//    for(int i =0; i<NumberOfParticles;i++) {
//        for(int j = 0;j<NumberOfParticles;j++) {
//            sum += myorbital[j]->laplacian(x.row(i),a)*SlaterInverse(j,i); //for some reason j and i are switched
//        }
//    }
//    return sum;
//}

//void OneDimensionalSlater::updateSlaterInverse(arma::mat x, int i) {
//    mat oldSlater = SlaterInverse;
//    for (int k = 0;k<NumberOfParticles; k++) {
//        for (int j = 0; j<NumberOfParticles; j++) {
//            if(j != i) {
//                double sum = 0;
//                for( int l = 0; l< NumberOfParticles; l++) {
//                    sum += oldSlater(l,j) * myorbital[l]->eval(x.row(i), a);
//                }
//                SlaterInverse(k,j) = oldSlater(k,j) - oldSlater(k,i) * sum / Rsd;
//            }
//            else {
//                SlaterInverse(k,j) = oldSlater(k,i)/ Rsd;
//            }
//        }
//    }
//}

//void OneDimensionalSlater::setUpForMetropolis(arma::mat &x) {
//    acceptanceCounter = 0;
//    arma_rng::set_seed_random(); // not sure if this helps
//    x = a*randn<mat>(this->NumberOfParticles,this->NumberOfDimensions);
//    this->getSlaterInverse(x);
//}

