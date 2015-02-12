#include <slater.h>
#include <armadillo>

using namespace std;

using arma::mat;
using arma::vec;
using arma::randn;
using arma::randu;
using arma::zeros;


Slater::Slater(vec a, int N, int Ndim) : WaveFunction(a,N,Ndim) {
    this->a = a(0);
    this->b = a(1);
    this->whichSlater = 0;
    this->splitSlater = 2;
    //    myorbital[0] = new HO0();
    //    myorbital[1] = new HO1();
    //    myorbital[2] = new HO2();
    myorbital[0] = new H1s();
    myorbital[1] = new H2s();
    myorbital[2] = new H2p(0);
    myorbital[3] = new H2p(1);
    myorbital[4] = new H2p(2);
    NumberOfParticles = N/2;
    NumberOfDimensions = Ndim;
    R = zeros<mat>(NumberOfParticles*2,NumberOfParticles*2);
    spins = R;
    correlationsMat = R;
    D = 0.5;
    switch(NumberOfParticles*2) {
    case(2) :
        stepSize = 0.0005;
        break;
    case(4) :
        stepSize = 0.02;
        break;
    case(10) :
        stepSize = 0.05;
        break;
    default :
        break;
    }
}

double Slater::laplacianLog(mat x)
{
    double lap1 = this->slaterAnalyticalLaplacianLog(x);
    double lap = this->slaterNumericalLaplacianLog(x);
    //    cout << lap1 -lap<< endl;
    return lap1;
}
bool Slater::newStep(mat &xnew, mat x,int &whichParticle)
{
    xnew = x;
    vec num  = randn<vec>(NumberOfDimensions);
    vec num2 = randu<vec>(2);
    whichParticle = round(num2(1)*(this->NumberOfParticles*2-1));
    whichSlater = whichParticle / NumberOfParticles;
    //xnew.col(whichParticle) = x.col(whichParticle) + this->stepSize/NumberOfParticles/NumberOfDimensions*a*num; //newstep(x)
    xnew.col(whichParticle) = x.col(whichParticle)
                                + num*sqrt(stepSize)
                                + this->quantumForceOld.col(whichParticle) * stepSize * D; //newstep(x)
    Rsd = 0;
    this->updateR(xnew,whichParticle);
    for(int j = 0;j<NumberOfParticles;j++)
    {
        Rsd += myorbital[j]->eval(xnew.col(whichParticle),a, R(whichParticle, whichParticle))*
                slaterInverse[whichSlater](j,whichParticle-whichSlater*NumberOfParticles);
    }
    bool accept = false;
    this->Rc = 1;
    double coeff = Rsd*Rsd*Rc*Rc * this->computeGreensFunction(xnew, x);
    if(coeff > num2(0)){
        accept = true;
        Rold = R;
        acceptanceCounter+=1;
    }
    else
    {
        R = Rold;
    }
    return accept;
}

double Slater::evaluateSlater(mat x) {
    return det(calculateSlater(x,0))*det(calculateSlater(x,1));   //   fix for general (non-split) slater!!
}

mat Slater::calculateSlater(mat x,int upordown) {
    mat Slater = zeros<mat>(NumberOfParticles,NumberOfParticles);
    for(int i = 0;i<NumberOfParticles;i++) {
        for(int j = 0; j<NumberOfParticles;j++) {

            Slater(i,j) = myorbital[j]->eval(x.col(i+upordown*NumberOfParticles),a,R(i+upordown*NumberOfParticles,i+upordown*NumberOfParticles));
        }
    }
    return Slater;
}

double Slater::slaterNumericalLaplacianLog(mat x)  // does not work with R !!!!!!!!!!!!
{
    double h= 0.0001;
    mat mod = zeros<mat>(NumberOfDimensions, NumberOfParticles*2);
    double sum = 0;
    for(int i = 0;i<NumberOfParticles*2;i++) {
        for (int j = 0; j<NumberOfDimensions; j++) {
            mod(j,i)+=h;
            sum += this->evaluateSlater(x+mod) +this->evaluateSlater(x-mod);
            mod(j,i)=0;
        }
    }
    sum/=this->evaluateSlater(x);
    sum -=NumberOfParticles*4*NumberOfDimensions;
    sum/=h*h;
    return sum;
}
void Slater::getSlaterInverse(mat x)
{
    for (int i = 0; i < splitSlater; i++) {
        slaterInverse[i] = this->calculateSlater(x,i);
        slaterInverse[i] = inv(slaterInverse[i]);
    }
}

void Slater::makeR(mat x)
{
    for(int i=0;i<NumberOfParticles*2;i++)
    {
        this->updateR(x,i);
    }
}

// Updates the distance matrix when we have just changed one coordinate of particle "i"
// (like we do when computing the derivative)*/
void Slater::updateR(mat x, int i){
    double sum;
    for(int k=0; k<i; k++) {
        sum = norm(x.col(i)- x.col(k));
        R(k,i) = sum; //R is the matrix of distances
    }
    for(int k=i+1;k<NumberOfParticles*2;k++){
        sum = norm(x.col(i)- x.col(k));
        R(i,k) = sum; //R is the matrix of distances
    }
    sum = norm(x.col(i));
    R(i,i) = sum;
    //    cout << norm(x.col(0))- R(0,0) << endl;
}

double Slater::slaterAnalyticalLaplacianLog(mat x) {
    // we should keep the laplacian wrt all the different particles separately and update only
    // the ones we have moved (or is that correct?)
    double sum = 0;
    for (int k = 0; k < splitSlater; k++) {
        for(int i =0; i<NumberOfParticles;i++) {
            for(int j = 0;j<NumberOfParticles;j++) {
//                double p = norm(x.col(i+k*NumberOfParticles)) - R(i+k*NumberOfParticles,i+k*NumberOfParticles);
//                if (abs(p) > 0.00001) {
//                    cout << x << R << endl;
//                    cout << p << endl;
//                    cout << i+k*NumberOfParticles << endl;
//                }

                sum += myorbital[j]->laplacian(x.col(i+k*NumberOfParticles),a, R(i+k*NumberOfParticles,i+k*NumberOfParticles))
                        *slaterInverse[k](j,i); //for some reason j and i are switched
            }
        }
    }
    return sum;
}

void Slater::updateSlaterInverse(mat x, int i)
{
    //    slaterInverse[whichSlater]= calculateSlater(x,whichSlater).i();
    i = i -whichSlater*NumberOfParticles;

    mat oldSlater = slaterInverse[whichSlater];
    for (int k = 0;k<NumberOfParticles; k++) {
        for (int j = 0; j<NumberOfParticles; j++) {
            if(j != i) {
                double sum = 0;
                for( int l = 0; l< NumberOfParticles; l++) {
                    sum += oldSlater(l,j) * myorbital[l]->eval(x.col(i+whichSlater*NumberOfParticles), a, R(i+whichSlater*NumberOfParticles,i+whichSlater*NumberOfParticles));
                }
                slaterInverse[whichSlater](k,j) = oldSlater(k,j) - oldSlater(k,i) * sum / Rsd;
            }
            else {
                slaterInverse[whichSlater](k,j) = oldSlater(k,i) / Rsd;
            }
        }
    }
}

/* Compute all the correlation factors. */
void Slater::fillCorrelationsMatrix() {
    for (int i = 0; i < NumberOfParticles*2; i++) {
        for (int j = i+1; j < NumberOfParticles*2; j++) {
            correlationsMat(i,j) = spins(i,j) * R(i,j) / (1 + b * R(i,j));
        }
    }
}

/* Updates the new correlations matrix. */
void Slater::updateCorrelationsMatrix(int particle) {
    int k = particle;
    for (int i = 0; i < k; i++) {
        correlationsMat(i,k) = spins(i,k) * R(i,k) / (1 + b * R(i,k));
    }
    for (int i = (k+1); i < NumberOfParticles*2; i++) {
        correlationsMat(k,i) = spins(k,i) * R(k,i) / (1 + b * R(k,i));
    }
}

//computes the gradient of the Slater part of the wavefunction
void Slater::updateSlaterGradient(mat x, int i) {
    int k = i %(NumberOfParticles); //the remainder
    vec sum = zeros<vec>(NumberOfDimensions);
    for (int j = 0; j < (NumberOfParticles); j++) {
        sum += slaterInverse[whichSlater](j,k) * myorbital[j]->gradient(x.col(i), a, R(i,i));
    }

    gradient.col(i) = 1.0/this->Rsd * sum;
}


void Slater::computeSlaterGradient(arma::mat x) {
    for (int i = 0; i < NumberOfParticles*2; i++) {
        whichSlater = i / NumberOfParticles;
        updateSlaterGradient(x, i);
    }
}


// Maybe make this more efficient (some day...)
double Slater::computeGreensFunction(arma::mat x, arma::mat xOld) {

    this->quantumForceOld = this->gradientOld;
    this->quantumForceNew = this->gradient;

    double greensFunction = 0;
    for (int j = 0; j < NumberOfParticles*2; j++) {
        for (int i = 0; i < NumberOfDimensions; i++) {
            greensFunction += 0.5* (quantumForceOld(i,j) +
                                    quantumForceNew(i,j)) *
                    (D * this->stepSize * 0.5* (quantumForceOld(i,j) -
                                    quantumForceNew(i,j)) -
                     x(i,j) + xOld(i,j));
        }
    }
    return exp(greensFunction);
}


/*
* Fills the matrix of spin values, which is used to compute the correlation
* part of the wavefunction.
*/
void Slater::fillSpinMatrix() {
    // Set which electrons have spins up, and which have spins down.
    vec electronSpins(NumberOfParticles*2);
    electronSpins.zeros();
    for (int i = 0; i < (NumberOfParticles); i++) {
        electronSpins(i) = 1;
    }
    // Set the electron-electron spin interaction matrix.
    for (int i = 0; i < NumberOfParticles*2; i++) {
        for (int j = i+1; j < NumberOfParticles*2; j++) {
            if (electronSpins(i) != electronSpins(j)) {
                spins(i,j) = 0.5;
            } else {
                spins(i,j) = 0.25;
            }
        }
    }
}


void Slater::setUpForMetropolis(mat &x) {
    acceptanceCounter = 0;
    arma::arma_rng::set_seed(200); // not sure if this helps
    x = a*randn<mat>(this->NumberOfDimensions, this->NumberOfParticles*2);
    this->makeR(x);
    this->fillSpinMatrix();
    this->fillCorrelationsMatrix();
    Rold = R;
    this->Rsd = 1;
    this->getSlaterInverse(x);
    this->gradient = zeros<mat>(this->NumberOfDimensions, 2*this->NumberOfParticles);
    this->computeSlaterGradient(x);
    this->quantumForceNew = gradient;
    this->quantumForceOld = gradient;
    this->gradientOld = gradient;

}
