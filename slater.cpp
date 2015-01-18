#include <slater.h>
#include <armadillo>

using namespace std;

using arma::mat;
using arma::vec;
using arma::randn;
using arma::randu;
using arma::zeros;


Slater::Slater(double a, int N,int Ndim) : WaveFunction(a,N,Ndim) {
    this->a = a;
    this->whichSlater = 0;
    this->splitSlater = 2;
    myorbital[0] = new HO0();
//    myorbital[1] = new HO1();
//    myorbital[2] = new HO2();
//    myorbital[0] = new H1s();

    NumberOfParticles = N/2;
    NumberOfDimensions = Ndim;
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
    xnew.col(whichParticle) = x.col(whichParticle) + 0.5/NumberOfParticles*a*num; //newstep(x)
    Rsd = 0;

    for(int j = 0;j<NumberOfParticles;j++)
    {
        Rsd += myorbital[j]->eval(xnew.col(whichParticle),a)*
                             slaterInverse[whichSlater](j,whichParticle-whichSlater*NumberOfParticles);
    }
    bool accept = false;
    if(Rsd*Rsd>num2(0)){
        accept = true;
        acceptanceCounter+=1;
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

            Slater(i,j) = myorbital[j]->eval(x.col(i+upordown*NumberOfParticles),a);
        }
    }
    return Slater;
}

double Slater::slaterNumericalLaplacianLog(mat x)
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

double Slater::slaterAnalyticalLaplacianLog(mat x) {
    // we should keep the laplacian wrt all the different particles separately and update only
    // the ones we have moved (or is that correct?)
    double sum = 0;
    for (int k = 0; k < splitSlater; k++) {
        for(int i =0; i<NumberOfParticles;i++) {
            for(int j = 0;j<NumberOfParticles;j++) {
                sum += myorbital[j]->laplacian(x.col(i+k*NumberOfParticles),a)
                                     *slaterInverse[k](j,i); //for some reason j and i are switched
            }
        }
    }
    return sum;
}

void Slater::updateSlaterInverse(mat x, int i)
{
    slaterInverse[whichSlater]= calculateSlater(x,whichSlater).i();
    //mat oldSlater = slaterInverse[whichSlater];

    /*mat old = zeros<mat>(4,4);
    mat NEW = zeros<mat>(4,4);
    mat oldSlater = zeros<mat>(4,4);
    mat newSlater = zeros<mat>(4,4);

    mat addd = zeros<mat>(4,4);

    addd.(2) += 1;


    for (int k = 0;k<4; k++) {
        for (int j = 0; j<4; j++) {
            old(k,j) = 1+k*24.12+j*5+j*j*7.1 - 5.245214;
        }
    }

    NEW = old+addd;


    oldSlater = old.i();
    newSlater = NEW.i();


    cout << "old:" << old << endl;
    cout << "new:" << NEW << endl;
    cout << "old.i():" << oldSlater << endl;
    cout << "mat.i():" << newSlater << endl;

*/


/*

    for (int k = 0;k<NumberOfParticles; k++) {
        for (int j = 0; j<NumberOfParticles; j++) {
            if(j != i) {
                double sum = 0;
                for( int l = 0; l< NumberOfParticles; l++) {
                    sum += oldSlater(l,j) * myorbital[l]->eval(x.col(i), a);
                }
                slaterInverse[whichSlater](k,j) = oldSlater(k,j) - oldSlater(k,i-whichSlater*NumberOfParticles) * sum / Rsd;
            }
            else {
                slaterInverse[whichSlater](k,j) = oldSlater(k,i-whichSlater*NumberOfParticles) / Rsd;
            }
        }
    }*/
}

void Slater::setUpForMetropolis(mat &x)
{
    acceptanceCounter = 0;
     // not sure if this helps
    x = a*randn<mat>(this->NumberOfDimensions, this->NumberOfParticles*2);
    this->getSlaterInverse(x);
}
