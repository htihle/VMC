#include "wavefunction.h"
#include <armadillo>

using namespace arma;

WaveFunction::WaveFunction(vec a, int N, int Ndim)
{
    this->a = a(0);
    this->acceptanceCounter = 0;
    NumberOfParticles = N;
    NumberOfDimensions = Ndim;
}


//double WaveFunction::laplacianLog(arma::mat x)
//{
//    cout << " error, in parent class 'Wavefunction' 1" << endl;
//    return 0;
//}

//bool WaveFunction::newStep(mat &xnew, mat x,int &WhichParticle)
//{
//    cout << " error, in parent class 'Wavefunction' 2" << endl;
//    return true;
//}

//void WaveFunction::getSlaterInverse(arma::vec x)
//{
//    cout << " error, in parent class 'Wavefunction' " << endl;
//}

//void WaveFunction::updateSlaterInverse(arma::mat x, int i)
//{
//    cout << " error, in parent class 'Wavefunction' 3" << endl;
//}

//void WaveFunction::setUpForMetropolis(arma::mat &x)
//{
//    cout << " error, in parent class 'Wavefunction' 4" << endl;
//}
