#include "metropolis.h"
#include <armadillo>
#include <expectationvalues.h>
#include <wavefunction.h>

using namespace arma;
using namespace std;

Metropolis::Metropolis(ExpectationValues *expect, WaveFunction *wave)
{
    this->expect = expect;
    this->wave = wave;
}

void Metropolis::Run(int n) {

    mat xnew,x;

    wave->setUpForMetropolis(x);

    for(int i =0;i<n; i++) {
        if(wave->newStep(xnew,x, WhichParticle)){
            expect->Sample(xnew,x,WhichParticle);
        } else {
            expect->ReSample(xnew,x);
        }
        if(i == (n/10)) expect->ev = zeros<vec>(4);
    }
    expect->ev /=((n/10)*9+1);
    cout << "Acceptance rate: " << double(wave->acceptanceCounter)/n << endl;
}

