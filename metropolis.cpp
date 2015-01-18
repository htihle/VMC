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
    }
    expect->ev /=(n+1);
    cout << "Acceptance rate: " << double(wave->acceptanceCounter)/n << endl;
}

