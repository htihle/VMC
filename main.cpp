#include <iostream>
#include <iomanip>
#include <metropolis.h>
#include <expectationvalues.h>
#include <Orbitals/orbital.h>
#include <Wavefunctions/onedimensionalslater.h>
#include <Wavefunctions/heliumwavefunction.h>
#include <Wavefunctions/slater.h>
#include <Hamiltonians/atom.h>

#include <Wavefunctions/wfTest.h>

using namespace std;
using namespace arma;




int main() {
    int N = 5;    //# of different a's
    int n = 9e5;  //# of iterations in metropolis
    int numberofpart = 2;
    bool interacting = true;

    vec avec = zeros<vec>(2);
    avec(0) = 1;    //alpha
    avec(1) = 0.35; //0.35;  //beta He: 0.35 Be: 0.098 Ne: 0.091
    Slater wave(avec,numberofpart,3,interacting); // (a, numParticls, numDims)
//    HeliumWaveFunction wave(avec,numberofpart,3,interacting);
    Atom ham(numberofpart); // Z
    ExpectationValues expect(4,&wave, &ham);
    Metropolis mysys(&expect, &wave);

    vec a = linspace(numberofpart-1,numberofpart+1,N);
    vec en;
    for(int i= 0; i<N ; i++) {
        wave.a = 1.85;//a(i);

        mysys.Run(n);



        en = expect.ev;
        cout << "a: "<< wave.a << endl;
        cout << "Variance in energy: " << en(1)-en(0)*en(0)<< endl;
        cout << setprecision(15) << "Energy: " << en(0) << endl;
        cout << "meanpos: " << en(2) << endl;
        cout << "Stddev: " << sqrt(en(3) -en(2)*en(2))<< endl << endl;
        expect.ev = zeros<vec>(4);

        cout << "snitt: " << wave.sum_snitt / n << endl;
        wave.sum_snitt = 0;
    }


    /*
     * Testing shit.
     */
//    HeliumWaveFunction wf1(avec,numberofpart,3,interacting);
//    Slater wf2(avec,numberofpart,3,interacting);
//    Atom   H(numberofpart);

//    wfTest(&wf1, &wf2, &H);

    return 0;
}






