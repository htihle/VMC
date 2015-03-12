#include <iostream>
#include <iomanip>
#include <metropolis.h>
#include <expectationvalues.h>
#include <orbital.h>
#include <onedimensionalslater.h>
#include <heliumwavefunction.h>
#include <slater.h>
#include <atom.h>

using namespace std;
using namespace arma;




int main()
{


    /*
     *
     *
     *
     *        \   /
     *        .   .
     *          v
     *      .        .
     *       .      .
     *         . . .
     *
     *
     * HÅVARD TROR DET FUNKER NÅ UTEN JASTROW MED SLATER SHIT. LOL, HAN ER NUB.
     *
     *
     *
     */

    int N = 15;    //# of different a's
    int n = 1e5;  //# of iterations in metropolis
    int numberofpart = 2;
    bool interacting = false;

    vec avec = zeros<vec>(2);
    avec(0) = 1;    //alpha
    avec(1) = 1e50; //0.35;  //beta He: 0.35 Be: 0.098 Ne: 0.091
    //Slater wave(avec,numberofpart,3,interacting); // (a, numParticls, numDims)
    HeliumWaveFunction wave(avec,numberofpart,3,interacting);
    Atom ham(numberofpart); // Z
    ExpectationValues expect(4,&wave, &ham);
    Metropolis mysys(&expect, &wave);

    vec a = linspace(numberofpart-1,numberofpart,N);
    vec en;

    for(int i= 0; i<N ; i++) {
        wave.a = a(i);

        mysys.Run(n);


        en = expect.ev;
        cout << "a: "<< wave.a << endl;
        cout << "Variance in energy: " << en(1)-en(0)*en(0)<< endl;
        cout << setprecision(15) << "Energy: " << en(0) << endl;
        cout << "meanpos: " << en(2) << endl;
        cout << "Stddev: " << sqrt(en(3) -en(2)*en(2))<< endl << endl;
    }
    return 0;
}

