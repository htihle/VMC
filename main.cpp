#include <iostream>
#include <iomanip>
#include <metropolis.h>
#include <expectationvalues.h>
#include <orbital.h>
#include <onedimensionalslater.h>
#include <slater.h>
#include <atom.h>

using namespace std;
using namespace arma;




int main()
{
    int N = 15;    //# of different a's
    int n = 1e6;  //# of iterations in metropolis
    int numberofpart = 2;
    vec avec = zeros<vec>(2);
    avec(0) = 1;    //alpha
    avec(1) = 0.5;  //beta
    Slater wave(avec,numberofpart,3); // (a, numParticls, numDims)
    Atom ham(numberofpart); // Z
    ExpectationValues expect(4,&wave, &ham);
    Metropolis mysys(&expect, &wave);

    vec a = linspace(numberofpart-0.3,numberofpart,N);
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

