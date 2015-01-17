#include <iostream>
#include <metropolis.h>
#include <expectationvalues.h>
#include <ho0.h>
#include <orbital.h>
#include <onedimensionalslater.h>

using namespace std;
using namespace arma;




int main()
{
    int N = 100;    //# of different a's
    int n = 100000;  //# of iterations in metropolis
    //WaveFunction wave(0.5,3);
    OneDimensionalSlater wave(0.5,3);
    ExpectationValues expect(4,&wave);
    Metropolis mysys(&expect, &wave);
    vec a = linspace(0.5,1.5,N);
    vec en;
    for(int i= 0; i<N ; i++) {
        wave.a = a(i);
        mysys.Run(n);
        en = expect.ev;
        cout << "a: "<< wave.a << endl;
        cout << "Variance in energy: " << en(1)-en(0)*en(0)<< endl;
        cout << "Energy: " << en(0) << endl;
        cout << "meanpos: " << en(2) << endl;
        cout << "Stddev: " << sqrt(en(3) -en(2)*en(2))<< endl;
    }
    return 0;
}

