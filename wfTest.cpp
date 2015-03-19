#include <wfTest.h>

using namespace std;
using namespace arma;


void wfTest(HeliumWaveFunction *wf1, Slater *wf2, Atom *H) {

    int numberOfParticles = wf1->NumberOfParticles;
    int numberOfDimensions = wf1->NumberOfDimensions;

    mat x = randn<mat>(numberOfParticles, numberOfDimensions);

    wf1->setUpForMetropolis(x);
    wf2->setUpForMetropolis(x);

    // Testing.
    cout << "\n\n\n                         **************************\n";
    cout <<       "                         ***   Start testing:   ***\n";
    cout <<       "                         **************************\n\n";

    cout << " * Test 1: Laplacian log:\n";
    cout << "-----------------------------\n";
    double lap1 = wf1->laplacianLog(x);
    double lap2 = wf2->laplacianLog(x);
    cout << "   - lapLog1: " << lap1 << endl;
    cout << "   - lapLog2: " << lap2 << endl;
    cout << "                                        Difference = " << fabs(lap1-lap2) << endl << endl;


    cout << " * Test 2: Quantum force:\n";
    cout << "-----------------------------\n";
    mat qf1 = wf1->QF(1,x).t();
    mat qf2 = wf2->quantumForceNew.col(1).t();
    cout << "   - qf1: " << qf1;
    cout << "   - qf2: " << qf2;
    cout << "                                        Difference = " << norm(qf1-qf2) << endl << endl;


    cout << " * Test 2: Wave function:\n";
    cout << "-----------------------------\n";
    double wave1 = wf1->wf();
    double wave2 = det(wf2->calculateSlater(x,0))*det(wf2->calculateSlater(x,1));
    cout << "   - wf1: " << wave1 << endl;
    cout << "   - wf2: " << wave2 << endl;
    cout << "                                        Difference = " << fabs(wave1-wave2) << endl << endl;

}
