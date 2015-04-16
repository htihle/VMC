#pragma once
#include <Wavefunctions/wavefunction.h>
#include <Hamiltonians/hamiltonian.h>
#include <armadillo>
#include <Wavefunctions/heliumwavefunction.h>
#include <Wavefunctions/slater.h>
#include <Hamiltonians/atom.h>

void wfTest(HeliumWaveFunction *wf1, Slater *wf2, Atom *H);
