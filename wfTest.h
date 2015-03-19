#pragma once
#include <wavefunction.h>
#include <hamiltonian.h>
#include <armadillo>
#include <heliumwavefunction.h>
#include <slater.h>
#include <atom.h>

void wfTest(HeliumWaveFunction *wf1, Slater *wf2, Atom *H);
