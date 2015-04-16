#include "harmonicoscillator.h"

using namespace arma;

HarmonicOscillator::HarmonicOscillator(int nParticles, double omega) {
    this->nParticles  = nParticles;
    this->omega       = omega;
}


double HarmonicOscillator::calculatePotential(arma::mat x) {
    double sum = 0;

    for (int i = 0; i < nParticles; i++) {
        double r = norm(x.row(i));
        sum +=0.5*omega * r * r;
    }

    return sum;
}
