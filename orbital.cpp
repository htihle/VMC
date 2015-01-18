#include "orbital.h"
#include <armadillo>

using namespace arma;

Orbital::Orbital() {
}


double Orbital::gradient(arma::vec x, double a) {

    int numberOfDimensions = x.size();

    for (int i = 0; i < numberOfDimensions; i ++) {

    }
    return 0;
}

double Orbital::laplacian(arma::vec x,double a) {
    int numberOfDimensions = x.size();
    double h = 0.001;
    arma::vec H = zeros<vec>(numberOfDimensions);
    double sum = 0;

    for (int i = 0; i < numberOfDimensions; i ++) {
        H(i) = h;
        sum += this->eval(x+H,a) + this->eval(x-H,a) - 2*this->eval(x,a);
        H(i) = 0;
    }
    return sum / (h*h);
}
