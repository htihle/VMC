#include "orbital.h"
#include <armadillo>

using namespace arma;

Orbital::Orbital() {
}


vec Orbital::gradient(arma::vec x, double a, double r) {

    int numberOfDimensions = x.size();
    double h = 0.00001;
    arma::vec H = zeros<vec>(numberOfDimensions);
    vec grad = zeros<vec>(numberOfDimensions);
    for (int i = 0; i < numberOfDimensions; i ++) {
        H(i) = h;
        grad(i)= (this->eval(x+H,a,norm(x+H)) - this->eval(x-H,a,norm(x-H)))/(2*h);
        H(i) = 0;
    }
    return grad;
}

double Orbital::laplacian(arma::vec x,double a, double r) {
    int numberOfDimensions = x.size();
    double h = 0.00001;
    arma::vec H = zeros<vec>(numberOfDimensions);
    double sum = 0;

    for (int i = 0; i < numberOfDimensions; i ++) {
        H(i) = h;
        sum += this->eval(x+H,a,norm(x+H)) + this->eval(x-H,a,norm(x-H)) - 2*this->eval(x,a,r);
        H(i) = 0;
    }
    return sum / (h*h);
}
