#include "h2p.h"

H2p::H2p(int state)
{
    this->state = state;
}

double H2p::eval(arma::vec x, double a,double r) {
    return x(state)*exp(-a * r/2);
}
