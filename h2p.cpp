#include "h2p.h"

H2p::H2p(int state)
{
    this->state = state;
}

double H2p::eval(arma::vec x, double a) {
    return x(state)*exp(-a * arma::norm(x)/2);
}
