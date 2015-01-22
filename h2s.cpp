#include "h2s.h"

double H2s::eval(arma::vec x, double a) {
    double r = arma::norm(x);
    return (1-a*r/2)*exp(-a/2 * r);
}
