#include "h2s.h"

double H2s::eval(arma::vec x, double a, double r) {
    return (1-a*r/2)*exp(-a/2 * r);
}
