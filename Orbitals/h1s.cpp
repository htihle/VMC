#include <h1s.h>

double H1s::eval(arma::vec x, double a, double r) {
    return exp(-a * r);
}
