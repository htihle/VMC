#include <h1s.h>

double H1s::eval(arma::vec x, double a) {
    return exp(-a * arma::norm(x));
}
