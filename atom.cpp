#include "atom.h"

using namespace arma;
Atom::Atom(double Z)
{
    this->Z = Z;
}

double Atom::calculatePotential(mat x)
{
    int NumberOfParticles = x.n_cols;
    double sum = 0;
    for(int i = 0;i<NumberOfParticles;i++)
    {
        sum -= Z/ (norm(x.col(i)));
        for(int j= i+1;j<NumberOfParticles;j++)
        {
            sum += 1/ (norm(x.col(i)-x.col(j)));
        }
    }
    return sum;
}
