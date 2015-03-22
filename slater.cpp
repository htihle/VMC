#include <slater.h>
#include <armadillo>

using namespace std;

using arma::mat;
using arma::vec;
using arma::randn;
using arma::randu;
using arma::zeros;


Slater::Slater(vec a, int N, int Ndim, bool interacting) : WaveFunction(a,N,Ndim) {
    this->interacting = interacting;
    this->a = a(0);
    this->b = a(1);
    this->whichSlater = 0;
    this->splitSlater = 2;
    //    myorbital[0] = new HO0();
    //    myorbital[1] = new HO1();
    //    myorbital[2] = new HO2();
    myorbital[0] = new H1s();
    myorbital[1] = new H2s();
    myorbital[2] = new H2p(0);
    myorbital[3] = new H2p(1);
    myorbital[4] = new H2p(2);
    NumberOfParticles = N/2;
    NumberOfDimensions = Ndim;
    R = zeros<mat>(NumberOfParticles*2,NumberOfParticles*2);
    spins = R;
    correlationsMatNew = R;
    correlationsMatOld = correlationsMatNew;
    D = 0.5;
    switch(NumberOfParticles*2) {
    case(2) :
        stepSize = 0.001;
        break;
    case(4) :
        stepSize = 0.001;   //nan when time step too large ?? Why? figure this out and you have found the problem
        break;
    case(10) :
        stepSize = 0.001;
        break;
    default :
        break;
    }
}

double Slater::laplacianLog(mat x)
{
    double lapS = this->slaterAnalyticalLaplacianLog(x);
    //double lap = this->slaterNumericalLaplacianLog(x); // does not work!!
    if(interacting) {
        return lapS + this->computeJastrowEnergy();
    } else {
        return lapS;
    }
}

bool Slater::newStep(mat &xnew, mat x,int &whichParticle) {

    double Rcoeff;
    xnew = x;
    vec num  = randn<vec>(NumberOfDimensions);
    vec num2 = randu<vec>(2);
    whichParticle = round(num2(1)*(this->NumberOfParticles*2-1));
    whichSlater = whichParticle / NumberOfParticles;
    //xnew.col(whichParticle) = x.col(whichParticle) + this->stepSize/NumberOfParticles/NumberOfDimensions*a*num; //newstep(x)

    xnew.col(whichParticle) = x.col(whichParticle)
            + num*sqrt(stepSize)
            + this->quantumForceOld.col(whichParticle) * stepSize * D; //newstep(x)
    Rsd = 0;
    this->updateR(xnew,whichParticle);
    for(int j = 0;j<NumberOfParticles;j++) // make a function get Rsd ??
    {
        Rsd += myorbital[j]->eval(xnew.col(whichParticle),a, R(whichParticle, whichParticle))*
                slaterInverse[whichSlater](j,whichParticle-whichSlater*NumberOfParticles);
    }
    bool accept = false;
    if(interacting){
        this->updateSlaterInverse(x,whichParticle);
        this->updateCorrelationsMatrix(whichParticle);
        this->computeJastrowGradient(whichParticle);
        this->makeJgradient(x); // implement updateJgradient
        this->computeJastrowLaplacian(whichParticle);
        this->computeRc(whichParticle);
        this->updateSlaterGradient(x, whichParticle);
        Rcoeff = Rsd*Rc;
    } else {
        this->updateSlaterInverse(x, whichParticle);
        this->updateSlaterGradient(x, whichParticle);
        Rcoeff = Rsd;
    }

    double coeff = Rcoeff*Rcoeff* this->computeGreensFunction(xnew, x);

    if(coeff > num2(0)){
        accept = true;

        Rold               = R;
        correlationsMatOld = correlationsMatNew;
        gradientOld        = gradient;
        JgradientOld       = Jgradient;
        quantumForceOld    = quantumForceNew;
        acceptanceCounter +=1;
    } else {
        R                  = Rold;
        correlationsMatNew = correlationsMatOld;
        gradient           = gradientOld;
        Jgradient          = JgradientOld;
        quantumForceNew    = quantumForceOld;
    }
    if (this->interacting) {
        this->quantumForceOld = 2 * (this->gradientOld + this->JgradientOld);
        this->quantumForceNew = 2 * (this->gradient    + this->Jgradient);
    } else {
        this->quantumForceOld = 2 * (this->gradientOld);
        this->quantumForceNew = 2 * (this->gradient);
    }
    return accept;
}

double Slater::evaluateSlater(mat x) {
    return det(calculateSlater(x,0))*det(calculateSlater(x,1));   //   fix for general (non-split) slater!!
}

mat Slater::calculateSlater(mat x,int upordown) {
    mat Slater = zeros<mat>(NumberOfParticles,NumberOfParticles);
    for(int i = 0;i<NumberOfParticles;i++) {
        for(int j = 0; j<NumberOfParticles;j++) {

            Slater(i,j) = myorbital[j]->eval(x.col(i+upordown*NumberOfParticles),a,R(i+upordown*NumberOfParticles,i+upordown*NumberOfParticles));
        }
    }
    return Slater;
}

double Slater::slaterNumericalLaplacianLog(mat x) { // does not work with R !!!!!!!!!!!!
    double h= 0.0001;
    mat mod = zeros<mat>(NumberOfDimensions, NumberOfParticles*2);
    double sum = 0;
    for(int i = 0;i<NumberOfParticles*2;i++) {
        for (int j = 0; j<NumberOfDimensions; j++) {
            mod(j,i)+=h;
            sum += this->evaluateSlater(x+mod) +this->evaluateSlater(x-mod);
            mod(j,i)=0;
        }
    }
    sum/=this->evaluateSlater(x);
    sum -=NumberOfParticles*4*NumberOfDimensions;
    sum/=h*h;
    return sum;
}
void Slater::getSlaterInverse(mat x)
{
    for (int i = 0; i < splitSlater; i++) {
        slaterInverse[i] = this->calculateSlater(x,i);
        slaterInverse[i] = inv(slaterInverse[i]);
    }
}

void Slater::makeR(mat x)
{
    for(int i=0;i<NumberOfParticles*2;i++)
    {
        this->updateR(x,i);
    }
}

// Updates the distance matrix when we have just changed one coordinate of particle "i"
// (like we do when computing the derivative)*/
void Slater::updateR(mat x, int i){
    double sum;
    for(int k=0; k<i; k++) {
        sum = norm(x.col(i)- x.col(k));
        R(k,i) = sum; //R is the matrix of distances
    }
    for(int k=i+1;k<NumberOfParticles*2;k++){
        sum = norm(x.col(i)- x.col(k));
        R(i,k) = sum; //R is the matrix of distances
    }
    sum = norm(x.col(i));
    R(i,i) = sum;
    //    cout << norm(x.col(0))- R(0,0) << endl;
}

double Slater::slaterAnalyticalLaplacianLog(mat x) {
    // we should keep the laplacian wrt all the different particles separately and update only
    // the ones we have moved (or is that correct?)
    double sum = 0;
    for (int k = 0; k < splitSlater; k++) {
        for(int i =0; i<NumberOfParticles;i++) {
            for(int j = 0;j<NumberOfParticles;j++) {
                //                double p = norm(x.col(i+k*NumberOfParticles)) - R(i+k*NumberOfParticles,i+k*NumberOfParticles);
                //                if (abs(p) > 0.00001) {
                //                    cout << x << R << endl;
                //                    cout << p << endl;
                //                    cout << i+k*NumberOfParticles << endl;
                //                }

                sum += myorbital[j]->laplacian(x.col(i+k*NumberOfParticles),a, R(i+k*NumberOfParticles,i+k*NumberOfParticles))
                        *slaterInverse[k](j,i); //for some reason j and i are switched
            }
        }
    }
    return sum;
}

void Slater::updateSlaterInverse(mat x, int i)
{
    //    slaterInverse[whichSlater]= calculateSlater(x,whichSlater).i();
    i = i -whichSlater*NumberOfParticles;

    mat oldSlater = slaterInverse[whichSlater];
    for (int k = 0;k<NumberOfParticles; k++) {
        for (int j = 0; j<NumberOfParticles; j++) {
            if(j != i) {
                double sum = 0;
                for( int l = 0; l< NumberOfParticles; l++) {
                    sum += oldSlater(l,j) * myorbital[l]->eval(x.col(i+whichSlater*NumberOfParticles), a, R(i+whichSlater*NumberOfParticles,i+whichSlater*NumberOfParticles));
                }
                slaterInverse[whichSlater](k,j) = oldSlater(k,j) - oldSlater(k,i) * sum / Rsd;
            }
            else {
                slaterInverse[whichSlater](k,j) = oldSlater(k,i) / Rsd;
            }
        }
    }
}

/* Compute all the correlation factors. */
void Slater::fillCorrelationsMatrix() {
    for (int i = 0; i < NumberOfParticles*2; i++) {
        for (int j = i+1; j < NumberOfParticles*2; j++) {
            correlationsMatNew(i,j) = spins(i,j) * R(i,j) / (1 + b * R(i,j));
        }
    }
}

/* Updates the new correlations matrix. */
void Slater::updateCorrelationsMatrix(int particle) {
    int k = particle;
    for (int i = 0; i < k; i++) {
        correlationsMatNew(i,k) = spins(i,k) * R(i,k) / (1 + b * R(i,k));
    }
    for (int i = (k+1); i < NumberOfParticles*2; i++) {
        correlationsMatNew(k,i) = spins(k,i) * R(k,i) / (1 + b * R(k,i));
    }
}

/* Computes the ratio of the new correlations to the old, Rc. */
void Slater::computeRc(int i) {
    this->Rc = 0;
    for (int j = 0; j < i; j++) {
        Rc += correlationsMatNew(j,i) - correlationsMatOld(j,i);
    }
    for (int j = i+1;j<NumberOfParticles*2; j++) {
        Rc += correlationsMatNew(i,j) - correlationsMatOld(i,j);
    }
    Rc = exp(Rc);
}


void Slater::computeJastrowGradient(int particle) {
    for (int i = 0; i<particle; i++) {
        double factor = 1+b*R(i,particle);
        jastrowGradient(i,particle) = spins(i,particle)/(factor*factor);
    }
    for (int i = particle+1; i < NumberOfParticles*2; i++) {
        double factor = 1+b*R(particle,i);
        jastrowGradient(particle,i) = spins(particle,i)/(factor*factor);
    }
}

void Slater::computeJastrowLaplacian(int particle) {
    for (int i = 0; i<particle; i++) {
        double factor = 1 + b*R(i,particle);
        jastrowLaplacian(i,particle) = -2*spins(i,particle)*b/(factor*factor*factor);
    }
    for (int i = particle+1; i < NumberOfParticles*2; i++) {
        double factor = 1 + b*R(particle,i);
        jastrowLaplacian(particle,i) = -2*spins(particle,i)*b/(factor*factor*factor);
    }
}

double Slater::computeJastrowEnergy() { //masse wastage!!
    double sum = 0.0;
    // 1/psi_J nabla^2 psi_J
    for (int k = 0; k<NumberOfParticles*2; k++) {
        sum += dot(Jgradient.col(k),Jgradient.col(k));
        for (int i = 0; i<k; i++) {
            sum +=(NumberOfDimensions -1)/R(i,k)*jastrowGradient(i,k) +jastrowLaplacian(i,k);
        }
        for (int i = k+1; i < NumberOfParticles*2; i++) {
            sum +=(NumberOfDimensions -1)/R(k,i)*jastrowGradient(k,i) +jastrowLaplacian(k,i);
        }
    }
    // 2 * (1/psi_J nabla psi_J dot 1/psi_S nabla psi_S)
    for (int k = 0; k<NumberOfParticles*2; k++) {
        if(Jgradient.col(k).n_elem != NumberOfDimensions){
            cout << "Noe er galt!!" << endl;
        }

        sum += 2*dot(Jgradient.col(k),gradient.col(k));
    }
    return sum;
}

void Slater::setUpJastrow() {

    Jgradient = zeros<mat>(NumberOfDimensions, NumberOfParticles*2);
    JgradientOld = zeros<mat>(NumberOfDimensions, NumberOfParticles*2);
    jastrowLaplacian = zeros<mat>(NumberOfParticles*2, NumberOfParticles*2);
    jastrowGradient = zeros<mat>(NumberOfParticles*2, NumberOfParticles*2);

    for (int i = 0; i < NumberOfParticles*2; i++) {
        this->computeJastrowGradient(i); // (R, jastrowGradient, i)
        this->computeJastrowLaplacian(i); // (R, jastrowLaplacian,i)
    }
}



//computes the gradient of the Slater part of the wavefunction
void Slater::updateSlaterGradient(mat x, int i) {
    int k = i %(NumberOfParticles); //the remainder
    vec sum = zeros<vec>(NumberOfDimensions);
    for (int j = 0; j < (NumberOfParticles); j++) {
        sum += slaterInverse[whichSlater](j,k) * myorbital[j]->gradient(x.col(i), a, R(i,i));
    }
    gradient.col(i) = 1.0/this->Rsd * sum;
}

void Slater::makeJgradient(mat x) {
    for (int k = 0; k < NumberOfParticles*2; k++) {
        for (int j =0 ; j < NumberOfDimensions; j++) {//we call the function jastrowgradient far too many times here, this should be done more effectively!!!!!


            double sum = 0;

            for (int i = 0; i < k; i++) {
                sum += (x(j,k)-x(j,i)) / R(i,k) * jastrowGradient(i,k);
            }
            for (int i = k+1; i < NumberOfParticles*2; i++) {
                sum -= (x(j,i)-x(j,k)) / R(k,i) * jastrowGradient(k,i);
            }

            Jgradient(j,k) = sum;
        }
    }
}

void Slater::computeSlaterGradient(arma::mat x) {
    for (int i = 0; i < NumberOfParticles*2; i++) {
        whichSlater = i / NumberOfParticles;
        updateSlaterGradient(x, i);
    }
}


// Maybe make this more efficient (some day...)
double Slater::computeGreensFunction(arma::mat x, arma::mat xOld) {
    this->makeJgradient(x);

    if (this->interacting) {
        this->quantumForceOld = 2 * (this->gradientOld + this->JgradientOld);
        this->quantumForceNew = 2 * (this->gradient    + this->Jgradient);
    } else {
        this->quantumForceOld = 2 * (this->gradientOld);
        this->quantumForceNew = 2 * (this->gradient);
    }

    double greensFunction = 0;
    for (int j = 0; j < NumberOfParticles*2; j++) {
        for (int i = 0; i < NumberOfDimensions; i++) {
            greensFunction += 0.5* (quantumForceOld(i,j) +
                                    quantumForceNew(i,j)) *
                    (D * this->stepSize * 0.5* (quantumForceOld(i,j) -
                                                quantumForceNew(i,j)) -
                     x(i,j) + xOld(i,j));
        }
    }
    return exp(greensFunction);
}


/*
* Fills the matrix of spin values, which is used to compute the correlation
* part of the wavefunction.
*/
void Slater::fillSpinMatrix() {
    // Set which electrons have spins up, and which have spins down.
    vec electronSpins(NumberOfParticles*2);
    electronSpins.zeros();
    for (int i = 0; i < (NumberOfParticles); i++) {
        electronSpins(i) = 1;
    }
    // Set the electron-electron spin interaction matrix.
    for (int i = 0; i < NumberOfParticles*2; i++) {
        for (int j = i+1; j < NumberOfParticles*2; j++) {
            if (electronSpins(i) != electronSpins(j)) {
                spins(i,j) = 0.5;
            } else {
                spins(i,j) = 0.25;
            }
        }
    }
}


void Slater::setUpForMetropolis(mat &x) {
    acceptanceCounter = 0;
    arma::arma_rng::set_seed(200); // not sure if this helps
    x = a/NumberOfParticles*randn<mat>(this->NumberOfDimensions, this->NumberOfParticles*2); //comment this out when testing (until we find a better way)
    this->makeR(x);
    this->fillSpinMatrix();
    this->fillCorrelationsMatrix();
    Rold = R;
    this->correlationsMatOld = this->correlationsMatNew;
    this->Rsd = 1;
    this->getSlaterInverse(x);
    this->gradient = zeros<mat>(this->NumberOfDimensions, 2*this->NumberOfParticles);
    this->computeSlaterGradient(x);
    this->setUpJastrow();
    this->makeJgradient(x);
    if(interacting) {
        this->quantumForceNew = 2*(gradient + Jgradient);
    } else {
        this->quantumForceNew = 2*gradient;
    }
    this->quantumForceOld = quantumForceNew;
    this->JgradientOld = Jgradient;
    this->gradientOld = gradient;
}
