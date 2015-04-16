#include <Wavefunctions/slater.h>
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

//    if (HF) {
//        // load from file.
//    } else {

//    }

    sum_snitt = 0;

    myorbital[0] = new H1s();
    myorbital[1] = new H2s();
    myorbital[2] = new H2p(0);
    myorbital[3] = new H2p(1);
    myorbital[4] = new H2p(2);
    NumberOfParticles = N;
    ParticlesInSlater = N/splitSlater;
    NumberOfDimensions = Ndim;
    R = zeros<mat>(NumberOfParticles,NumberOfParticles);
    spins = R;
    correlationsMatNew = R;
    correlationsMatOld = correlationsMatNew;
    D = 0.5;
    switch(NumberOfParticles) {
    case(2) :
        stepSize = 0.001;
        break;
    case(4) :
        stepSize = 0.0005;
        break;
    case(10) :
        stepSize = 0.0002;  //this works well with grnfnc, if not to small or too large timestep (aim for 0.999 or close) increase number of steps for small step-sizes
        break;
    default :
        stepSize = 0.0005;
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
    whichParticle = round(num2(1)*(this->NumberOfParticles-1));
    whichSlater = whichParticle / ParticlesInSlater;
    xnew.col(whichParticle) = x.col(whichParticle)
            + num*sqrt(stepSize)
            + this->quantumForceOld.col(whichParticle) * stepSize * D; //newstep(x)

    Rsd = 0;
    this->updateR(xnew,whichParticle);
    for(int j = 0;j<ParticlesInSlater;j++) // make a function get Rsd ??
    {
        Rsd += myorbital[j]->eval(xnew.col(whichParticle),a, R(whichParticle, whichParticle))*
                slaterInverse[whichSlater](j,whichParticle-whichSlater*ParticlesInSlater);
    }
    bool accept = false;
    if(interacting){
        this->updateSlaterGradient(xnew, whichParticle); //dobbelsjekk rekkefølge og logikk her
        this->updateSlaterInverse(xnew,whichParticle);
        this->updateCorrelationsMatrix(whichParticle);
        this->computeJastrowGradient(whichParticle);
        this->makeJgradient(xnew); // implement updateJgradient
        this->computeJastrowLaplacian(whichParticle);
        this->computeRc(whichParticle);
        Rcoeff = Rsd*Rc;
    } else {
        this->updateSlaterInverse(xnew, whichParticle);
        this->updateSlaterGradient(xnew, whichParticle);
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
    mat Slater = zeros<mat>(ParticlesInSlater,ParticlesInSlater);
    for(int i = 0;i<ParticlesInSlater;i++) {
        for(int j = 0; j<ParticlesInSlater;j++) {
            Slater(i,j) = myorbital[j]->eval(x.col(i+upordown*ParticlesInSlater),a,R(i+upordown*ParticlesInSlater,i+upordown*ParticlesInSlater));
        }
    }
    return Slater;
}

double Slater::slaterNumericalLaplacianLog(mat x) { // does not work with R !!!!!!!!!!!!
    double h= 0.0001;
    mat mod = zeros<mat>(NumberOfDimensions, NumberOfParticles);
    double sum = 0;
    for(int i = 0;i<NumberOfParticles;i++) {
        for (int j = 0; j<NumberOfDimensions; j++) {
            mod(j,i)+=h;
            sum += this->evaluateSlater(x+mod) +this->evaluateSlater(x-mod);
            mod(j,i)=0;
        }
    }
    sum/=this->evaluateSlater(x);
    sum -=ParticlesInSlater*4*NumberOfDimensions;
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
    for(int i=0;i<NumberOfParticles;i++)
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
    for(int k=i+1;k<NumberOfParticles;k++){
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
        for(int i =0; i<ParticlesInSlater;i++) {
            for(int j = 0;j<ParticlesInSlater;j++) {
                sum += myorbital[j]->laplacian(x.col(i+k*ParticlesInSlater),a, R(i+k*ParticlesInSlater,i+k*ParticlesInSlater))
                        *slaterInverse[k](j,i); //for some reason j and i are switched
            }
        }
    }
    return sum;
}

void Slater::updateSlaterInverse(mat x, int i)
{
    //    slaterInverse[whichSlater]= calculateSlater(x,whichSlater).i();
    i = i -whichSlater*ParticlesInSlater;

    mat oldSlater = slaterInverse[whichSlater];
    for (int k = 0;k<ParticlesInSlater; k++) {
        for (int j = 0; j<ParticlesInSlater; j++) {
            if(j != i) {
                double sum = 0;
                for( int l = 0; l< ParticlesInSlater; l++) {
                    sum += oldSlater(l,j) * myorbital[l]->eval(x.col(i+whichSlater*ParticlesInSlater), a, R(i+whichSlater*ParticlesInSlater,i+whichSlater*ParticlesInSlater));
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
    for (int i = 0; i < NumberOfParticles; i++) {
        for (int j = i+1; j < NumberOfParticles; j++) {
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
    for (int i = (k+1); i < NumberOfParticles; i++) {
        correlationsMatNew(k,i) = spins(k,i) * R(k,i) / (1 + b * R(k,i));
    }
}

/* Computes the ratio of the new correlations to the old, Rc. */
void Slater::computeRc(int i) {
    this->Rc = 0;
    for (int j = 0; j < i; j++) {
        Rc += correlationsMatNew(j,i) - correlationsMatOld(j,i);
    }
    for (int j = i+1;j<NumberOfParticles; j++) {
        Rc += correlationsMatNew(i,j) - correlationsMatOld(i,j);
    }
    Rc = exp(Rc);
}


void Slater::computeJastrowGradient(int particle) {
    for (int i = 0; i<particle; i++) {
        double factor = 1+b*R(i,particle);
        jastrowGradient(i,particle) = spins(i,particle)/(factor*factor);
    }
    for (int i = particle+1; i < NumberOfParticles; i++) {
        double factor = 1+b*R(particle,i);
        jastrowGradient(particle,i) = spins(particle,i)/(factor*factor);
    }
}

void Slater::computeJastrowLaplacian(int particle) {
    for (int i = 0; i<particle; i++) {
        double factor = 1 + b*R(i,particle);
        jastrowLaplacian(i,particle) = -2*spins(i,particle)*b/(factor*factor*factor);
    }
    for (int i = particle+1; i < NumberOfParticles; i++) {
        double factor = 1 + b*R(particle,i);
        jastrowLaplacian(particle,i) = -2*spins(particle,i)*b/(factor*factor*factor);
    }
}

double Slater::computeJastrowEnergy() { //masse wastage!!
    double sum = 0.0;
    // 1/psi_J nabla^2 psi_J
    for (int k = 0; k<NumberOfParticles; k++) {
        sum += dot(Jgradient.col(k),Jgradient.col(k));
        for (int i = 0; i<k; i++) {
            sum +=(NumberOfDimensions -1)/R(i,k)*jastrowGradient(i,k) +jastrowLaplacian(i,k);
        }
        for (int i = k+1; i < NumberOfParticles; i++) {
            sum +=(NumberOfDimensions -1)/R(k,i)*jastrowGradient(k,i) +jastrowLaplacian(k,i);
        }
    }
    // 2 * (1/psi_J nabla psi_J dot 1/psi_S nabla psi_S)
    for (int k = 0; k<NumberOfParticles; k++) {
        sum += 2*dot(Jgradient.col(k),gradient.col(k));
    }
    //cout << sum << endl;
    sum_snitt += sum;
    return sum;
}

void Slater::setUpJastrow() {

    Jgradient = zeros<mat>(NumberOfDimensions, NumberOfParticles);
    JgradientOld = zeros<mat>(NumberOfDimensions, NumberOfParticles);
    jastrowLaplacian = zeros<mat>(NumberOfParticles, NumberOfParticles);
    jastrowGradient = zeros<mat>(NumberOfParticles, NumberOfParticles);

    for (int i = 0; i < NumberOfParticles; i++) {
        this->computeJastrowGradient(i); // (R, jastrowGradient, i)
        this->computeJastrowLaplacian(i); // (R, jastrowLaplacian,i)
    }
}



//computes the gradient of the Slater part of the wavefunction
void Slater::updateSlaterGradient(mat x, int i) {
    int k = i %(ParticlesInSlater); //the remainder
    vec sum = zeros<vec>(NumberOfDimensions);
    for (int j = 0; j < (ParticlesInSlater); j++) {
        sum += slaterInverse[whichSlater](j,k) * myorbital[j]->gradient(x.col(i), a, R(i,i));
    }
    gradient.col(i) = 1.0/this->Rsd * sum;
}

void Slater::makeJgradient(mat x) {
    for (int k = 0; k < NumberOfParticles; k++) {
        for (int j =0 ; j < NumberOfDimensions; j++) {//we call the function jastrowgradient far too many times here, this should be done more effectively!!!!!


            double sum = 0;

            for (int i = 0; i < k; i++) {
                sum += (x(j,k)-x(j,i)) / R(i,k) * jastrowGradient(i,k);
            }
            for (int i = k+1; i < NumberOfParticles; i++) {
                sum -= (x(j,i)-x(j,k)) / R(k,i) * jastrowGradient(k,i);
            }

            Jgradient(j,k) = sum;
        }
    }
}

void Slater::computeSlaterGradient(arma::mat x) {
    for (int i = 0; i < NumberOfParticles; i++) {
        whichSlater = i / ParticlesInSlater;
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
    for (int j = 0; j < NumberOfParticles; j++) {
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
    vec electronSpins(NumberOfParticles);
    electronSpins.zeros();
    for (int i = 0; i < (ParticlesInSlater); i++) {
        electronSpins(i) = 1;
    }
    // Set the electron-electron spin interaction matrix.
    for (int i = 0; i < NumberOfParticles; i++) {
        for (int j = i+1; j < NumberOfParticles; j++) {
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
    x = a/ParticlesInSlater*randn<mat>(this->NumberOfDimensions, this->NumberOfParticles); //comment this out when testing (until we find a better way)
    this->makeR(x);
    this->fillSpinMatrix();
    this->fillCorrelationsMatrix();
    Rold = R;
    this->correlationsMatOld = this->correlationsMatNew;
    this->Rsd = 1;
    this->getSlaterInverse(x);
    this->gradient = zeros<mat>(this->NumberOfDimensions, 2*this->ParticlesInSlater);
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
