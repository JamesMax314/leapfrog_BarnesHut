#ifndef POISSON
#define POISSON

#include <bodies.h>
#include <fftw3.h>
#include "bodies.h"

using namespace std;

class grid{
    double spacing;
    double dim;
    fftw_plan fwrd;
    fftw_plan bwrd; //[3];
public:
    int numPts;
    double G = 6.674e-11;
    double pi = 3.14159;
    double* realPot;
    double* realField[3];
    double* realField1;
    fftw_complex* comp[3];
    fftw_complex* compFFTRho;
    fftw_complex* compFFTRho1;

    grid(double gridSpacing, double dim);

    void updateGrid(vector<body>* bods);
    void solveField();
    vector<double> vecPos(int i, int j, int k);
    vector<vector<int>> meshPos(vector<double> pos);
    double w(vector<int> vec, body& bod);
    void diff(double scale);
    void interp(vector<body>* bods);

    void ctor();

    /// Pybinding getters ///
    double *data() { return realPot; }
    size_t size()  { return numPts; }
};

void compMultFFT(fftw_complex v1, fftw_complex v2, fftw_complex out);



#endif //POISSON