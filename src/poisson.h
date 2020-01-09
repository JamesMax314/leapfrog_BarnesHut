#ifndef POISSON
#define POISSON

#include <bodies.h>
#include <fftw3.h>
#include "bodies.h"

using namespace std;

class grid{
    double spacing;
    double dim;
    int numPts;
    fftw_plan fwrd[3];
    fftw_plan bwrd[3];
public:
    double G = 6.674e-11;
    double pi = 3.14159;
    double* realPot;
    double* realField[3];
    fftw_complex* comp[3];

    grid(double gridSpacing, double dim);

    void updateGrid(vector<body>* bods);
    void solveField();
    vector<double> vecPos(int i, int j, int k);
    double w(vector<int> vec, body& bod);
    void interp(vector<body>* bods);
};


#endif //POISSON