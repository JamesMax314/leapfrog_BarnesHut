#ifndef POISSON
#define POISSON

#include <bodies.h>
#include <fftw3.h>
#include "bodies.h"

using namespace std;

class grid{
    double spacing;
    int numPts;
    fftw_plan fwrd;
    fftw_plan bwrd;
public:
    double G = 6.674e-11;
    double pi = 3.14159;
    double* real;
    fftw_complex* comp;

    grid(double gridSpacing, double dim, vector<body>& bodies);

    void updateGrid(vector<body>& bods);
    void solveFiled();
    vector<double> vecPos(int i, int j, int k);
    double w(vector<int> vec, body& bod);
    void interp(vector<body>& bods);
};


#endif //POISSON