#ifndef POISSON
#define POISSON

#include <bodies.h>
#include <fftw3.h>
#include <cmath>
#include "bodies.h"

using namespace std;

class grid{
    fftw_plan fwrd{};
    fftw_plan bwrd[3]{};
public:
    double spacing{};
    vector<double> dim;
    vector<int> numPts;
    double G = 6.674e-11;
    double pi = M_PI;
    fftw_complex* realPot{};
    double* realField[3]{};
    double* realField1{};
    int* keys{};
    fftw_complex* comp[3]{};
    fftw_complex* cField[3]{};
    fftw_complex* compFFTRho{};
    fftw_complex* compFFTRho1{};
    vector<int> activeBods;

    grid(double gridSpacing, double dim);
    grid(double gridSpacing, vector<int> numPs);
    grid(grid const &g);
    grid();

    void updateGrid(vector<body>* bods);
    void solveField();
    vector<double> vecPos(int i, int j, int k);
    vector<vector<int>> meshPos(vector<double> pos);
    double w(vector<int> vec, body& bod);
    void diff(double scale);
    void interpW(vector<body>* bods, bool resetForce);
    void interp(vector<body>* bods);
    void initActiveBods(vector<body>* bods);

    void ctor(fftw_complex* arr, double* out);
    void magF();
    vector<vector<vector<double>>> getF(int indx);

    /// Pybinding getters ///
    double *data() { return realField1; }
    vector<int> size()  { return numPts; }
};

void compMultFFT(fftw_complex v1, fftw_complex v2, fftw_complex out);



#endif //POISSON