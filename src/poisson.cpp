#include <complex>
#include <fftw3.h>
#include <vector>
#include <cmath>
#include "poisson.h"

grid::grid(double gridSpacing, double dim, vector<body>& bodies): spacing(gridSpacing) {
    numPts = (int)(dim/gridSpacing);
    fwrd = fftw_plan_dft_r2c_3d(numPts,numPts,numPts,real,comp,FFTW_MEASURE);
    bwrd = fftw_plan_dft_c2r_3d(numPts,numPts,numPts,comp,real,FFTW_MEASURE);
}

void grid::updateGrid(vector<body>& bods){
    real = (double*) fftw_malloc(sizeof(double) * pow(numPts, 3));

    /// Row major format
    for (int i=0; i<numPts; i++){
        for (int j=0; j<numPts; j++){
            for (int k=0; k<numPts; k++){
                real[int(i*pow(numPts, 2) + j*numPts + k)] = 0;
                for (auto body : bods){
                    real[int(i*pow(numPts, 2) + j*numPts + k)] = w({i, j, k}, body);
                }
            }
        }
    }
}

void grid::solveFiled(){
    fftw_execute(fwrd);
    for (int i=0; i<numPts; i++){
        for (int j=0; j<numPts; j++){
            for (int k=0; k<numPts; k++){
                for (int l=0; l<2; l++)
                    comp[int(i*pow(numPts, 2) + j*numPts + k)][l] *= -4*pi*G;
            }
        }
    }
    fftw_execute(bwrd);
}

void grid::interp(vector<body>& bods){
    for (auto body : bods){
        int count = 0;
        vector
        for (int i=0; i<numPts; i++){
            for (int j=0; j<numPts; j++){
                for (int k=0; k<numPts; k++){
                    bool inCube = true;
                    for (int l=0; l<2; l++){
                        if (abs(body.pos.back()[i] - vecPos(i, j, k)[l]))
                            inCube = false;
                    }
                    if (inCube){

                    }
                        real[int(i*pow(numPts, 2) + j*numPts + k)];
                }
            }
        }
    }
}

vector<double> grid::vecPos(int i, int j, int k){
    return {i*spacing, j*spacing, k*spacing};
}


double grid::w(vector<int> vec, body& bod){
    double out = 0;
    for (int i=0; i<3; i++){
        double dist = abs(bod.pos.back()[i] - vecPos(vec[0], vec[1], vec[2])[i]);
        if (dist < 1)
            out += 1-dist;
        else
            return 0;
    }
    return out;
}