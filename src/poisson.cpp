#include <fftw3.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <omp.h>
#include "vecMaths.h"
#include "poisson.h"

grid::grid(double gridSpacing, double dim): dim(dim), spacing(gridSpacing) {
    numPts = (int) (dim / gridSpacing);
    realPot = (double*) fftw_malloc(sizeof(double) * pow(numPts, 3));
    for (int i=0; i<3; i++) {
        comp[i] = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * pow(numPts, 3));
        realField[i] = (double *) fftw_malloc(sizeof(double) * pow(numPts, 3));
    }
    for (int i = 0; i < 3; i++) {
        fwrd[i] = fftw_plan_dft_r2c_3d(numPts, numPts, numPts, realPot, comp[i], FFTW_MEASURE);
        bwrd[i] = fftw_plan_dft_c2r_3d(numPts, numPts, numPts, comp[i], realField[i], FFTW_MEASURE);
    }
}

void grid::updateGrid(vector<body>* bods){
    /// Row major format
//#pragma omp parallel for //default(none) shared(bods)
    for (int i=0; i<numPts; i++){
        for (int j=0; j<numPts; j++){
            for (int k=0; k<numPts; k++){
                realPot[int(i*pow(numPts, 2) + j*numPts + k)] = 0;
                for (auto & body : (*bods)){
                    realPot[int(i*pow(numPts, 2) + j*numPts + k)] += w({i, j, k}, body) * body.mass.back();
//                    cout << realPot[int(i*pow(numPts, 2) + j*numPts + k)] << endl;
                }
            }
        }
    }
//    cout << endl;
//    cout << "[ ";
//    for (int i=0; i<numPts; i++){
//        cout << "[ ";
//        for (int j=0; j<numPts; j++){
//            cout << "[ ";
//            for (int k=0; k<numPts; k++){
//                cout << realPot[int(i*pow(numPts, 2) + j*numPts + k)] << ", ";
//            }
//            cout << " ], ";
//        }
//        cout << " ], ";
//    }
//    cout << " ]";
//    cout << endl;
}

void grid::solveField(){
    for (int axis=0; axis<3; axis++) {
        fftw_execute(fwrd[axis]);
//#pragma omp parallel for // default(none) shared(axis)
        for (int i = 0; i < numPts; i++) {
            for (int j = 0; j < numPts; j++) {
                for (int k = 0; k < numPts; k++) {
                    vector<double> kVec = vecPos(i, j, k);
                    for (int l = 0; l < 2; l++) {
                        comp[axis][int(i * pow(numPts, 2) + j * numPts + k)][l] *= -4 * pi * G * kVec[axis] / spacing;
//                        cout << comp[axis][int(i * pow(numPts, 2) + j * numPts + k)][l] << endl;
                    }
                }
            }
        }
        fftw_execute(bwrd[axis]);
    }
//    cout << endl;
//    cout << "[ ";
//    for (int i=0; i<numPts; i++){
//        cout << "[ ";
//        for (int j=0; j<numPts; j++){
//            cout << "[ ";
//            for (int k=0; k<numPts; k++){
//                cout << comp[0][int(i*pow(numPts, 2) + j*numPts + k)][0] << " + i" << comp[0][int(i*pow(numPts, 2) + j*numPts + k)][1] << ", ";
//            }
//            cout << " ], ";
//        }
//        cout << " ], ";
//    }
//    cout << " ]";
//    cout << endl;
//
//    cout << endl;
//    cout << "[ ";
//    for (int i=0; i<numPts; i++){
//        cout << "[ ";
//        for (int j=0; j<numPts; j++){
//            cout << "[ ";
//            for (int k=0; k<numPts; k++){
//                cout << realField[0][int(i*pow(numPts, 2) + j*numPts + k)] << ", ";
//            }
//            cout << " ], ";
//        }
//        cout << " ], ";
//    }
//    cout << " ]";
//    cout << endl;
}

void grid::interp(vector<body>* bods){
    for (auto & body : (*bods)){
        int count = 0;
        vector<double> f(3, 0);
//#pragma omp parallel for
        for (int i=0; i<numPts; i++){
            for (int j=0; j<numPts; j++){
                for (int k=0; k<numPts; k++){
                    bool inCube = true;
                    for (int axis=0; axis<3; axis++){
                        if (abs(body.pos.back()[axis] - vecPos(i, j, k)[axis]) > spacing) {
                            inCube = false;
//                            cout << body.pos.back()[axis] << endl;
                        }
                    }
                    if (inCube){
//                        cout << "in cube" << endl;
                        count += 1;
                        for (int axis=0; axis<3; axis++) {
                            double scaling = abs(body.pos.back()[axis] - vecPos(i, j, k)[axis])/spacing;
                            f[axis] += scaling * realField[axis][int(i * pow(numPts, 2) + j * numPts + k)] * 1e20;
//                            cout << abs(body.pos.back()[axis] - vecPos(i, j, k)[axis])  << ", " << realField[axis][int(i * pow(numPts, 2) + j * numPts + k)] << endl;
                        }
                    }
                    if (count == 8)
                        i = numPts;
                        j = numPts;
                        k = numPts;
                }
            }
        }
        body.acc.emplace_back(f);
//        printVec(body.acc.back());
    }
}

vector<double> grid::vecPos(int i, int j, int k){
    return {i*spacing - dim/2, j*spacing - dim/2, k*spacing - dim/2};
}



double grid::w(vector<int> vec, body& bod){
    double out = 0;
//    cout << "bodPos, vecPos: "; printVec(bod.pos.back()) ;cout << " "; printVec(vecPos(vec[0], vec[1], vec[2])); cout << endl;
    for (int i=0; i<3; i++){
        double dist = abs(bod.pos.back()[i] - vecPos(vec[0], vec[1], vec[2])[i]);
//        cout << "dist" << dist << endl;
        if (dist < spacing) {
            out += spacing - dist;
        } else
            return 0;
    }
//    cout << "ok" << endl;

    return out;
}

