#include <fftw3.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <omp.h>
#include "vecMaths.h"
#include "poisson.h"

grid::grid(double gridSpacing, double dim): dim(dim), spacing(gridSpacing) {
    numPts = (int) (dim / gridSpacing);
    cout << "numPts: " << numPts << endl;
    realPot = (double*) fftw_malloc(sizeof(double) * pow(numPts, 3));
    compFFTRho = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * pow(numPts, 3));
    for (int i=0; i<3; i++) {
        comp[i] = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * pow(numPts, 3));
        realField[i] = (double *) fftw_malloc(sizeof(double) * pow(numPts, 3));
    }
    fwrd = fftw_plan_dft_r2c_3d(numPts, numPts, numPts, realPot, compFFTRho, FFTW_MEASURE);
    for (int i = 0; i < 3; i++) {
        bwrd[i] = fftw_plan_dft_c2r_3d(numPts, numPts, numPts, comp[i], realField[i], FFTW_MEASURE);
    }
}

void grid::updateGrid(vector<body>* bods){
    /// Row major format
//    realPot[int(i*pow(numPts, 2) + j*numPts + k)] = 0;
//#pragma omp parallel for //default(none) shared(bods)
    for (int i = 0; i < numPts; i++) {
        for (int j = 0; j < numPts; j++) {
            for (int k = 0; k < numPts; k++) {
                realPot[int(i*pow(numPts, 2) + j*numPts + k)] = 0;
            }
        }
    }

    for (auto & body : (*bods)){
        vector<vector<int>> mPos = meshPos(body.pos.back());
        for (auto vec : mPos) {
//            cout << "vec: " << vec[0] << ", " << vec[1] << ", " << vec[2] << endl;
//            cout << "realPos: " << realPot[int(vec[0]*pow(numPts, 2) + vec[1]*numPts + vec[2])] << endl;
            realPot[int(vec[0]*pow(numPts, 2) + vec[1]*numPts + vec[2])] +=
                    4 * pi * G * w({vec[0], vec[1], vec[2]}, body) * body.mass.back()/pow(spacing, 3);
        }
//                    cout << realPot[int(i*pow(numPts, 2) + j*numPts + k)] << endl;
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
    fftw_execute(fwrd);
    for (int axis=0; axis<3; axis++) {
//#pragma omp parallel for // default(none) shared(axis)
        for (int i = 0; i < numPts; i++) {
            for (int j = 0; j < numPts; j++) {
                for (int k = 0; k < numPts; k++) {
                    vector<double> indxVec = {double(i), double(j), double(k)};
                    fftw_complex scale;
                    if (i==0 && j==0 && k==0) {
                        scale[0] = 0;
                        scale[1] = 0;
                    } else {
                        scale[0] = 0;
                        scale[1] = indxVec[axis]/(2*pi*dim*(i*i + j*j + k*k));
                    }
                    compMultFFT(compFFTRho[int(i * pow(numPts, 2) + j * numPts + k)], scale,
                            comp[axis][int(i * pow(numPts, 2) + j * numPts + k)]);
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
        vector<double> f(3, 0);
        vector<vector<int>> mPos = meshPos(body.pos.back());
//#pragma omp parallel for
        if (mPos.size() == 8) { /// Check inside mesh
            double xd = (body.pos.back()[0] - vecPos(mPos[0][0], mPos[0][1], mPos[0][2])[0]) /
                    (vecPos(mPos[4][0], mPos[4][1], mPos[4][2])[0] -
                    vecPos(mPos[0][0], mPos[0][1], mPos[0][2])[0]);
            double yd = (body.pos.back()[1] - vecPos(mPos[0][0], mPos[0][1], mPos[0][2])[1]) /
                    (vecPos(mPos[2][0], mPos[2][1], mPos[2][2])[1] -
                    vecPos(mPos[0][0], mPos[0][1], mPos[0][2])[1]);
            double zd = (body.pos.back()[2] - vecPos(mPos[0][0], mPos[0][1], mPos[0][2])[2]) /
                        (vecPos(mPos[1][0], mPos[1][1], mPos[1][2])[2] -
                         vecPos(mPos[0][0], mPos[0][1], mPos[0][2])[2]);
            for (int axis = 0; axis < 3; axis++) {
                double C00 = realField[axis][int(mPos[0][0] * pow(numPts, 2) + mPos[0][1] * numPts + mPos[0][2])] *
                        (1 - xd) + realField[axis][int(mPos[4][0] * pow(numPts, 2) + mPos[4][1] * numPts + mPos[4][2])];
                double C01 = realField[axis][int(mPos[1][0] * pow(numPts, 2) + mPos[1][1] * numPts + mPos[1][2])] *
                             (1 - xd) + realField[axis][int(mPos[5][0] * pow(numPts, 2) + mPos[5][1] * numPts + mPos[5][2])];
                double C10 = realField[axis][int(mPos[2][0] * pow(numPts, 2) + mPos[2][1] * numPts + mPos[2][2])] *
                             (1 - xd) + realField[axis][int(mPos[6][0] * pow(numPts, 2) + mPos[6][1] * numPts + mPos[6][2])];
                double C11 = realField[axis][int(mPos[3][0] * pow(numPts, 2) + mPos[3][1] * numPts + mPos[3][2])] *
                             (1 - xd) + realField[axis][int(mPos[7][0] * pow(numPts, 2) + mPos[7][1] * numPts + mPos[7][2])];

                double C0 = C00 * (1-yd) + C10*yd;
                double C1 = C01 * (1-yd) + C11*yd;

                double C = C0 * (1-zd) + C1*zd;

                f[axis] = C;

            }

//            for (auto vec : mPos) {
//                for (int axis = 0; axis < 3; axis++) {
//                    double scaling =
//                            abs(body.pos.back()[axis] - vecPos(vec[0], vec[1], vec[2])[axis]) / (spacing * spacing);
//                    f[axis] += scaling * realField[axis][int(vec[0] * pow(numPts, 2) + vec[1] * numPts + vec[2])];
////              cout << abs(body.pos.back()[axis] - vecPos(i, j, k)[axis])  << ", " << realField[axis][int(i * pow(numPts, 2) + j * numPts + k)] << endl;
//                }
//            }
        }
        body.acc.emplace_back(f);
//      printVec(body.acc.back());
    }
}

vector<double> grid::vecPos(int i, int j, int k){
    return {i*spacing - dim/2, j*spacing - dim/2, k*spacing - dim/2};
}

vector<vector<int>> grid::meshPos(vector<double> pos) {
    vector<vector<int>> meshPoints;
    vector<int> pt(3, 0);
    for (int i = 0; i < 3; i++) {
        pt[i] = floor((pos[i] + dim/2) / spacing);
//        cout <<  "compute: " << (pos[i] + dim) << endl;
    }
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                if (pt[0]+i<numPts && pt[1]+j<numPts && pt[2]+k<numPts && pt[0]+i>0 && pt[1]+j>0 && pt[2]+k>0)
                    meshPoints.emplace_back(vecAdd(pt, {i, j, k}));
            }
        }
    }
    return meshPoints;
}

double grid::w(vector<int> vec, body& bod){
    double out = 1;
//    cout << "bodPos, vecPos: "; printVec(bod.pos.back()) ;cout << " "; printVec(vecPos(vec[0], vec[1], vec[2])); cout << endl;
    for (int i=0; i<3; i++){
        double dist = abs(bod.pos.back()[i] - vecPos(vec[0], vec[1], vec[2])[i]);
//        cout << "dist" << dist << endl;
        if (dist < spacing) {
            out *= spacing - dist/spacing;
        } else
            return 0;
    }
//    cout << "ok" << endl;

    return out;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void compMultFFT(fftw_complex v1, fftw_complex v2, fftw_complex out) {
    double real = v1[0]*v2[0] - v1[1]*v2[1];
    double imag = v1[1]*v2[0] + v1[0]*v2[1];
//    cout << v1[0] << v1[1] << endl;
    out[0] = real,
    out[1] = imag;
}

