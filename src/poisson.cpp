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
    realPot = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * pow(numPts, 3));
    realField1 = (double *) fftw_malloc(sizeof(double) * pow(numPts, 3));
    compFFTRho = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * pow(numPts, 3));
    compFFTRho1 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * pow(numPts, 3));
    for (auto & i : realField) {
        i = new double[int(pow(numPts, 3))];
    }
    for (int i=0; i<3; i++) {
        comp[i] = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * pow(numPts, 3));
        cField[i] = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * pow(numPts, 3));
//        realField[i] = (double *) fftw_malloc(sizeof(double) * pow(numPts, 3));
    }
    fwrd = fftw_plan_dft_3d(numPts, numPts, numPts, realPot, compFFTRho, FFTW_FORWARD, FFTW_MEASURE);
//    bwrd = fftw_plan_dft_3d(numPts, numPts, numPts, compFFTRho1, realPot, FFTW_BACKWARD, FFTW_MEASURE);
    for (int i = 0; i < 3; i++) {
        bwrd[i] = fftw_plan_dft_3d(numPts, numPts, numPts, comp[i], cField[i], FFTW_BACKWARD,FFTW_MEASURE);
    }
}

void grid::updateGrid(vector<body>* bods){
    /// Row major format
//    realPot[int(i*pow(numPts, 2) + j*numPts + k)] = 0;
//#pragma omp parallel for //default(none) shared(bods)
    for (int i = 0; i < numPts; i++) {
        for (int j = 0; j < numPts; j++) {
            for (int k = 0; k < numPts; k++) {
                realPot[int(i*pow(numPts, 2) + j*numPts + k)][0] = 0;
            }
        }
    }

    for (auto & body : (*bods)){
        vector<vector<int>> mPos = meshPos(body.pos.back());
        for (auto vec : mPos) {
            realPot[int(vec[0]*pow(numPts, 2) + vec[1]*numPts + vec[2])][0] +=
                    w({vec[0], vec[1], vec[2]}, body) * body.mass.back(); // /pow(spacing, 3);
        }
    }
}

void grid::solveField(){
    fftw_execute(fwrd);
//#pragma omp parallel for // default(none) shared(axis)
    for (int axis=0; axis<3; axis++) {
        for (int i = 0; i < numPts; i++) { //floor(numPts/2) + 1
            for (int j = 0; j < numPts; j++) {
                for (int k = 0; k < numPts; k++) {
                    /// Accounting for weird fft structure
                    int kx = (i <= numPts / 2) ? i : i - numPts;
                    int ky = (j <= numPts / 2) ? j : j - numPts;
                    int kz = (k <= numPts / 2) ? k : k - numPts;
                    vector<int> ks = {kx, ky, kz};

                    fftw_complex scale;
                    fftw_complex one;

                    one[0] = 1;
                    one[1] = 1;

                    scale[1] = (i == 0 && j == 0 && k == 0) ? 0 : dim * ks[axis] * 4 * pi * G /
                            (pow(numPts, 3) * (kx * kx + ky * ky + kz * kz)); //dim*4*pi*G
                    scale[0] = 0;

    //                compMultFFT(compFFTRho[int(i * pow(numPts, 2) + j * numPts + k)], scale,
    //                            compFFTRho1[int(i * pow(numPts, 2) + j * numPts + k)]);
                    compMultFFT(compFFTRho[int(i * pow(numPts, 2) + j * numPts + k)], scale,
                                comp[axis][int(i * pow(numPts, 2) + j * numPts + k)]);
                }
            }
        }
        fftw_execute(bwrd[axis]);
        ctor(cField[axis], realField[axis]);
    }
//    diff(-1);
}

void grid::interp(vector<body>* bods){
    for (auto & body : (*bods)){
        vector<double> f(3, 0);
        vector<vector<int>> mPos = meshPos(body.pos.back());
//#pragma omp parallel for
        if (mPos.size() == 8) { /// Check inside mesh
            double xd = (body.pos.back()[0] - vecPos(mPos[0][0], mPos[0][1], mPos[0][2])[0]) / spacing;
            double yd = (body.pos.back()[1] - vecPos(mPos[0][0], mPos[0][1], mPos[0][2])[1]) / spacing;
            double zd = (body.pos.back()[2] - vecPos(mPos[0][0], mPos[0][1], mPos[0][2])[2]) / spacing;
            for (int axis = 0; axis < 3; axis++) {
                double C00 = realField[axis][int(mPos[0][0] * pow(numPts, 2) + mPos[0][1] * numPts + mPos[0][2])] *
                        (1 - xd) + xd*realField[axis][int(mPos[4][0] * pow(numPts, 2) + mPos[4][1] * numPts + mPos[4][2])];
                double C01 = realField[axis][int(mPos[1][0] * pow(numPts, 2) + mPos[1][1] * numPts + mPos[1][2])] *
                        (1 - xd) + xd*realField[axis][int(mPos[5][0] * pow(numPts, 2) + mPos[5][1] * numPts + mPos[5][2])];
                double C10 = realField[axis][int(mPos[2][0] * pow(numPts, 2) + mPos[2][1] * numPts + mPos[2][2])] *
                        (1 - xd) + xd*realField[axis][int(mPos[6][0] * pow(numPts, 2) + mPos[6][1] * numPts + mPos[6][2])];
                double C11 = realField[axis][int(mPos[3][0] * pow(numPts, 2) + mPos[3][1] * numPts + mPos[3][2])] *
                        (1 - xd) + xd*realField[axis][int(mPos[7][0] * pow(numPts, 2) + mPos[7][1] * numPts + mPos[7][2])];

                double C0 = C00 * (1-yd) + C10*yd;
                double C1 = C01 * (1-yd) + C11*yd;

                double C = C0 * (1-zd) + C1*zd;

                f[axis] = C / body.mass.back();

            }
        }
        body.acc.emplace_back(f);
//        printVec(body.acc.back());
    }
}

void grid::diff(double scale) {
//#pragma omp parallel for
    for (int i=0; i<numPts; i++){
        for (int j=0; j<numPts; j++){
            for (int k=0; k<numPts; k++){
                if (i != 0 && i != numPts-1 && j != 0 && j != numPts-1 && k != 0 && k != numPts-1) {
                    realField[0][int(i * pow(numPts, 2) + j * numPts + k)] =
                            scale*(realPot[int((i - 1) * pow(numPts, 2) + j * numPts + k)][0] -
                                    realPot[int((i + 1) * pow(numPts, 2) + j * numPts + k)][0]) / (2 * spacing);
                    realField[1][int(i * pow(numPts, 2) + j * numPts + k)] =
                            scale*(realPot[int(i * pow(numPts, 2) + (j - 1) * numPts + k)][0] -
                                    realPot[int(i * pow(numPts, 2) + (j + 1) * numPts + k)][0]) / (2 * spacing);
                    realField[2][int(i * pow(numPts, 2) + j * numPts + k)] =
                            scale*(realPot[int(i * pow(numPts, 2) + j * numPts + k - 1)][0] -
                                    realPot[int(i * pow(numPts, 2) + j * numPts + k + 1)][0]) / (2 * spacing);
                } else {
                    for (auto & axis : realField){
                        axis[int(i * pow(numPts, 2) + j * numPts + k)] = 0;
                    }
                }
            }
        }
    }
}

vector<double> grid::vecPos(int i, int j, int k){
    return {i*spacing - dim/2, j*spacing - dim/2, k*spacing - dim/2};
}

vector<vector<int>> grid::meshPos(vector<double> pos) {
//    pos = {dim/2, dim/2, dim/2};
    vector<vector<int>> meshPoints;
    vector<int> pt(3, 0);
    for (int i = 0; i < 3; i++) {
        pt[i] = floor((pos[i] + dim/2) / spacing);
//        cout <<  "compute: " << (pos[i] + dim) << endl;
    }
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                if (pt[0]+i<numPts && pt[1]+j<numPts && pt[2]+k<numPts && pt[0]+i>=0 && pt[1]+j>=0 && pt[2]+k>=0)
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
            out *= abs((spacing - dist)/spacing);
        } else
            return 0;
    }
//    cout << "ok" << endl;

    return out;
}

void grid::ctor(fftw_complex* arr, double* out){
    for (int i=0; i<numPts; i++) {
        for (int j = 0; j < numPts; j++) {
            for (int k = 0; k < numPts; k++) {
                out[int(i * pow(numPts, 2) + j * numPts + k)] =
                        arr[int(i * pow(numPts, 2) + j * numPts + k)][0];
            }
        }
    }
}

void grid::magF(){
    for (int i=0; i<numPts; i++) {
        for (int j = 0; j < numPts; j++) {
            for (int k = 0; k < numPts; k++) {
//                cout << realField1[int(i * pow(numPts, 2) + j * numPts + k)] << endl;
                realField1[int(i * pow(numPts, 2) + j * numPts + k)] =
                        pow(realField[0][int(i * pow(numPts, 2) + j * numPts + k)], 2) +
                        pow(realField[1][int(i * pow(numPts, 2) + j * numPts + k)], 2) +
                        pow(realField[2][int(i * pow(numPts, 2) + j * numPts + k)], 2);
            }
        }
    }
}

vector<vector<vector<double>>> grid::getF(int indx){
    vector<vector<vector<double>>> out;
    for (int i=0; i<numPts; i++) {
        vector<vector<double>> tmp;
        for (int j = 0; j < numPts; j++) {
            vector<double> tmp1;
            for (int k = 0; k < numPts; k++) {
                tmp1.emplace_back(realField[indx][int(i * pow(numPts, 2) + j * numPts + k)] + 1e-9);
            }
            tmp.emplace_back(tmp1);
        }
        out.emplace_back(tmp);
    }
    return out;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void compMultFFT(fftw_complex v1, fftw_complex v2, fftw_complex out) {
    double real = v1[0]*v2[0] - v1[1]*v2[1];
    double imag = v1[1]*v2[0] + v1[0]*v2[1];
    out[0] = real,
    out[1] = imag;
//    if (real){
//        cout << "(" << v1[0] << " + i" << v1[1] << ")(" <<  v2[0] << " + i" << v2[1] << ") = " << out[0] << " + i" << out[1] << endl;
//        cout << "v1[0]*v2[0] = " << v1[0]*v2[0] << endl;}

}



