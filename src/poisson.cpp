#include <fftw3.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <omp.h>
#include "vecMaths.h"
#include "poisson.h"

grid::grid(double gridSpacing, double dim): dim({dim, dim, dim}), spacing(gridSpacing) {
    for (int axis=0; axis<3; axis++)
        numPts.emplace_back((int) (dim / gridSpacing));
//    cout << "numPts: " << numPts << endl;
    realPot = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * numPts[0]*numPts[1]*numPts[2]);
    realField1 = (double *) fftw_malloc(sizeof(double) * numPts[0]*numPts[1]*numPts[2]);
    compFFTRho = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numPts[0]*numPts[1]*numPts[2]);
    compFFTRho1 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numPts[0]*numPts[1]*numPts[2]);
    keys = new int[numPts[0]*numPts[1]*numPts[2]];
    for (auto & i : realField) {
        i = new double[int(numPts[0]*numPts[1]*numPts[2])];
    }
    for (int i=0; i<3; i++) {
        comp[i] = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numPts[0]*numPts[1]*numPts[2]);
        cField[i] = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numPts[0]*numPts[1]*numPts[2]);
//        realField[i] = (double *) fftw_malloc(sizeof(double) * pow(numPts, 3));
    }
    fwrd = fftw_plan_dft_3d(numPts[0], numPts[1], numPts[2], realPot, compFFTRho, FFTW_FORWARD, FFTW_MEASURE);
//    bwrd = fftw_plan_dft_3d(numPts, numPts, numPts, compFFTRho1, realPot, FFTW_BACKWARD, FFTW_MEASURE);
    for (int i = 0; i < 3; i++) {
        bwrd[i] = fftw_plan_dft_3d(numPts[0], numPts[1], numPts[2], comp[i], cField[i], FFTW_BACKWARD, FFTW_MEASURE);
    }
}

grid::grid(double gridSpacing, vector<int> numPs): spacing(gridSpacing), numPts(numPs) {
    for (int axis=0; axis<3; axis++)
        dim.emplace_back((int) (numPs[axis] * gridSpacing));
    realPot = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * numPts[0]*numPts[1]*numPts[2]);
    realField1 = (double *) fftw_malloc(sizeof(double) * numPts[0]*numPts[1]*numPts[2]);
    compFFTRho = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numPts[0]*numPts[1]*numPts[2]);
    compFFTRho1 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numPts[0]*numPts[1]*numPts[2]);
    keys = new int[numPts[0]*numPts[1]*numPts[2]];
    for (auto & i : realField) {
        i = new double[int(numPts[0]*numPts[1]*numPts[2])];
    }
    for (int i=0; i<3; i++) {
        comp[i] = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numPts[0]*numPts[1]*numPts[2]);
        cField[i] = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * numPts[0]*numPts[1]*numPts[2]);
    }
    fwrd = fftw_plan_dft_3d(numPts[0], numPts[1], numPts[2], realPot, compFFTRho, FFTW_FORWARD, FFTW_MEASURE);
    for (int i = 0; i < 3; i++) {
        bwrd[i] = fftw_plan_dft_3d(numPts[0], numPts[1], numPts[2], comp[i], cField[i], FFTW_BACKWARD, FFTW_MEASURE);
    }
}

grid::grid(const grid & g) : numPts(g.numPts) {
    for (int axis=0; axis<3; axis++)
        dim.emplace_back(g.dim[axis]);
    for (auto & i : realField) {
        i = new double[int(numPts[0]*numPts[1]*numPts[2])];
    }
}

void grid::updateGrid(vector<body>* bods){
    /// Row major format
//    realPot[int(i*pow(numPts, 2) + j*numPts + k)] = 0;
//#pragma omp parallel for //default(none) shared(bods)
    for (int i = 0; i < numPts[0]; i++) {
        for (int j = 0; j < numPts[1]; j++) {
            for (int k = 0; k < numPts[2]; k++) {
                realPot[int(i*numPts[2]*numPts[1] + j*numPts[2] + k)][0] = 0;
            }
        }
    }

    for (auto & bIndx : activeBods){
        vector<vector<int>> mPos = meshPos((*bods)[bIndx].pos.back());
        for (auto vec : mPos) {
            realPot[int(vec[0]*numPts[2]*numPts[1] + vec[1]*numPts[2] + vec[2])][0] +=
                   w({vec[0], vec[1], vec[2]}, (*bods)[bIndx]) * (*bods)[bIndx].mass.back() / spacing;
        }
    }
}

void grid::solveField(){
    fftw_execute(fwrd);
//#pragma omp parallel for // default(none) shared(axis)
    for (int axis=0; axis<3; axis++) {
        for (int i = 0; i < numPts[0]; i++) { //floor(numPts/2) + 1
            for (int j = 0; j < numPts[1]; j++) {
                for (int k = 0; k < numPts[2]; k++) {
                    /// Accounting for weird fft structure
                    int kx = (i <= numPts[0] / 2) ? i : i - numPts[0];
                    int ky = (j <= numPts[1] / 2) ? j : j - numPts[1];
                    int kz = (k <= numPts[2] / 2) ? k : k - numPts[2];
                    vector<int> ks = {kx, ky, kz};

                    fftw_complex scale;

                    /// Normalize FFT with 1/N^3 the inde must be scaled to give k;
                    /// The k vectors are given by index /(spacing * N);
                    /// The Greens function for force is -ik-> / k^2;
                    /// The overall scaling is given by: spacing * N / N^3 = spacing / N^2.
                    scale[1] = (i == 0 && j == 0 && k == 0) ? 0 : spacing * ks[axis] / pow(numPts[0], 2) *
                            4 * pi * G /
                            (numPts[0]*numPts[1]*numPts[2] *
                            (kx * kx / pow(numPts[0], 2) +
                            ky * ky / pow(numPts[1], 2) +
                            kz * kz / pow(numPts[2], 2))); //dim*4*pi*G
                    scale[0] = 0;

    //                compMultFFT(compFFTRho[int(i * pow(numPts, 2) + j * numPts + k)], scale,
    //                            compFFTRho1[int(i * pow(numPts, 2) + j * numPts + k)]);
                    compMultFFT(compFFTRho[int(i * numPts[2]*numPts[1] + j * numPts[2] + k)], scale,
                                comp[axis][int(i * numPts[2]*numPts[1] + j * numPts[2] + k)]);
                }
            }
        }
        fftw_execute(bwrd[axis]);
        ctor(cField[axis], realField[axis]);
    }
//    diff(-1);
}

void grid::interpW(vector<body>* bods, bool resetForce){
    for (auto & bIndx : activeBods){
        vector<double> f(3, 0);
        vector<vector<int>> mPos = meshPos((*bods)[bIndx].pos.back());
//#pragma omp parallel for
        for (auto vec : mPos) {
            for (int axis = 0; axis < 3; axis++) {
                f[axis] += w({vec[0], vec[1], vec[2]}, (*bods)[bIndx]) *
                        realField[axis][int(vec[0] * numPts[2]*numPts[1] + vec[1] * numPts[2] + vec[2])] /
                        (*bods)[bIndx].mass.back();
            }
        }
        if (resetForce)
            (*bods)[bIndx].acc.emplace_back(f);
        else
            (*bods)[bIndx].acc.back() = vecAdd((*bods)[bIndx].acc.back(), f);
//        printVec(body.acc.back());
    }
}

/// !!! Conatains legacy code; no actvie bodies list !!!
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
                double C00 = realField[axis][int(mPos[0][0] * numPts[2]*numPts[1] + mPos[0][1] * numPts[2] + mPos[0][2])] *
                        (1 - xd) + xd*realField[axis][int(mPos[4][0] * numPts[2]*numPts[1] + mPos[4][1] * numPts[2] + mPos[4][2])];
                double C01 = realField[axis][int(mPos[1][0] * numPts[2]*numPts[1] + mPos[1][1] * numPts[2] + mPos[1][2])] *
                        (1 - xd) + xd*realField[axis][int(mPos[5][0] * numPts[2]*numPts[1] + mPos[5][1] * numPts[2] + mPos[5][2])];
                double C10 = realField[axis][int(mPos[2][0] * numPts[2]*numPts[1] + mPos[2][1] * numPts[2] + mPos[2][2])] *
                        (1 - xd) + xd*realField[axis][int(mPos[6][0] * numPts[2]*numPts[1] + mPos[6][1] * numPts[2] + mPos[6][2])];
                double C11 = realField[axis][int(mPos[3][0] * numPts[2]*numPts[1] + mPos[3][1] * numPts[2] + mPos[3][2])] *
                        (1 - xd) + xd*realField[axis][int(mPos[7][0] * numPts[2]*numPts[1] + mPos[7][1] * numPts[2] + mPos[7][2])];

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
    for (int i=0; i<numPts[0]; i++){
        for (int j=0; j<numPts[1]; j++){
            for (int k=0; k<numPts[2]; k++){
                if (i != 0 && i != numPts[0]-1 && j != 0 && j != numPts[1]-1 && k != 0 && k != numPts[2]-1) {
                    realField[0][int(i * numPts[2]*numPts[1] + j * numPts[2] + k)] =
                            scale*(realPot[int((i - 1) * numPts[2]*numPts[1] + j * numPts[2] + k)][0] -
                                    realPot[int((i + 1) * numPts[2]*numPts[1] + j * numPts[2] + k)][0]) / (2 * spacing);
                    realField[1][int(i * numPts[2]*numPts[1] + j * numPts[2] + k)] =
                            scale*(realPot[int(i * numPts[2]*numPts[1] + (j - 1) * numPts[2] + k)][0] -
                                    realPot[int(i * numPts[2]*numPts[1] + (j + 1) * numPts[2] + k)][0]) / (2 * spacing);
                    realField[2][int(i * numPts[2]*numPts[1] + j * numPts[2] + k)] =
                            scale*(realPot[int(i * numPts[2]*numPts[1] + j * numPts[2] + k - 1)][0] -
                                    realPot[int(i * numPts[2]*numPts[1] + j * numPts[2] + k + 1)][0]) / (2 * spacing);
                } else {
                    for (auto & axis : realField){
                        axis[int(i * numPts[2]*numPts[1] + j * numPts[2] + k)] = 0;
                    }
                }
            }
        }
    }
}

vector<double> grid::vecPos(int i, int j, int k){
    return {i*spacing - dim[0]/2, j*spacing - dim[1]/2, k*spacing - dim[2]/2};
}

vector<vector<int>> grid::meshPos(vector<double> pos) {
//    pos = {dim/2, dim/2, dim/2};
    vector<vector<int>> meshPoints;
    vector<int> pt(3, 0);
    for (int i = 0; i < 3; i++) {
        pt[i] = floor((pos[i] + dim[i]/2) / spacing);
//        cout <<  "compute: " << (pos[i] + dim) << endl;
    }
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                if (pt[0]+i<numPts[0] && pt[1]+j<numPts[1] && pt[2]+k<numPts[2] && pt[0]+i>=0 && pt[1]+j>=0 && pt[2]+k>=0)
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
    for (int i=0; i < numPts[0]; i++) {
        for (int j = 0; j < numPts[1]; j++) {
            for (int k = 0; k < numPts[2]; k++) {
                out[int(i * numPts[2]*numPts[1] + j * numPts[2] + k)] =
                        arr[int(i * numPts[2]*numPts[1] + j * numPts[2] + k)][0];
            }
        }
    }
}

void grid::magF(){
    for (int i=0; i<numPts[0]; i++) {
        for (int j = 0; j < numPts[1]; j++) {
            for (int k = 0; k < numPts[2]; k++) {
//                cout << realField1[int(i * pow(numPts, 2) + j * numPts + k)] << endl;
                realField1[int(i * numPts[2]*numPts[1] + j * numPts[2] + k)] =
                        pow(realField[0][int(i * numPts[2]*numPts[1] + j * numPts[2] + k)], 2) +
                        pow(realField[1][int(i * numPts[2]*numPts[1] + j * numPts[2] + k)], 2) +
                        pow(realField[2][int(i * numPts[2]*numPts[1] + j * numPts[2] + k)], 2);
            }
        }
    }
}

vector<vector<vector<double>>> grid::getF(int indx){
    vector<vector<vector<double>>> out;
    for (int i=0; i<numPts[0]; i++) {
        vector<vector<double>> tmp;
        for (int j = 0; j < numPts[1]; j++) {
            vector<double> tmp1;
            for (int k = 0; k < numPts[2]; k++) {
                tmp1.emplace_back(realField[indx][int(i * numPts[2]*numPts[1] + j * numPts[2] + k)] + 1e-9);
            }
            tmp.emplace_back(tmp1);
        }
        out.emplace_back(tmp);
    }
    return out;
}

void grid::initActiveBods(vector<body>* bods){
    vector<int> tmp;
    for (int i=0; i<(*bods).size(); i++){
        tmp.emplace_back(i);
    }
    activeBods = tmp;
}

grid::grid() = default;
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



