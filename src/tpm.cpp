#include <map>
#include <omp.h>
#include <utility>
#include <iostream>
#include <cmath>
#include "poisson.h"
#include "vecMaths.h"
#include "bodies.h"
#include "trees.h"
#include "tpm.h"
#include "leapfrog.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
sg_seed::sg_seed(int key, vector<int> & index) : key(key) {
    indices.emplace_back(index);
    den = 0;
}

void sg_seed::maxMin(vector<int> & vec){
    for (int axis=0; axis<3; axis++){
        if (vec[axis] < min[axis])
            min[axis] = vec[axis];
        else if (vec[axis] > max[axis])
            max[axis] = vec[axis];
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

vector<int> sub_grid::getSubIndx(vector<int> index) {
    vector<int> out(3, 0);
    for (int axis=0; axis<3; axis++){
        out[axis] = index[axis] - min[axis] + int(numPts[axis]/3);
        if (out[axis] < int(numPts[axis]/3) || out[axis] >= max[axis] - min[axis] + int(numPts[axis]/3))
            out = {-1, -1, -1};
    }
    return out;
}

sub_grid::sub_grid(const grid& g, const sg_seed& seed, vector<int> pts, double timeStep) : grid(g.spacing, move(pts)){
    dt = timeStep;
    min = seed.min;
    max = seed.max;
    mainPoints = seed.indices;
    activeBods = {};
}

sub_grid::~sub_grid() {
}

sub_grid::sub_grid(const sub_grid& sg): grid(sg) {
    dt = sg.dt;
    min = sg.min;
    max = sg.max;
    activeBods = sg.activeBods;
    mainPoints = sg.mainPoints;
}

sub_grid::sub_grid() = default;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

tree_PM::tree_PM(vector<body>& bods, double gridSpacing, double dimension, double density, double timeStep, double time) :
        g(grid(gridSpacing, dimension)), cg(comp_grid(g)) {
    dt = timeStep;
    dim = {dimension, dimension, dimension};
    gridSpace = gridSpacing;
    den = density;
    bodies = &bods;
    nu = 0;
    c = 1.066e-12;
    t = time;
}

tree_PM::tree_PM(vector<body>& bods, double gridSpacing, double dimension, double density, double timeStep) :
        g(grid(gridSpacing, dimension)), cg(comp_grid(g)) {
    dt = timeStep;
    dim = {dimension, dimension, dimension};
    gridSpace = gridSpacing;
    den = density;
    bodies = &bods;
    nu = 0;
    c = 1.066e-12;
    t = 1e9*365*24*3600;
}

tree_PM::tree_PM() : g(), cg() {
    nu = 0;
    c = 1.066e-12;
}

void tree_PM::genSeeds() {
    /// Compute densities
    g.updateGrid(bodies);
    /// keys are used to determine whether points are adjacent
    vector<int> keys;
    /// Stores the indices of grid points that surpass the density threshold
    vector<vector<int>> indices;
    int count = 0;
//#pragma omp parallel for
    /// Iterate through each grid point and check the density
    for (int i = 0; i < g.numPts[0]; i++) {
        for (int j = 0; j < g.numPts[1]; j++) {
            for (int k = 0; k < g.numPts[2]; k++) {
                /// Set Key value to -1
                g.keys[i*g.numPts[2]*g.numPts[1] + j*g.numPts[2] + k] = -1;
                if (g.realPot[int(i*g.numPts[2]*g.numPts[1] + j*g.numPts[2] + k)][0] > den){
                    /// Above density threshold
//                    cout << "den: " << g.realPot[int(i*g.numPts[2]*g.numPts[1] + j*g.numPts[2] + k)][0] << endl;
                    /// Record index in main grid
                    vector<int> index = {i, j, k};
                    indices.emplace_back(index);
                    /// Generate new key
                    keys.emplace_back(count);
                    g.keys[i*g.numPts[2]*g.numPts[1] + j*g.numPts[2] + k] = count;
                    count ++;
                }
            }
        }
    }
//    cout << "indices len: " << indices.size() << endl;
    for (int i=0; i<keys.size(); i++){
        for (int j=0; j<keys.size(); j++) {
            bool adjacent = true; /// Assume adjacent
            for (int axis=0; axis<3; axis++){
                if (pow(indices[i][axis] - indices[j][axis], 2.) > 1)
                    adjacent = false; /// if not adjacent in one dimension
            }
            if (adjacent){
                /// Set key values equal
//                cout << "adjacent: " << keys[i] << keys[j] << endl;
                keys[j] = keys[i];
                g.keys[indices[j][0]*g.numPts[2]*g.numPts[1] + indices[j][1]*g.numPts[2] + indices[j][2]] = keys[i];
            }
        }
    }
//    cout << "keys: " << endl;
//    for (int i=0; i<keys.size(); i++){
//        cout << keys[i] << " " << endl;
//    }
    vector<sg_seed> seeds;
    /// Add indices (on main grid) to seeds based on key
    for (int i=0; i<keys.size(); i++){
        int seedIndx = -1;
        for (int s=0; s<seeds.size(); s++){
             if (seeds[s].key == keys[i])
                 seedIndx = s;
        }
        if (seedIndx < 0) {
            /// Initialize a new seed
            seeds.emplace_back(sg_seed(keys[i], indices[i]));
            seeds.back().min = indices[i];
            seeds.back().max = indices[i];
            seeds.back().den = g.realPot[int(indices[i][0]*g.numPts[2]*g.numPts[1] +
                                         indices[i][1]*g.numPts[2] +
                                         indices[i][2])][0];
        } else {
            /// Add current main grid index to seed
            seeds[seedIndx].indices.emplace_back(indices[i]);
            seeds[seedIndx].maxMin(indices[i]);
            /// Update the maximum density in the subgrid
            if (seeds[seedIndx].den < g.realPot[int(indices[i][0]*g.numPts[2]*g.numPts[1] +
                                                    indices[i][1]*g.numPts[2] +
                                                    indices[i][2])][0])
                seeds[seedIndx].den = g.realPot[int(indices[i][0]*g.numPts[2]*g.numPts[1] +
                                                    indices[i][1]*g.numPts[2] +
                                                    indices[i][2])][0];
        }
    }
    seedVec = seeds;
}

void tree_PM::genSubGrids(){
    /// Uses the seeds to initialize sub grids
    map<int, sub_grid> sgs;
    for (auto & seed : seedVec) {
        vector<int> points(3, 0);
        for (int axis = 0; axis < 3; axis++) {
            /// Compute the buffering required for non periodicity
            points[axis] = 3 * (seed.max[axis] - seed.min[axis] + 1);
        }
        /// Set variable time step
        double subTimeStep = dt; // / (den - seed.den + 1);
        /// Add subgrid to sgs
        sub_grid tmp = sub_grid(g, seed, points, subTimeStep);
        sgs.insert(pair<int, sub_grid>(seed.key, tmp));
//        cout << sgs[seed.key].realPot[0][0] << endl;
    }
    sgVec = sgs;
}

void tree_PM::classiftBods() {
    /// Activates bodies in appropriate sub grid
    g.activeBods.erase(g.activeBods.begin(), g.activeBods.end());
    for (int i=0; i<(*bodies).size(); i++){
        /// Get the <=8 nearest grid points
        vector<vector<int>> pos = g.meshPos((*bodies)[i].pos.back());
        int currKey = -1;
        /// For each mesh pos
        for (auto & po : pos) {
            int k = g.keys[po[1] * g.numPts[2] * g.numPts[1] + po[1] * g.numPts[2] + po[1]];
//            cout << "k: " << k << endl;
            /// If k != -1 then particle is in high density region
            if (k != -1) {
                /// If HD
                currKey = k;
                break;
            }
        }
        if (currKey >= 0) {
            sgVec[currKey].activeBods.emplace_back(i); /// Adds body index to active bodies list in subgrid
//            cout << "body: " << i << endl;
        }
        else
            g.activeBods.emplace_back(i); /// Adds body to main grid
    }
}

void tree_PM::runTrees() {
    g.solveField();
    update_a_t();
    /// For efficiency this should be changed to the centre and size of the sub grids
    vector<double> width = g.dim;
    vector<double> centre = {0, 0, 0};
//#pragma omp parallel for
//    cout << "num subGs " << sgVec.size() << endl;
    /// For each subgrid
    double partAvrgM = 0;
    double partAvrgR = 0;
    for (auto &subG : sgVec) {
        vector<int> tmp = vecAdd(subG.second.max, scalMult(-1, subG.second.min));
        partAvrgR += 1 + m_modulus({double(tmp[0]), double(tmp[1]), double(tmp[2])}, false);
        if (subG.second.activeBods.size() > 0) {
//            cout << "Running tree" << endl;
            partAvrgM += double(subG.second.activeBods.size());
//            cout << subG.second.activeBods.size() << endl;
            /// Initialize tree section
            barnesHut bh = barnesHut(*bodies, width, centre);
            bh.activeBods = subG.second.activeBods;
//        cout << subG.second.realPot[0][0] << endl;
//        cout << "sub iters: " << int(dt / subG.second.dt) << endl;
            /// For each time step
            for (int j = 0; j < int(dt / subG.second.dt); j++) {
                /// Update potentials in sub grid then solve field
                subG.second.updateGrid(bodies);
                subG.second.solveField();
                /// Remove the subgrid forces form the main grid and assign to comp grid
                cg.updateCompGrid(subG.second);
                /// Build tree and compute forces on particles
                treeMake(bh);
                interaction(bh);
                /// Add forces from the grid points
                cg.interpW(bodies, false);
                /// Update the pos and vel
                bodiesUpdate(bodies, subG.second.activeBods, t, subG.second.dt);
                treeBreak(bh);
            }
            /// Put particles into the correct location with PBCs
            /// i.e. Account for the particle moving outside the boundary box
            PBC(bodies, subG.second.activeBods, g.dim);
        }
    }
    /// Store the meta date
    if (sgVec.size() == 0) {
        avrgM.emplace_back(0.);
        avrgR.emplace_back(0.);
    } else {
        avrgM.emplace_back(partAvrgM / sgVec.size());
        avrgR.emplace_back(partAvrgR / sgVec.size());
    }
    numTrees.emplace_back(sgVec.size());
    /// Update outside particles
    g.solveField();
    g.interpW(bodies, true);
    bodiesUpdate(bodies, g.activeBods, t, dt, g.dim);
}

void tree_PM::update_a_t() {
    /// On the first run through the average comoving density is computed
    if (nu == 0){
        for (int i=0; i<g.numPts[0]*g.numPts[1]*g.numPts[2]; i++){
            nu += g.realPot[i][0];
        }
        nu = nu / g.numPts[0]*g.numPts[1]*g.numPts[2];
    }
    /// The scale factor is calculated according to the Freidmann equation
//    a =  // pow(8*g.pi*g.G*nu/3, 1/3) * pow(t*3/2, 2/3);
//    ad = c*2/3*pow(t, -1./3.);
//    add = -c*2/3*1/3*pow(t, -4./3.);
//    cout << "dv: " << dt*add/a*dim[0] << ", V: " << ad/a*dim[0] << endl;
//    for (int i=0; i<3; i++) { dim[i] += ad*dt/a; }
    t += dt;
}

double tree_PM::get_a(double tc) {
    return c*pow(tc, 2./3.);
}

double tree_PM::get_ad(double tc) {
    return c*2/3*pow(tc, -1./3.);
}

double tree_PM::get_add(double tc) {
    return c*2/3*-1/3*pow(tc, -4./3.);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

comp_grid::comp_grid(const grid &g) : grid(g, true), mainG(g) {}
comp_grid::comp_grid() {}

void comp_grid::updateCompGrid(sub_grid & sg) {
    /// Remove regional forces from the main grid to generate a comp grid
    activeBods = sg.activeBods;
    for (int i = 0; i < numPts[0]; i++) {
        for (int j = 0; j < numPts[1]; j++) {
            for (int k = 0; k < numPts[2]; k++) {
                vector<int> sI = sg.getSubIndx({i, j, k});
                for (int axis=0; axis<3; axis++) {
                    if (sI[axis] >= 0) {
                        realField[axis][int(i * numPts[2] * numPts[1] + j * numPts[2] + k)] =
                              mainG.realField[axis][int(i * numPts[2] * numPts[1] + j * numPts[2] + k)] -
                                sg.realField[axis][int(sI[0] * sg.numPts[2] * sg.numPts[1] + sI[1] * sg.numPts[2] + sI[2])];
                    } else {
                        realField[axis][int(i * numPts[2] * numPts[1] + j * numPts[2] + k)] =
                                mainG.realField[axis][int(i * numPts[2] * numPts[1] + j * numPts[2] + k)];
                    }
                }
            }
        }
    }
}


