#include <map>
#include <utility>
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

sub_grid::sub_grid(const grid& g, const sg_seed& seed, vector<int> pts, int timeStep) : grid(g.spacing, std::move(pts)){
    dt = timeStep;
    min = seed.min;
    max = seed.max;
    mainPoints = seed.indices;
}

sub_grid::sub_grid() = default;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

tree_PM::tree_PM(vector<body>& bods, double gridSpacing, double dimension, double density, double timeStep) :
        g(grid(gridSpacing, dimension)), cg(comp_grid(g)) {
    dt = timeStep;
    dim = {dimension, dimension, dimension};
    gridSpace = gridSpacing;
    den = density;
    bodies = &bods;

}

void tree_PM::genSeeds() {
    vector<int> keys;
    vector<vector<int>> indices;
    int count = 0;
    for (int i = 0; i < g.numPts[0]; i++) {
        for (int j = 0; j < g.numPts[1]; j++) {
            for (int k = 0; k < g.numPts[2]; k++) {
                g.keys[i*g.numPts[2]*g.numPts[1] + j*g.numPts[2] + k] = -1; /// Set Key value to -1
                if (g.realPot[int(i*g.numPts[2]*g.numPts[1] + j*g.numPts[2] + k)][0] > den){
                    /// Above density threshold
                    vector<int> index = {i, j, k}; /// Record index in main grid
                    indices.emplace_back(index);
                    keys.emplace_back(count); /// Generate new key
                    count ++;
                }
            }
        }
    }
    for (int i=0; i<keys.size(); i++){
        for (int j=0; j<keys.size(); j++) {
            bool adjacent = true; /// Assume adjacent
            for (int axis=0; axis<3; axis++){
                if (abs(indices[i][axis] - indices[j][axis]) > 1)
                    adjacent = false; /// if not adjacent in one dimension
            }
            if (adjacent){
                /// Set key values equal
                keys[j] = keys[i];
                g.keys[indices[j][0]*g.numPts[2]*g.numPts[1] + indices[j][1]*g.numPts[2] + indices[j][2]] = keys[i];
            }
        }
    }
    vector<sg_seed> seeds;
    /// Add indices (on main grid) to seeds based on key
    for (int i=0; i<keys.size(); i++){
        int seedIndx = -1;
        for (int s=0; s<seeds.size(); s++){
             if (seeds[s].key == keys[i])
                 seedIndx = s;
        }
        if (seedIndx < 0) {
            seeds.emplace_back(sg_seed(keys[i], indices[i]));
            seeds.back().min = indices[i];
            seeds.back().max = indices[i];
            seeds.back().den = g.realPot[int(indices[i][0]*g.numPts[2]*g.numPts[1] +
                                         indices[i][1]*g.numPts[2] +
                                         indices[i][2])][0];
        } else {
            seeds[seedIndx].indices.emplace_back(indices[i]);
            seeds[seedIndx].maxMin(indices[i]);
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
    map<int, sub_grid> sgs;
    for (auto & seed : seedVec) {
        vector<int> points(3, 0);
        for (int axis = 0; axis < 3; axis++) {
            /// Compute tha buffering required for non periodicity
            points[axis] = 3 * (seed.max[axis] - seed.min[axis] + 1);
        }
        int subTimeStep = dt / (den - seed.den + 1);
        sgs.insert(pair<int, sub_grid>(seed.key, sub_grid(g, seed, points, subTimeStep)));
    }
    sgVec = sgs;
}

void tree_PM::classiftBods() {
    for (int i=0; i<(*bodies).size(); i++){
        vector<vector<int>> pos = g.meshPos((*bodies)[i].pos.back());
        int currKey = -1;
        for (auto & po : pos){
            int k = g.keys[po[1]*g.numPts[2]*g.numPts[1] + po[1]*g.numPts[2] + po[1]];
            if (k == -1 || (k != currKey && currKey != -1))
                break;
            else
                currKey = k; /// If body is in a sub grid according to given mes pos
        }
        if (currKey >= 0)
            sgVec[currKey].activeBods.emplace_back(i); /// Adds body index to active bodies list
        else
            g.activeBods.emplace_back(i);
    }
}

void tree_PM::runTrees() {
    g.solveField();
    /// For efficiency this should be changed to the centre and size of the sub grids
    vector<double> width = g.dim;
    vector<double> centre = {0, 0, 0};
    for (auto & subG : sgVec) {
        barnesHut bh = barnesHut(*bodies, width, centre);
        for(int j=0; j<int(dt / subG.second.dt); j++) {
            subG.second.updateGrid(bodies);
            subG.second.solveField();
            cg.updateCompGrid(subG.second);
            treeMake(bh);
            interaction(bh);
            cg.interpW(bodies, false);
            bodiesUpdate(bodies, subG.second.activeBods, subG.second.dt);
            treeBreak(bh);
        }
        PBC(bodies, subG.second.activeBods, g.dim); /// Put particles into the correct location with PBCs
    }
    /// Update outside particles
    g.updateGrid(bodies);
    g.solveField();
    g.interpW(bodies, true);
    bodiesUpdate(bodies, g.activeBods, dt, g.dim);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

comp_grid::comp_grid(const grid &g) : grid(g), mainG(g) {}

void comp_grid::updateCompGrid(sub_grid & sg) {
    activeBods = sg.activeBods;
    for (int i = 0; i < numPts[0]; i++) {
        for (int j = 0; j < numPts[1]; j++) {
            for (int k = 0; k < numPts[2]; k++) {
                vector<int> sI = sg.getSubIndx({i, j, k});
                for (int axis=0; axis<3; axis++) {
                    if (sI[axis] > 0) {
                        realField[axis][int(i * numPts[2] * numPts[1] + j * numPts[2] + k)] =
                                mainG.realField[axis][int(i * numPts[2] * numPts[1] + j * numPts[2] + k)] -
                                sg.realField[axis][int(sI[0] * numPts[2] * numPts[1] + sI[1] * numPts[2] + sI[2])];
                    } else {
                        realField[axis][int(i * numPts[2] * numPts[1] + j * numPts[2] + k)] =
                                mainG.realField[axis][int(i * numPts[2] * numPts[1] + j * numPts[2] + k)];
                    }
                }
            }
        }
    }
}

