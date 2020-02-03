#include "poisson.h"
#include "vecMaths.h"
#include "bodies.h"
#include "trees.h"
#include "tpm.h"


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
sg_seed::sg_seed(int key, vector<int> & index) : key(key) {
    indices.emplace_back(index);
}

void sg_seed::maxMin(vector<int> & vec){
    for (int axis=0; axis<3; axis++){
        if (vec[axis] < min[axis])
            min[axis] = vec[axis];
        else if (vec[axis] > max[axis])
            max[axis] = vec[axis];
    }
}

vector<int> sub_grid::getSubIndx(vector<int> index) {
    vector<int> out(3, 0);
    for (int axis=0; axis<3; axis++){
        out[axis] = index[axis] - min[axis] + int(numPts[axis]/3);
    }
    return out;
}

sub_grid::sub_grid(const grid& g, const sg_seed& seed, vector<int> pts) : grid(g.spacing, pts){
    min = seed.min;
    max = seed.max;
    mainPoints = seed.indices;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

vector<sg_seed> genSeeds(grid & g, double den) {
    vector<int> keys;
    vector<vector<int>> indices;
    int count = 0;
    for (int i = 0; i < g.numPts[0]; i++) {
        for (int j = 0; j < g.numPts[1]; j++) {
            for (int k = 0; k < g.numPts[2]; k++) {
                if (g.realPot[int(i*g.numPts[2]*g.numPts[1] + j*g.numPts[2] + k)][0] > den){
                    vector<int> index = {i, j, k};
                    indices.emplace_back(index);
                    keys.emplace_back(count);
                    count ++;
                }
            }
        }
    }
    for (int i=0; i<keys.size(); i++){
        for (int j=0; j<keys.size(); j++) {
            for (int axis=0; axis<3; axis++){
                if (abs(indices[i][axis] - indices[j][axis]) <= 1)
                    keys[j] = keys[i];
            }
        }
    }
    vector<sg_seed> seeds;
    for (int i=0; i<keys.size(); i++){
        int seedIndx = -1;
        for (int s=0; s<seeds.size(); s++){
             if (seeds[s].key == keys[i])
                 seedIndx = s;
        }
        if (seedIndx < 0) {
            seeds.emplace_back(sg_seed(keys[i], indices[i]));
            seeds[0].min = indices[i];
            seeds[0].max = indices[i];
        } else
            seeds[seedIndx].indices.emplace_back(indices[i]);
    }
    return seeds;
}

sub_grid makeSG(grid & g, sg_seed seed){
    vector<int> points(3, 0);
    for (int axis=0; axis<3; axis++){
        points[axis] = 3*(seed.max[axis] - seed.min[axis]);
    }
    sub_grid sg = sub_grid(g, seed, points);
}

void classifyBods(vector<body> &bods, vector<sub_grid> & subGs){
    for (int i=0; i<bods.size(); i++){
        if (){

        }
    }
}