#ifndef TPM
#define TPM

#include <fftw3.h>
#include <cmath>
#include <map>
#include "bodies.h"
#include "poisson.h"

using namespace std;

class sg_seed{
public:
    int key;
    double den;
    vector<vector<int>> indices;
    vector<int> min;
    vector<int> max;
    void maxMin(vector<int> & vec);
    sg_seed(int key, vector<int> & index);
};

class sub_grid : public grid{
public:
//    vector<int> padding;
    double dt{};
    vector<int> min;
    vector<int> max;
//    vector<int> activeBods;
    vector<vector<int>> mainPoints;

    //    vector<vector<int>> subPoints;
    vector<int> getSubIndx(vector<int> index);
    sub_grid(const grid& g, const sg_seed& seed, vector<int> pts, double timeStep);
    sub_grid();
    sub_grid(const sub_grid &sg);
    ~sub_grid();
};

class comp_grid : public grid{
public:
    grid mainG;
    void updateCompGrid(sub_grid & sg);
    explicit comp_grid(const grid& g);
};

class tree_PM{
public:
    double dt;
    double den;
    double gridSpace;
    vector<double> dim;
    vector<sg_seed> seedVec;
    map<int, sub_grid> sgVec;
    vector<body>* bodies;

    grid g;
    comp_grid cg;

    void genSeeds();
    void genSubGrids();
    void classiftBods();
    void runTrees();

    tree_PM(vector<body>& bods, double gridSpacing, double dim, double density, double timeStep);

};



#endif //TPM