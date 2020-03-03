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
    comp_grid();
};

class tree_PM{
public:
    double nu;
    double c;
    double a;
    double ad;
    double add;
    double t;
    double dt{};
    double den{};
    double gridSpace{};
    vector<double> dim;
    vector<sg_seed> seedVec;
    map<int, sub_grid> sgVec;
    vector<body>* bodies{};

    grid g;
    comp_grid cg;

    void genSeeds();
    void genSubGrids();
    void classiftBods();
    void runTrees();
    void update_a_t();

    double get_a(double t);
    double get_ad(double t);
    double get_add(double t);

    tree_PM(vector<body>& bods, double gridSpacing, double dim, double density, double timeStep, double time);
    tree_PM(vector<body>& bods, double gridSpacing, double dim, double density, double timeStep);
    tree_PM();

};



#endif //TPM