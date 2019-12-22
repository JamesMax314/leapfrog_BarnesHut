#ifndef PYINTERFACE
#define PYINTERFACE

#include "trees.h"

class progbar{
public:
    double max;
    double length;

    progbar(double max, double barLength);
    void update(int current);
};

vector<body> basicRun(vector<body>&, vector<double> centre, vector<double> dim, int numIter, double dt);

#endif //PYINTERFACE