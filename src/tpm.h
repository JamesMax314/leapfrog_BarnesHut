#ifndef TPM
#define TPM

#include <fftw3.h>
#include <cmath>
#include "bodies.h"
#include "poisson.h"

class sub_grid : public grid{
public:
    vector<int> padding;
    vector<vector<int>> subPoints;
    vector<vector<int>> mainPoints;
    sub_grid(grid g, int pts);
};

class sg_seed{
public:
    int key;
    vector<vector<int>> indices;
    vector<int> min;
    vector<int> max;
    void maxMin(vector<int> & vec);
    sg_seed(int key, vector<int> & index);
};

#endif TPM