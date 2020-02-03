#ifndef TPM
#define TPM

#include <fftw3.h>
#include <cmath>
#include "bodies.h"
#include "poisson.h"

class sg_seed{
public:
    int key;
    vector<vector<int>> indices;
    vector<int> min;
    vector<int> max;
    void maxMin(vector<int> & vec);
    sg_seed(int key, vector<int> & index);
};

class sub_grid : public grid{
public:
//    vector<int> padding;
    vector<int> min;
    vector<int> max;
    vector<int> activeBods;
//    vector<vector<int>> subPoints;
    vector<vector<int>> mainPoints;
    vector<int> getSubIndx(vector<int> index);
    sub_grid(const grid& g, const sg_seed& seed, vector<int> pts);
};



#endif TPM