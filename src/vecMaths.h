#ifndef VECMATHS
#define VECMATHS

#include <vector>

using namespace std;

double modulus(vector<double>, bool);
double dot(vector<double>&, vector<double>&, bool);
vector<double> vecAdd(vector<double>, vector<double>);
vector<double> scalMult(double, vector<double>&);

#endif //VECMATHS