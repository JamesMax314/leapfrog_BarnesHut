#ifndef VECMATHS
#define VECMATHS

#include <vector>

using namespace std;

double m_modulus(vector<double>, bool);
double dot(vector<double>&, vector<double>&, bool);
vector<double> vecAdd(vector<double>, vector<double>);
vector<int> vecAdd(vector<int>, vector<int>);
vector<double> scalMult(double, vector<double>);
vector<int> scalMult(int, vector<int>);
vector<double> compMult(vector<double>, vector<double>);
void printVec(vector<double> vec);

#endif //VECMATHS