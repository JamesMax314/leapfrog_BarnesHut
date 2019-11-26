#ifndef BODIES
#define BODIES

#include <vector>

using namespace std;

struct body{
	double mass{};
	vector<double> pos;
	vector<double> vel;
	vector<double> acc;
	double ek{};
	double ep{};

    body();
    body(double&, vector<double>&, vector<double>&, vector<double>&, double&, double&);
};

// Macros to retrieve body data; x is a pointer
#define Vel(x) (((body*) (x))->vel) // cast x to type body* then retrieve vel
#define Acc(x) (((body*) (x))->acc)
#define Ek(x) (((body*) (x))->ek)
#define Ep(x) (((body*) (x))->ep)

#endif //BODIES