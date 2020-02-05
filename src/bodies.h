#ifndef BODIES
#define BODIES

#include <vector>

using namespace std;

struct body{
    double softening = 0;
	vector<double> mass;
	vector<vector<double>> pos;
	vector<vector<double>> vel;
	vector<vector<double>> acc;
    vector<double> time;
    vector<double> ek;
	vector<double> ep;
	vector<bool> active;

	void setPos(const vector<double>&);
    void setAcc(const vector<double>&);
    void setVel(const vector<double>&);
    void setMass(double);
    void setSoftening(double);
    vector<vector<double>> getPos();
    vector<vector<double>> getAcc();
    vector<vector<double>> getVel();
    vector<double> getMass();
    double getSoftening();

    body();
    body(double& mass, vector<double>& pos, vector<double>& vel, vector<double>& acc);
};

// Macros to retrieve body data; x is a pointer
#define Vel(x) (((body*) (x))->vel) // cast x to type body* then retrieve vel
#define Acc(x) (((body*) (x))->acc)
#define Ek(x) (((body*) (x))->ek)
#define Ep(x) (((body*) (x))->ep)

#endif //BODIES