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
    vector<double> time = {0};
    vector<double> ek;
	vector<double> ep;
	vector<bool> active;

	/// Meta Data: store in the first body for easy access in python
    vector<int> numTrees;
    vector<double> avrgR;
    vector<double> avrgM;

	void setPos(const vector<double>&);
    void setAcc(const vector<double>&);
    void setVel(const vector<double>&);
    void setMass(double);
    void setSoftening(double);
    vector<vector<double>> getPos();
    vector<vector<double>> getAcc();
    vector<vector<double>> getVel();
    vector<double> getMass();
    vector<int> get_numTrees();
    vector<double> get_avrgR();
    vector<double> get_avrgM();
    vector<int> set_numTrees();
    vector<double> set_avrgR();
    vector<double> set_avrgM();
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