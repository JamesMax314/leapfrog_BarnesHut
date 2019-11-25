#ifndef BODIES
#define BODIES

struct body{
	double mass;
	vector<double> pos;
	vector<double> vel;
	vector<double> acc;
	double ek;
	double ep;
};

body::body(){};

body::body(double m, vector<double> p, vector<double> v,
 vector<double> a, double k, double pe){
	mass = m;
	pos = p;
	vel = v;
	acc = a;
	ek = k;
	ep = pe;
};

// Macros to retrieve body data; x is a pointer
#define Vel(x) (((body*) (x))->vel) // cast x to type body* then retrieve vel
#define Acc(x) (((body*) (x))->acc)
#define Ek(x) (((body*) (x))->ek)
#define Ep(x) (((body*) (x))->ep)

#endif //BODIES