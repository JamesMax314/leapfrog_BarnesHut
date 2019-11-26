#include "bodies.h"

/// Body constructors
body::body() = default;

body::body(double &m, vector<double> &p, vector<double> &v,
           vector<double> &a, double &k, double &pe){
    mass = m;
    pos = p;
    vel = v;
    acc = a;
    ek = k;
    ep = pe;
}