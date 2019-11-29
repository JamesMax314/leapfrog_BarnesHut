#include <iostream>
#include <utility>
#include "bodies.h"

using namespace std;

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

void body::setPos(vector<double> p){
    pos = std::move(p);
}
void body::setAcc(vector<double> a){
    acc = std::move(a);
}
void body::setVel(vector<double> v){
    acc = std::move(v);
}
void body::setMass(double m){
    mass = m;
}
vector<double> body::getPos(){
    return pos;
}
vector<double> body::getAcc(){
    return acc;
}
vector<double> body::getVel(){
    return vel;
}
double body::getMass(){
    return mass;
}
