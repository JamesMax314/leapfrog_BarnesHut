#include <iostream>
#include <utility>
#include "bodies.h"

using namespace std;

/// Body constructors
body::body() = default;

body::body(double &m, vector<double> &p, vector<double> &v, vector<double> &a){
    mass.emplace_back(m);
    pos.emplace_back(p);
    vel.emplace_back(v);
    acc.emplace_back(a);
    active.emplace_back(true);
}

void body::setPos(const vector<double>& p){
    pos.emplace_back(p);
}
void body::setAcc(const vector<double>& a){
    acc.emplace_back(a);
}
void body::setVel(const vector<double>& v){
    vel.emplace_back(v);
}
void body::setMass(double m){
    mass.emplace_back(m);
}
vector<vector<double>> body::getPos(){
    return pos;
}
vector<vector<double>> body::getAcc(){
    return acc;
}
vector<vector<double>> body::getVel(){
    return vel;
}
vector<double> body::getMass(){
    return mass;
}
