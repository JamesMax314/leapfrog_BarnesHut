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
void body::setSoftening(double s){
    softening = s;
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
double body::getSoftening(){
    return softening;
}

vector<int> body::get_numTrees() {
    return numTrees;
}

vector<double> body::get_avrgR() {
    return avrgR;
}

vector<double> body::get_avrgM() {
    return avrgM;
}

vector<int> body::set_numTrees() {
    return vector<int>();
}

vector<double> body::set_avrgR() {
    return vector<double>();
}

vector<double> body::set_avrgM() {
    return vector<double>();
}
