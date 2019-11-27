#include <iostream>
#include "bodies.h"
#include "trees.h"
#include "leapfrog.h"

using namespace std;

int main(){
    double mass = 1e24;
    vector<double> centre = {0, 0, 0};
    vector<double> pos = {1.1e10, 0, 0};
    vector<double> pos1 = {1e10, 0, 0};
    vector<double> vel = {0, 0, 0};
    vector<double> acc = {0, 0, 0};
    double ek = 0;
    double ep = 0;
    vector<body> bodies;
    bodies.push_back(body(mass, pos, vel, acc, ek, ep));
    bodies.push_back(body(mass, pos1, vel, acc, ek, ep));
    vector<double> dim = {10e10, 10e10, 10e10};
    barnesHut bh = barnesHut(bodies, dim, centre);
    treeMake(bh);
    cout << "tree made" << endl;
    interaction(bh);
    cout << "acceleration: ";
    for(int i=0; i<3; i++){
        cout << (*bh.bodies)[1].acc[i] << ", ";
    }
    cout << endl;
    cout << "acceleration: ";
    for(int i=0; i<3; i++){
        cout << bodies[1].acc[i] << ", ";
    }
    cout << endl;
}