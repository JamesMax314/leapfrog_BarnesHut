#include <iostream>
#include <fstream>
#include "save.h"
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
    vector<vector<body>> vecBodies;
    bodies.emplace_back(mass, pos, vel, acc, ek, ep);
    bodies.emplace_back(mass, pos1, vel, acc, ek, ep);
    vector<double> dim = {10e10, 10e10, 10e10};
    barnesHut bh = barnesHut(bodies, dim, centre);

    for(int j=0; j<10000; j++) {
        //cout << "making tree..." << endl;
        treeMake(bh);
        //cout << "tree made" << endl;
        interaction(bh);
//        cout << "acceleration: ";
//        for (int i = 0; i < 3; i++) {
//            cout << (*bh.bodies)[1].acc[i] << ", ";
//        }
//        cout << endl;
//        cout << "acceleration: ";
//        for (int i = 0; i < 3; i++) {
//            cout << bodies[1].acc[i] << ", ";
//        }
//        cout << endl;
//
//        cout << "updating bodies..." << endl;
        bodiesUpdate(bh, 1e3);
//        for (auto b : bodies) {
//            for (int i = 0; i < 3; i++) {
//                cout << b.pos[i] << ", ";
//            }
//            cout << endl;
//        }

//        cout << "breaking tree..." << endl;
        treeBreak(bh);
//        for (int i = 0; i < 8; i++) {
//            cout << bh.root->num << endl;
//        }
        vecBodies.emplace_back(bodies);
    }
    cout << "done" << endl;
    cout << "saving data..." << endl;
    ofstream file("data.txt");
}