#include <vector>
#include <algorithm>
#include <iostream>
#include "bodies.h"
#include "trees.h"
#include "vecMaths.h"
#include "leapfrog.h"

/// Leapfrog functions
void treeMake(barnesHut &hut){
    hut.treeBuild();
}

void treeBreak(barnesHut &hut){
//    cout << "breaking..." << endl;
	hut.treeChop(hut.root);
}

void interaction(barnesHut &hut){
//    cout << "interacting..." << endl;
    hut.acceleration(hut.root);
}

void bodiesUpdate(barnesHut &hut, double dt){
//    cout << "updating..." << endl;
    node* tree = hut.root;
    vector<body>* bodies = hut.bodies;
    partialUpdate(tree, bodies, dt);
}

void partialUpdate(node* tree, vector<body>* bodies, double dt){
//    cout << "num Children: " << tree->liveChildren.size() << endl;
    if (tree->liveChildren.size() != 0){
        for(auto i: tree->liveChildren){
//            cout << "i: " << i << endl;
            partialUpdate(tree->children[i], bodies, dt);
        }
    } else{
        //cout << "Updating leaf node" << endl;
        auto i = tree->bodyindx;
//        cout << "Body num: " << i << endl;
//        cout << "vel length: " << (*bodies)[i].vel.size() << endl;
//        cout << "pos length: " << (*bodies)[i].pos.size() << endl;
        (*bodies)[i].vel.emplace_back(vecAdd((*bodies)[i].vel.back(), scalMult(dt, (*bodies)[i].acc.back())));
//        cout << "new vx: " << (*bodies)[i].vel[0][0] << endl;
        (*bodies)[i].pos.emplace_back(vecAdd((*bodies)[i].pos.back(), scalMult(dt, (*bodies)[i].vel.back())));
//        cout << "new x: " << (*bodies)[i].pos[0][0] << endl;
    }
}