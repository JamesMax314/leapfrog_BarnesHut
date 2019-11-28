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
	hut.treeChop(hut.root);
}

void interaction(barnesHut &hut){
    hut.acceleration(hut.root);
}

void bodiesUpdate(barnesHut &hut, double dt){
    node* tree = hut.root;
    vector<body>* bodies = hut.bodies;
    partialUpdate(tree, bodies, dt);
}

void partialUpdate(node* tree, vector<body>* bodies, double dt){
    if (tree->numChildren != 0){
        for(auto i: tree->liveChildren){
            partialUpdate(tree->children[i], bodies, dt);
        }
    } else{
        //cout << "Updating leaf node" << endl;
        auto i = tree->bodyindx;
        (*bodies)[i].vel = vecAdd((*bodies)[i].vel, scalMult(dt, (*bodies)[i].acc));
        //cout << "new vx: " << (*bodies)[i].vel[0] << endl;
        (*bodies)[i].pos = vecAdd((*bodies)[i].pos, scalMult(dt, (*bodies)[i].vel));
        //cout << "new x: " << (*bodies)[i].pos[0] << endl;
    }
}