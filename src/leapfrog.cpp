#include <vector>
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

void bodiesUpdate(barnesHut &hut, double dt){
    node* tree = hut.root;
    vector<body> bodies = hut.bodies;
    partialUpdate(tree, bodies, dt);
}

void partialUpdate(node* tree, vector<body> &bodies, double dt){
    if (tree->children){
        node* child = tree->children;
        for(int i=0; i<tree->numChildren; i++){
            partialUpdate(child, bodies, dt);
            child++;
        }
    } else{
        body b = bodies[tree->bodyindx];
        b.vel = vecAdd(b.vel, scalMult(dt, b.acc));
        b.pos = vecAdd(b.pos, scalMult(dt, b.vel));
    }
}