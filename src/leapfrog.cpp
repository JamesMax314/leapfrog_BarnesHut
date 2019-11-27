#include <vector>
#include <algorithm>
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
    if (tree->numChildren != 1){
        for(int i=0; i<tree->numChildren; i++){
            partialUpdate(tree->children[tree->liveChildren[i]], bodies, dt);
        }
    } else{
        body b = (*bodies)[tree->bodyindx];
        b.vel = vecAdd(b.vel, scalMult(dt, b.acc));
        b.pos = vecAdd(b.pos, scalMult(dt, b.vel));
    }
}