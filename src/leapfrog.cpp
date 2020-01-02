#include <vector>
#include <algorithm>
#include <iostream>
#include <omp.h>

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
    if (tree->liveChildren.size() != 0){
        for(auto i: tree->liveChildren){
            partialUpdate(tree->children[i], bodies, dt);
        }
    } else{
        auto i = tree->bodyindx;
        (*bodies)[i].vel.emplace_back(vecAdd((*bodies)[i].vel.back(), scalMult(dt, (*bodies)[i].acc.back())));
        (*bodies)[i].pos.emplace_back(vecAdd((*bodies)[i].pos.back(), scalMult(dt, (*bodies)[i].vel.back())));
    }
}

void boundaryInteract(barnesHut& bh, vector<body>& boundary){
    vector<body>* bodies = bh.bodies;
#pragma omp parallel for
    for(auto body : *bodies){
#pragma omp parallel for
        for(auto boundBod : boundary){
            vector<double> dist = vecAdd(body.pos.back(), scalMult(-1, boundBod.pos.back()));
            if (abs(dist[0]) > bh.width[0] && abs(dist[1]) > bh.width[1] && abs(dist[2]) > bh.width[2]) {
                body.acc.back() = vecAdd(body.acc.back(),
                                         bh.ngl(body.pos.back(),
                                                 boundBod.pos.back(), boundBod.mass.back(), 0));
            }
        }
    }
}