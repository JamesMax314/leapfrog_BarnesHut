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

//void bodiesUpdate(barnesHut &hut, double dt){
//    node* tree = hut.root;
//    vector<body>* bodies = hut.bodies;
//    partialUpdate(tree, bodies, dt);
//}

void bodiesUpdate(vector<body>* bodies, const vector<int>& activeBods, double dt, vector<double> dim){
    for (auto & bIndx : activeBods) {
        auto v = vecAdd((*bodies)[bIndx].vel.back(), scalMult(dt, (*bodies)[bIndx].acc.back()));
        auto p = vecAdd((*bodies)[bIndx].pos.back(), scalMult(dt, (*bodies)[bIndx].vel.back()));
        for (int axis=0; axis<3; axis++){
            if (p[axis] < -dim[axis]/2) {
                p[axis] += dim[axis];
//                v[axis] *= -1;
            }
            if (p[axis] > dim[axis]/2) {
                p[axis] -= dim[axis];
//                v[axis] *= -1;
            }
        }
        (*bodies)[bIndx].vel.emplace_back(v);
        (*bodies)[bIndx].pos.emplace_back(p);
//        cout << body.vel.back()[0] << endl;
    }
}

void PBC(vector<body>* bodies, const vector<int>& activeBods, vector<double> dim){
    for (auto & bIndx : activeBods) {
        auto p = (*bodies)[bIndx].pos.back();
        for (int axis=0; axis<3; axis++){
            if (p[axis] < -dim[axis]/2) {
                p[axis] += dim[axis];
//                v[axis] *= -1;
            }
            if (p[axis] > dim[axis]/2) {
                p[axis] -= dim[axis];
//                v[axis] *= -1;
            }
        }
        (*bodies)[bIndx].pos.back() = p;
//        cout << body.vel.back()[0] << endl;
    }
}

void bodiesUpdate(vector<body>* bodies, const vector<int>& activeBods, double dt){
    for (auto & bIndx : activeBods) {
        (*bodies)[bIndx].vel.emplace_back(vecAdd((*bodies)[bIndx].vel.back(),
                scalMult(dt, (*bodies)[bIndx].acc.back())));
        (*bodies)[bIndx].pos.emplace_back(vecAdd((*bodies)[bIndx].pos.back(),
                scalMult(dt, (*bodies)[bIndx].vel.back())));
    }
}

//void partialUpdate(node* tree, vector<body>* bodies, double dt){
//    if (tree->liveChildren.size() != 0){
//        for(auto i: tree->liveChildren){
//            partialUpdate(tree->children[i], bodies, dt);
//        }
//    } else{
//        auto i = tree->bodyindx;
//        (*bodies)[i].vel.emplace_back(vecAdd((*bodies)[i].vel.back(), scalMult(dt, (*bodies)[i].acc.back())));
//        (*bodies)[i].pos.emplace_back(vecAdd((*bodies)[i].pos.back(), scalMult(dt, (*bodies)[i].vel.back())));
//    }
//}

void boundaryInteract(barnesHut& bh, vector<body>& boundary){
    vector<body>* bodies = bh.bodies;
    for(auto body : *bodies){
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