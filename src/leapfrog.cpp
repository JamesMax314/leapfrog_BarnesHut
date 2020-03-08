#include <vector>
#include <algorithm>
#include <iostream>
#include <omp.h>
#include <cmath>

#include "bodies.h"
#include "trees.h"
#include "vecMaths.h"
#include "leapfrog.h"
#include "tpm.h"

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
#pragma omp parallel for
    for (int i=0; i<activeBods.size(); i++) {
        auto bIndx = activeBods[i];
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
        (*bodies)[bIndx].time.emplace_back((*bodies)[bIndx].time.back() + dt);
//        cout << body.vel.back()[0] << endl;
    }
}

void bodiesUpdate(vector<body>* bodies, const vector<int>& activeBods, double t, double dt, vector<double> dim){
    tree_PM tpm_instance = tree_PM();
#pragma omp parallel for
    for (int i=0; i<activeBods.size(); i++) {
        double a_0 = tpm_instance.get_a(t);
        double a_1 = tpm_instance.get_a(t);
        double ad_0 = tpm_instance.get_ad(t);
        double ad_1 = tpm_instance.get_ad(t);
        double add_0 = tpm_instance.get_add(t);
        double add_1 = tpm_instance.get_add(t);
        auto bIndx = activeBods[i];
        /// step 1 in computing the acceleration in cc
        auto acc = vecAdd(scalMult(1/a_0, (*bodies)[bIndx].acc.back()), scalMult(-ad_0/a_0, (*bodies)[bIndx].vel.back()));
//        Compute the cc acceleration
//        auto ud = vecAdd(scalMult(-ad_0, (*bodies)[bIndx].vel.back()), acc);
        /// Update the cc velocity
        auto xd = vecAdd((*bodies)[bIndx].vel.back(), scalMult(dt, acc));
        /// Compute the real position
        auto x = vecAdd((*bodies)[bIndx].pos.back(), scalMult(dt, xd));
        /// Apply periodic BC's
        for (int axis=0; axis<3; axis++){
            if (x[axis] < -dim[axis]/2) {
                x[axis] += dim[axis];
//                v[axis] *= -1;
            }
            if (x[axis] > dim[axis]/2) {
                x[axis] -= dim[axis];
//                v[axis] *= -1;
            }
        }
        /// Append results to body
        (*bodies)[bIndx].vel.emplace_back(xd);
        (*bodies)[bIndx].pos.emplace_back(x);
        (*bodies)[bIndx].time.emplace_back((*bodies)[bIndx].time.back() + dt);
//        cout << body.vel.back()[0] << endl;
    }
}

void PBC(vector<body>* bodies, const vector<int>& activeBods, vector<double> dim){
#pragma omp parallel for
    for (int i=0; i<activeBods.size(); i++) {
        auto bIndx = activeBods[i];
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
#pragma omp parallel for
    for (int i=0; i<activeBods.size(); i++) {
        auto bIndx = activeBods[i];
        (*bodies)[bIndx].vel.emplace_back(vecAdd((*bodies)[bIndx].vel.back(),
                scalMult(dt, (*bodies)[bIndx].acc.back())));
        (*bodies)[bIndx].pos.emplace_back(vecAdd((*bodies)[bIndx].pos.back(),
                scalMult(dt, (*bodies)[bIndx].vel.back())));
        (*bodies)[bIndx].time.emplace_back((*bodies)[bIndx].time.back() + dt);
    }
}

void bodiesUpdate(vector<body>* bodies, const vector<int>& activeBods, double t, double dt){
    tree_PM tpm_instance = tree_PM();
#pragma omp parallel for
    for (int i=0; i<activeBods.size(); i++) {
        double a_0 = tpm_instance.get_a(t);
        double a_1 = tpm_instance.get_a(t);
        double ad_0 = tpm_instance.get_ad(t);
        double ad_1 = tpm_instance.get_ad(t);
        double add_0 = tpm_instance.get_add(t);
        double add_1 = tpm_instance.get_add(t);
        auto bIndx = activeBods[i];
        /// step 1 in computing the acceleration in cc
        auto acc = vecAdd(scalMult(1/pow(a_0, 3), (*bodies)[bIndx].acc.back()), scalMult(-ad_0/a_0, (*bodies)[bIndx].vel.back()));
//        Compute the cc acceleration
//        auto ud = vecAdd(scalMult(-ad_0, (*bodies)[bIndx].vel.back()), fi);
        /// Update the cc velocity
        auto xd = vecAdd((*bodies)[bIndx].vel.back(), scalMult(dt, acc));
        /// Compute the real position
        auto x = vecAdd((*bodies)[bIndx].pos.back(), scalMult(dt, xd));
        /// Apply periodic BC's
        (*bodies)[bIndx].vel.emplace_back(xd);
        (*bodies)[bIndx].pos.emplace_back(x);
        (*bodies)[bIndx].time.emplace_back((*bodies)[bIndx].time.back() + dt);
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