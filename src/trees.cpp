#include <cmath>
#include <algorithm>
#include <iostream>
#include <omp.h>
#include <thread>
#include "vecMaths.h"
#include "trees.h"
#include "treeShow.h"

/// node constructors

node::node(){
    pos = vector<double>(3);
    width = vector<double>(3);
    centre = vector<double>(3);
}

// node constructors
node::node(vector<double> w, vector<double> &c){
    num = 0;
    width = w;
    centre = c;
    parent = nullptr;
}

node::node(node* tree){
    num = 0;
    parent = tree;
    pos = vector<double>(3);
    width = vector<double>(3);
    centre = {0, 0, 0};
    for(int i=0; i<3; i++){
        width[i] = tree->width[i]/2;
    }
}

node::node(node* tree, int chldIndx){
    num = 0;
    childIndx = chldIndx;
    parent = tree;
    pos = vector<double>(3);
    width = vector<double>(3);
    centre = {0, 0, 0};
    for(int i=0; i<3; i++){
        width[i] = tree->width[i]/2;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// barnesHut constructor

barnesHut::barnesHut(vector<body> &bods, vector<double> dim): bodies(&bods) {
    width = dim;
    centre = {dim[0]/2, dim[1]/2, dim[2]/2};
    root = new node(width, centre);
}

barnesHut::barnesHut(vector<body> &bods, vector<double> &dim, vector<double> &cent): bodies(&bods){
    width = dim;
    centre = cent;
    root = new node(width, centre);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// barnesHut functions

void barnesHut::initActiveBods(){
    vector<int> tmp;
    for (int i=0; i<(*bodies).size(); i++){
        tmp.emplace_back(i);
    }
    activeBods = tmp;
}

bool barnesHut::inNode(const vector<double>& pos, node* nod){
    vector<double> displacement = vecAdd(pos, scalMult(-1, nod->centre));
    for (int j=0; j<3; j++) {
        /// Check in node for each coordinate
        if (abs(displacement[j]) > nod->width[j]/2) {
            return false;
        } else if (abs(displacement[j]) == nod->width[j]){
            vector<double> nPos = vecAdd(pos, scalMult(nod->width[j]/10, pos));
            return inNode(nPos, nod);
        }
    }
	return true;
}

node* barnesHut::whichChild(node* tree, int i){
    node* child = nullptr;
    for(int j=0; j<8; j++) {
        /// For each child in parent
        child = tree->children[j];
        if (inNode((*bodies)[activeBods[i]].pos.back(), child)) {
            /// activeBods[i] replaced with i
            /// If in child j
            break;
        }
    }
    if (!child)
        cout << "not in child" << endl;
    return child;
}


/// Tree building functions
void barnesHut::treeBuild(){
    if(!root){
        root = new node(width, centre);
    }
//    vector<vector<int>> segments = segment(root, (*bodies));
//#pragma omp parallel for
    for(int i=0; i<int(activeBods.size()); i++){
        if(inNode((*bodies)[activeBods[i]].pos.back(), root) && (*bodies)[activeBods[i]].active.back()){
            treeInsert(root, i);
        }
    }
//    updateRoot();
}

/// Initialise children for parent node
void barnesHut::addChildren(node* tree){
    tree->children = vector<node*>(8);
    for(int i=0; i<8; i++){
        tree->children[i] = new node(tree, i);
    }
    int child = 0;
    vector<double> w = tree->width;
    for(int i=-1; i<=1; i=i+2){
        for(int j=-1; j<=1; j=j+2){
            for(int k=-1; k<=1; k=k+2){
                /// Setup centre of child
                vector<double> shift{i*w[0]/4, j*w[1]/4, k*w[2]/4};
                tree->children[child]->centre = vecAdd(tree->centre, shift);
                child += 1;
            }
        }
    }
}

/// Segment bodies according to octant !!! Contains legacy code; no activeBods !!!
vector<vector<int>> barnesHut::segment(node* nod, vector<body> bods){
    vector<vector<int>> out = {{}, {}, {}, {}, {}, {}, {}, {}};
    if (!bods.empty()) {
        if (nod->liveChildren.empty()) {
            addChildren(nod);
        }
        for (int i = 0; i < bods.size(); i++) {
            for (int j=0; j < nod->children.size(); j++){
                /// for each body check each child node
                if (inNode(bods[i].pos.back(), nod->children[j]) && bods[i].active[0])
                    out[j].emplace_back(i);
            }
        }
    }
/*
    cout << "segments: {";
    for (int i=0; i<nod->children.size(); i++){
        cout << "{";
        for (int j=0; j<out[i].size(); j++){
            cout << out[i][j];
            if (j<out[i].size()-1)
                cout << ", ";
        }
        cout << "}";
        if (i<nod->children.size()-1)
            cout << ", ";
    }
    cout << "}" << endl;
*/
    return out;
}

/// Recursive insertion of bodies into tree
void barnesHut::treeInsert(node* tree, int i){
    node* current = tree;
    while (current->num !=0){
        /// Update current pos and COM
        current->num += 1;
        current->pos = vecAdd(scalMult(current->mass, current->pos),
                              scalMult((*bodies)[activeBods[i]].mass.back(), (*bodies)[activeBods[i]].pos.back()));
        current->mass += (*bodies)[activeBods[i]].mass.back();
        current->pos = scalMult(1 / current->mass, current->pos);
        if (current->num == 2){
            /// Add children to saturated node
            addChildren(current);

            /// Move original body
            node* child = whichChild(current, current->bodyindx);
            current->liveChildren.emplace_back(child->childIndx);
            child->num += 1;
            child->bodyindx = current->bodyindx;
            child->pos = (*bodies)[activeBods[current->bodyindx]].pos.back();
    		child->mass = (*bodies)[activeBods[current->bodyindx]].mass.back();

    		/// Update current
            current = whichChild(current, i);
        } else if (current->num > 2) {
            /// Children already exist so update current
            current = whichChild(current, i);
        }
    }

    /// Fill in leaf data
    current->num += 1;
    current->bodyindx = activeBods[i];
    current->pos = (*bodies)[activeBods[i]].pos.back();
    current->mass = (*bodies)[activeBods[i]].mass.back();
    if (current->parent != nullptr) {
        current->parent->liveChildren.emplace_back(current->childIndx);
    }
//    printTree(root, 5);
}

/// Fill in root node
void barnesHut::updateRoot(){
    vector<double> cm = {0.,0.,0.};
    for (auto i : root->liveChildren){
        root->num += root->children[i]->num;
        root->mass += root->children[i]->mass;
        cm = vecAdd(cm, scalMult(root->children[i]->mass, root->children[i]->pos));
    }
    root->pos = scalMult(1/root->mass, cm);
}

/// Clear heap for next iteration
void barnesHut::treeChop(node* tree){
    if (tree->liveChildren.size() != 0){
        vector<int> children = tree->liveChildren;
        /// Deallocate memory
        for(int i : children){
            treeChop(tree->children[i]);
        }
    }
    if(tree->parent != nullptr){
        delete tree;
    } else{
        tree->liveChildren = {};
        tree->num = 0;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Kinematic functions

void barnesHut::acceleration(node* tree){
    for(int i=0; i<int(activeBods.size()); i++){
        (*bodies)[activeBods[i]].acc.emplace_back(treeAcc(tree, i));
    }
}

/// Tree traversal for acceleration
vector<double> barnesHut::treeAcc(node* tree, int i){
	vector<double> f(3, 0);
	if (tree->num == 1 && tree->bodyindx != i){
	    /// Particle-particle interaction
		f = ngl((*bodies)[activeBods[i]].pos.back(), tree->pos, tree->mass,
		        (*bodies)[activeBods[i]].softening + (*bodies)[activeBods[tree->bodyindx]].softening);
	} else{
		double distance = m_modulus(vecAdd(scalMult(-1, (*bodies)[activeBods[i]].pos.back()), tree->pos), false);
		vector<double> w = tree->width;
		double minWidth = *min_element(w.begin(), w.end());
        /// Checking approximation...
		if (minWidth/distance < theta){
		    /// Bulk node computation
			f = ngl((*bodies)[activeBods[i]].pos.back(), tree->pos, tree->mass,
			        (*bodies)[activeBods[i]].softening + (*bodies)[activeBods[tree->bodyindx]].softening);
		} else{
		    vector<double> fs = {0, 0, 0};
            for(int j=0; j<tree->liveChildren.size(); j++){
                /// For each sub node
				f = vecAdd(f, treeAcc(tree->children[tree->liveChildren[j]], i));
			}
		}
	}
    return f;
}

/// Newtons Gravitational Law
vector<double> barnesHut::ngl(vector<double> &r1, vector<double> &r2, double mass, double softening){
	vector<double> out;
	vector<double> delta = vecAdd(scalMult(-1, r1), r2);
	out = scalMult(mass*G/pow(pow(m_modulus(delta, false), 2.) + pow(softening, 2.), 3./2.), delta);
	return out;
}