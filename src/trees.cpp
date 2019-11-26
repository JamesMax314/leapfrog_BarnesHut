#include <cmath>
#include "vecMaths.h"
#include "trees.h"

vector<body> bodies;
using namespace barnesHut;

// Node functions
bool inNode(vector<double> pos, node* nod){
	for(int i=0; i<3; i++){
	    vector<double> displacement = vecAdd(pos, scalMult(-1, nod->centre));
	    double distance = modulus(displacement, false);
		if (distance > pow(nod->width, 2)){
			return false;
		};
	};
	return true;
};

void addChildren(node* tree){
	tree->children = new node[8]();
    tree->numChildren = 8;
    for(int i=0; i<8; i++){
	    tree->children[i] = node(tree);
	}
	int child = 0;
	for(int i=-1; i<1; i=i+2){
		for(int j=-1; j<1; j=j+2){
			for(int k=-1; k<1; k=k+2){
			    double w = tree->width;
				vector<double> shift{i*w, j*w, k*w};
				tree->children[child].centre = vecAdd(tree->children[child].centre, shift);
				child += 1;
			};
		};
	};
};

void treePrune(node* tree){
	if (tree->children){
		node* child = tree->children;
		for(int i=0; i<8; i++){
			treePrune(child);
			child++;
		};
	} else if (tree->num == 0){
	    tree->parent->numChildren -= 1;
		delete tree;
	};
};

void treeChop(node* tree){
	if (tree->children){
		node* child = tree->children;
		for(int i=0; i<8; i++){
			treeChop(child);
			child++;
		};
	} else{
	    tree->parent->numChildren -= 1;
		delete tree;
	};
};


// Tree building functions
void treeBuild(node* root){ // double width, double centre){
	//node* tree = new node(width, centre);
	addChildren(root);
	for(int i=0; i<bodies.size(); i++){
		treeInsert(root, i);
	}
	treePrune(root);
}

void treeInsert(node* tree, int i){
	if (tree->num > 1){
		node* child = tree->children;
		for(int j=0; j<8; j++) {
            if (inNode(bodies[i].pos, child)) {
                treeInsert(child, i);

                tree->num += 1;

                // Update CM and mass of node
                tree->pos = vecAdd(scalMult(tree->mass, tree->pos), scalMult(bodies[i].mass, bodies[i].pos));
                tree->mass += bodies[i].mass;
                tree->pos = scalMult(1 / tree->mass, tree->pos);
                break;
            }
            child++;
        }
	} else if (tree->num == 1){
		addChildren(tree);
		int primarybody = tree->bodyindx;
		treeInsert(tree, primarybody);
		treeInsert(tree, i);

		tree->num += 1;

		// Update CM and mass of node
        tree->pos = vecAdd(scalMult(tree->mass, tree->pos), scalMult(bodies[i].mass, bodies[i].pos));
		tree->mass += bodies[i].mass;
        tree->pos = scalMult(1 / tree->mass, tree->pos);
	} else{
		tree->bodyindx = i;
		tree->num = 1;
		tree->pos = bodies[i].pos;
		tree->mass += bodies[i].mass;
	};
};


// Kinematic functions
void acceleration(node* tree){
	for(int i=0; i<bodies.size(); i++){
		bodies[i].acc = treeAcc(tree, i);
	};
};

vector<double> treeAcc(node* tree, int i){
	vector<double> f;
	if (tree->num == 1 && tree->bodyindx != i){
		f = ngl(bodies[i].pos, tree->pos, tree->mass);
	} else{
		double distance = modulus(vecAdd(scalMult(-1, bodies[i].pos), tree->centre), false);
		if (distance / (tree->width*2) < theta){
			f = ngl(bodies[i].pos, tree->pos, tree->mass);
		} else{
			node* child = tree->children;
			for(int j=0; j<8; j++){
				f = vecAdd(f, treeAcc(child, i));
				child++;
			}
		}
	}
}

vector<double> ngl(vector<double> &r1, vector<double> &r2, double mass){
	vector<double> out;
	vector<double> delta = vecAdd(scalMult(-1, r1), r2);
	out = scalMult(-mass*G*pow(modulus(delta, false), 3), delta);
	return out;
};