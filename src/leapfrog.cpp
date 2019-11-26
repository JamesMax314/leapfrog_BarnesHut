#include <vector>
#include "bodies.h"
#include "trees.h"
#include "vecMaths.h"

node* treeMake(vector<body> &bodies, vector<double> &dim){
	double width = dim[0];
	vector<double> centre = {dim[0]/2, dim[0]/2};
	node* tree = new node(width, centre);
	return tree;
};

void treeBreak(node* tree){
	treeChop(tree);
}

void bodiesUpdate(node* tree, vector<body> &bodies, double dt){
	if (tree->children){
		node* child = tree->children;
		for(int i=0; i<tree->numChildren; i++){
			bodiesUpdate(child, bodies, dt);
			child++;
		}
	} else{
		body b = bodies[tree->bodyindx];
		b.vel = vecAdd(b.vel, scalMult(dt, b.acc));
		b.pos = vecAdd(b.pos, scalMult(dt, b.vel));
	}
}