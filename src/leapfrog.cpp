#include <vector>
#include "bodies.h"
#include "trees.h"

node* treeMake(body* bodies, double theta, vector<double> dim){
	double width = dim[0];
	double centre = dim[0]/2;
	node* tree = new node(width, centre);
	return tree;
};

void treeBreak(node* tree){
	treeChop(tree);
}

void bodiesUpdate(node* tree, double dt){
	if (tree->children){
		node* child = tree->children;
		for(int i=0; i<; i++){
			bodiesUpdate(child, dt);
			child++;
		};
	} else{
		body* b = &bodyptr[tree->bodyindx];
		b->vel = vecAdd(b->vel, scalMult(dt, b->acc));
		b->pos = vecAdd(b->pos, scalMult(dt, b->vel));
	};
};