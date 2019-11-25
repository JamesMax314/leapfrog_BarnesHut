#include <math.h>
#include "vecMath.h"
#include "trees.h"

body* bodyptr;

// Node functions
bool inNode(vector<double> pos, node* nod){
	for(int i=0; i<3; i++){
		if (pow(pos[i]-nod->centre, 2) > pow(nod->width, 2)){
			return false;
		};
	};
	return true;
};

void addChildren(node* tree){
	tree->children = new node(tree)[8];
	int child = 0;
	for(int i=-1; i<1; i=i+2){
		for(int j=-1; i<1; j=j+2){
			for(int k=-1; i<1; k=k+2){
				vector<double> shift = {i, j, k};
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
			treePrune(ptr);
			child++;
		};
	} else if (tree->num == 0){
		delete tree;
	};
};


// Tree building functions
void treeBuild(node* root, ){ // double width, double centre){
	//node* tree = new node(width, centre);
	addChildren(tree);
	for(int i=0; i<n; i++){
		treeInsert(tree, i, root);
	};
	treePrune(tree);
};

void treeInsert(node* tree, int i){
	if (tree->num > 1){
		node* child = tree->children[j];
		for(int j=0; j<8; j++){
			if inNode(bodyptr[i]->pos, child){
				break;
			};
			child++
		};
		treeInsert(tree->children[j], i);

		tree->num += 1;

		// Update CM and mass of node
		tree->pos = tree->pos*tree->mass + bodyptr[i].pos*bodyptr[i].mass;
		tree->mass += bodyptr[i].mass;
		tree->pos = tree->pos / tree->mass;
	} else if (tree->num = 1){
		addChildren(tree);
		primarybody = tree->bodyindx;
		treeInsert(tree, primarybody);
		treeInsert(tree, i);

		tree->num += 1;

		// Update CM and mass of node
		tree->pos = tree->pos*tree->mass + bodyptr[i].pos*bodyptr[i].mass;
		tree->mass += bodyptr[i].mass;
		tree->pos = tree->pos / tree->mass;
	} else{
		tree->bodyindx = i;
		tree->num = 1;
		tree->pos = bodyptr[i].pos;
		tree->mass += bodyptr[i].mass;
	};
};


// Kinematic functions
void acceleration(node* tree){
	for(int i=0; i<n; i++){
		bodyptr[i].acceleration = treeAcc(tree, i);
	};
};

vector<double> treeAcc(node* tree, int i){
	vector<double> f = {0, 0, 0};
	if (tree->num == 1 && tree->bodyindx != i){
		f = ngl(podyptr[i].pos, tree->pos, tree->mass);
	} else{
		double distance = modulus(vecAdd(scalMult(-1, bodyptr[i].pos), tree->centre));
		if (distance / (tree->width*2) < theta){
			f = ngl(podyptr[i].pos, tree->pos, tree->mass);
		} else{
			node* child = tree->children;
			for(int i=0; i<8; i++){
				f = vecAdd(f, treeAcc(child));
				child++;
			};
		};
	};
};

void ngl(vector<double> r1, vector<double> r2, double mass){
	vector<double> out;
	vector<double> delta = vecAdd(scalMult(-1, r1), r2);
	out = scalMult(-G*pow(modulus(deta, 3), delta))
	return out;
};