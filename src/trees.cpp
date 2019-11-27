#include <cmath>
#include <algorithm>
#include <iostream>
#include "vecMaths.h"
#include "trees.h"

/// node constructors
node::node(){
    pos = vector<double>(3);
    width = vector<double>(3);
    centre = vector<double>(3);
}

// node constructors
node::node(vector<double> w, vector<double> &c){
    num = 0;
    numChildren = 0;
    width = w;
    centre = c;
}

node::node(node* tree){
    num = 0;
    numChildren = 0;
    parent = tree;
    children = nullptr;
    pos = vector<double>(3);
    width = vector<double>(3);
    centre = {0, 0, 0};
    for(int i=0; i<3; i++){
        width[i] = tree->width[i]/2;
    }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// barnesHut constructor
barnesHut::barnesHut(vector<body> &bods, vector<double> dim): bodies(bods) {
    width = dim;
    centre = {dim[0]/2, dim[1]/2, dim[2]/2};
    root = new node(width, centre);
}

barnesHut::barnesHut(vector<body> &bods, vector<double> &dim, vector<double> &cent): bodies(bods){
    width = dim;
    centre = cent;
    root = new node(width, centre);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// barnesHut functions
bool barnesHut::inNode(const vector<double>& pos, node* nod){
    cout << "position, node centre: " << pos[0] << " " << nod->centre[0] << endl;
    cout << "width: " <<  nod->width[0] << endl;
    vector<double> displacement = vecAdd(pos, scalMult(-1, nod->centre));
    cout << "displacement: ";
    for(int k=0; k<3; k++){
        cout << displacement[k] << ", ";
    }
    cout << endl;
    double distance = modulus(displacement, false);
    cout << "distance: " << distance << endl;
    for (int j=0; j<3; j++) {
        if (distance >= nod->width[j]) {
            return false;
        }
    }
	return true;
}

void barnesHut::addChildren(node* tree){
	tree->children = new node[8];
	cout << "Children initialized at: " << tree->children << endl;
    tree->numChildren = 8;
    for(int i=0; i<8; i++){
	    tree->children[i] = node(tree);
	}
	int child = 0;
    vector<double> w = tree->width;
    for(int i=-1; i<=1; i=i+2){
		for(int j=-1; j<=1; j=j+2){
			for(int k=-1; k<=1; k=k+2){
				vector<double> shift{i*w[0]/4, j*w[1]/4, k*w[2]/4};
				tree->children[child].centre = vecAdd(tree->children[child].centre, shift);
				child += 1;
			}
		}
	}
}

void barnesHut::treePrune(node* tree){
    //cout << "tree->children: " << tree->children << endl;
	if (tree->numChildren != 0){
		node* child = tree->children;
		for(int i=0; i<8; i++){
			treePrune(child);
			child++;
		}
	} else if (tree->num == 0){
	    tree->parent->numChildren -= 1;
	    cout << "number of root children: " << root->numChildren << endl;
		delete[] tree->parent->children;
		cout << "ok" << endl;
	}
}

void barnesHut::treeChop(node* tree){
	if (tree->children){
		node* child = tree->children;
		for(int i=0; i<8; i++){
			treeChop(child);
			child++;
		}
	} else{
	    tree->parent->numChildren -= 1;
		delete tree;
	}
}


// Tree building functions
void barnesHut::treeBuild(){ // double width, double centre){
    if(!root){
        root = new node(width, centre);
        cout << "root initialised at: " << root << endl;
    }
	//addChildren(root);
	for(int i=0; i<bodies.size(); i++){
		treeInsert(root, i);
		cout << "root num: " << root->num << endl;
		cout << "body "<< i << " has been inserted" << endl;
	}
	cout << "pruning..." << endl;
	treePrune(root);
}

node* barnesHut::whichChild(node* tree, int i){
    node* child = tree->children;
    int ia = 0;
    for(int j=0; j<8; j++) {
        if (inNode(bodies[i].pos, child)) {
            cout << "in child " << j << endl;
            break;
        }
        child++;
        ia += 1;
    }
    cout << "ia: " << ia << endl;
    return child;
}

void barnesHut::treeInsert(node* tree, int i){
    cout << "i: " << i << endl;
    cout << "tree->num: " << tree->num << endl;
	if (tree->num > 1){
	    node* child = whichChild(tree, i);
        cout << "child insertion of i: " << i << endl;
        treeInsert(child, i);

        // Update CM and mass of node
        tree->pos = vecAdd(scalMult(tree->mass, tree->pos), scalMult(bodies[i].mass, bodies[i].pos));
        tree->mass += bodies[i].mass;
        tree->pos = scalMult(1 / tree->mass, tree->pos);

	} else if (tree->num == 1){
	    cout << "adding children" << endl;
		addChildren(tree);
		int primarybody = tree->bodyindx;
        tree->num += 1;
        cout << "pushing primary" << endl;
        treeInsert(tree, primarybody);
		cout << "pushing secondary" << endl;
		treeInsert(tree, i);

		// Update CM and mass of node
        tree->pos = vecAdd(scalMult(tree->mass, tree->pos), scalMult(bodies[i].mass, bodies[i].pos));
		tree->mass += bodies[i].mass;
        tree->pos = scalMult(1 / tree->mass, tree->pos);
	} else{
	    cout << "adding first element" << endl;
		tree->bodyindx = i;
		tree->num = 1;
		tree->pos = bodies[i].pos;
		tree->mass = bodies[i].mass;
	}
}

// Kinematic functions
void barnesHut::acceleration(node* tree){
	for(int i=0; i<bodies.size(); i++){
		bodies[i].acc = treeAcc(tree, i);
	}
}

vector<double> barnesHut::treeAcc(node* tree, int i){
	vector<double> f;
	if (tree->num == 1 && tree->bodyindx != i){
		f = ngl(bodies[i].pos, tree->pos, tree->mass);
	} else{
		double distance = modulus(vecAdd(scalMult(-1, bodies[i].pos), tree->centre), false);
		vector<double> w = tree->width;
		double minWidth = *min_element(w.begin(), w.end());
		if (minWidth/distance < theta){
			f = ngl(bodies[i].pos, tree->pos, tree->mass);
		} else{
			node* child = tree->children;
			for(int j=0; j<8; j++){
				f = vecAdd(f, treeAcc(child, i));
				child++;
			}
		}
	}
    return f;
}

vector<double> barnesHut::ngl(vector<double> &r1, vector<double> &r2, double mass){
	vector<double> out;
	vector<double> delta = vecAdd(scalMult(-1, r1), r2);
	out = scalMult(-mass*G*pow(modulus(delta, false), 3), delta);
	return out;
}