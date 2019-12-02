#include <cmath>
#include <algorithm>
#include <iostream>
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
    numChildren = 0;
    width = w;
    centre = c;
    parent = nullptr;
}

node::node(node* tree){
    num = 0;
    numChildren = 0;
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
    numChildren = 0;
    childIndx = chldIndx;
    parent = tree;
    //tree->liveChildren.push_back(childIndx);
    //cout << "length: " << tree->liveChildren.size() << endl;
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
bool barnesHut::inNode(const vector<double>& pos, node* nod){
    //cout << "position, node centre: " << pos[0] << " " << nod->centre[0] << endl;
    //cout << "width: " <<  nod->width[0] << endl;
//    cout << "nod->centre: "; printVec(nod->centre);
    vector<double> displacement = vecAdd(pos, scalMult(-1, nod->centre));
//    cout << "pos: "; printVec(pos);
//    cout << "displacement: "; printVec(displacement);
//    cout << "width: "; printVec(nod->width);
//    for(int k=0; k<3; k++){
//        cout << displacement[k] << ", ";
//    }
    //cout << endl;
    //double distance = modulus(displacement, false);
    //cout << "distance: " << distance << endl;
    for (int j=0; j<3; j++) {
        if (abs(displacement[j]) > nod->width[j]/2) {
//            cout << "false" << endl;
            return false;
        } else if (abs(displacement[j]) == nod->width[j]){
            vector<double> nPos = vecAdd(pos, scalMult(nod->width[j]/10, pos));
//            cout << "node: " << nod;
            return inNode(nPos, nod);
        }
    }
	return true;
}

void barnesHut::addChildren(node* tree){
	tree->children = vector<node*>(8);
    //cout << "Children initialized at: " << tree->children[0] << endl;
    tree->numChildren = 8;
    for(int i=0; i<8; i++){
        tree->children[i] = new node(tree, i);
    //    cout << "added: " << i << endl;
    }
	int child = 0;
    vector<double> w = tree->width;
    for(int i=-1; i<=1; i=i+2){
		for(int j=-1; j<=1; j=j+2){
			for(int k=-1; k<=1; k=k+2){
				vector<double> shift{i*w[0]/4, j*w[1]/4, k*w[2]/4};
				tree->children[child]->centre = vecAdd(tree->centre, shift);
				child += 1;
			}
		}
	}
}

void barnesHut::treePrune(node* tree){
    //cout << "tree->children: " << tree->children << endl;
    //cout << "tree parent: " << tree->parent << endl;
    if (tree->numChildren != 0){
		for(int i=0; i<8; i++){
			treePrune(tree->children[i]);
		}
	} else if (tree->num == 0 && tree->parent != nullptr){
//        cout << "liveChildren: " << endl;
        //cout << "not null tree parent: " << tree->parent << endl;
        //for(int i : tree->parent->liveChildren){
        //    cout << i << ", ";
        //}
        //cout << endl;
        tree->parent->numChildren -= 1;
        //cout << "Erasing vector: " << tree->childIndx << endl;
        auto indx = find(tree->parent->liveChildren.begin(), tree->parent->liveChildren.end(), tree->childIndx);
        //cout << "index: " << *indx << endl;
        tree->parent->liveChildren.erase(indx);
//        cout << "number of root children: " << root->numChildren << endl;
		delete tree;
        //cout << "ok" << endl;
	}
}

void barnesHut::treeChop(node* tree){
    //cout << "treeChop" << endl;
    //cout << "numChildren:" << tree->numChildren << endl;
	if (tree->numChildren != 0){
	    vector<int> children = tree->liveChildren;
		for(int i : children){
			treeChop(tree->children[i]);
		}
		tree->num = 0;
	} else if(tree->parent != nullptr){
	    tree->parent->numChildren -= 1;
        auto indx = find(tree->parent->liveChildren.begin(), tree->parent->liveChildren.end(), tree->childIndx);
        tree->parent->liveChildren.erase(indx);
		delete tree;
	} else{
	    tree->num = 0;
	}
}


// Tree building functions
void barnesHut::treeBuild(){ // double width, double centre){
//    cout << "building..." << endl;

    if(!root){
//        cout << "initialising root" << endl;
        root = new node(width, centre);
        //cout << "root initialised at: " << root << endl;
    }
	//addChildren(root);
//	cout << "size: " << ((*bodies).size()) << endl;
	for(int i=0; i<int((*bodies).size()); i++){
	    cout << i << " in root: " << inNode((*bodies)[i].pos.back(), root) << endl;
	    if(inNode((*bodies)[i].pos.back(), root) && (*bodies)[i].active[0]){
            treeInsert(root, i);
            cout << "body "<< i << " has been inserted" << endl;
        }
        //cout << i << endl;
        //cout << "root num: " << root->num << endl;
	}
    cout << "pruning..." << endl;
	treePrune(root);
}

node* barnesHut::whichChild(node* tree, int i){
    node* child = nullptr;
    for(int j=0; j<8; j++) {
        child = tree->children[j];
        if (inNode((*bodies)[i].pos.back(), child)) {
//            cout << "in child " << j << endl;
            break;
        }
    }
    if (!child)
        cout << "not in child" << endl;
    return child;
}

void barnesHut::treeInsert(node* tree, int i){
    node* current = tree;
    while (current->num !=0){
        cout << "current: " << current << endl;
        if (current->num == 1){
            cout << "if" << endl;
            addChildren(current);
            node* child = whichChild(current, current->bodyindx);
            current->num += 1;
            current->liveChildren.emplace_back(child->childIndx);
            child->num += 1;
            child->bodyindx = current->bodyindx;
            child->pos = (*bodies)[current->bodyindx].pos.back();
    		child->mass = (*bodies)[current->bodyindx].mass.back();
        } else{
            cout << "else" << endl;
            current = whichChild(current, i);
            cout << "current->num: " << current->num << endl;
            cout << "current index: " << current->childIndx << endl;
            cout << "body centre: "; printVec((*bodies)[i].pos.back()); cout << endl;
        }
        printTree(root, 0);
    }
    current->parent->liveChildren.emplace_back(i);
    current->num += 1;
    current->bodyindx = i;
    current->pos = (*bodies)[i].pos.back();
    current->mass = (*bodies)[i].mass.back();
    printTree(root, 0);
}
//
//void barnesHut::treeInsert(node* tree, int i){
//    cout << "i: " << i << endl;
//    cout << "tree->num: " << tree->num << endl;
//	if (tree->num > 1){
//	    node* child = whichChild(tree, i);
//        cout << "insertion of i: " << i << endl;
//        treeInsert(child, i);
//        cout << "inserted" << endl;
//
//        // Update CM and mass of node
//        tree->pos = vecAdd(scalMult(tree->mass, tree->pos),
//                scalMult((*bodies)[i].mass.back(), (*bodies)[i].pos.back()));
//        tree->mass += (*bodies)[i].mass.back();
//        tree->pos = scalMult(1 / tree->mass, tree->pos);
//        tree->num += 1;
//
//	} else if (tree->num == 1){
//        cout << "adding children" << endl;
//		addChildren(tree);
//		int primarybody = tree->bodyindx;
//        tree->num += 1;
//        cout << "pushing primary" << primarybody << endl;
//        cout << "tree num" << tree->num << endl;
//        treeInsert(tree, primarybody);
//        cout << "pushing secondary" << endl;
//		treeInsert(tree, i);
//        tree->num -= 2;
//
//        // Update CM and mass of node
//        tree->pos = vecAdd(scalMult(tree->mass, tree->pos),
//                scalMult((*bodies)[i].mass.back(), (*bodies)[i].pos.back()));
//		tree->mass += (*bodies)[i].mass.back();
//        tree->pos = scalMult(1 / tree->mass, tree->pos);
//	} else{
//        cout << "adding first element" << endl;
//		tree->bodyindx = i;
//		tree->num = 1;
//		tree->pos = (*bodies)[i].pos.back();
//		tree->mass = (*bodies)[i].mass.back();
//	}
//}

// Kinematic functions
void barnesHut::acceleration(node* tree){
	for(int i=0; i<int((*bodies).size()); i++){
//	    cout << "acc: " << i << endl;
        (*bodies)[i].acc.emplace_back(treeAcc(tree, i));
	}
}

vector<double> barnesHut::treeAcc(node* tree, int i){
//    cout << "treeAcc" << endl;
	vector<double> f(3, 0);
	if (tree->num == 1 && tree->bodyindx != i){
//	    cout << "if" << endl;
		f = ngl((*bodies)[i].pos.back(), tree->pos, tree->mass);
	} else{
//	    cout << "else" << endl;
//	    cout << "pos size: " << (*bodies)[i].pos.size() << endl;
		double distance = modulus(vecAdd(scalMult(-1, (*bodies)[i].pos.back()), tree->pos), false);
//		cout << "distance: " << distance << endl;
		vector<double> w = tree->width;
		double minWidth = *min_element(w.begin(), w.end());
//		cout << "checking approx..." << endl;
		if (minWidth/distance < theta){
//		    cout << "computing force on " << i << " from bulk node" << endl;
//		    cout << "tree->pos: ";
//		    printVec(tree->pos);
			f = ngl((*bodies)[i].pos.back(), tree->pos, tree->mass);
		} else{
			for(int j=0; j<tree->numChildren; j++){
//			    cout << "j: " << j << endl;
//                cout << "tree->liveChildren[j]" << tree->liveChildren[j] << endl;
				f = vecAdd(f, treeAcc(tree->children[tree->liveChildren[j]], i));
			}
		}
	}
//    cout << "f: " << f[0] << endl;
    return f;
}

vector<double> barnesHut::ngl(vector<double> &r1, vector<double> &r2, double mass){
//    cout << "ngl" << endl;
	vector<double> out;
	vector<double> delta = vecAdd(scalMult(-1, r1), r2);
	out = scalMult(mass*G/pow(modulus(delta, false), 3), delta);
//    cout << "ngl out: " << out[0] << endl;
	return out;
}