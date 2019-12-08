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

//node::~node(){
//    cout << "deconstructing" << endl;
//    for (node* child : children)
//        delete child;
//}

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
    vector<double> displacement = vecAdd(pos, scalMult(-1, nod->centre));
    for (int j=0; j<3; j++) {
        if (abs(displacement[j]) > nod->width[j]/2) {
            return false;
        } else if (abs(displacement[j]) == nod->width[j]){
            vector<double> nPos = vecAdd(pos, scalMult(nod->width[j]/10, pos));
            return inNode(nPos, nod);
        }
    }
	return true;
}

void barnesHut::addChildren(node* tree){
	tree->children = vector<node*>(8);
    //cout << "Children initialized at: " << tree->children[0] << endl;
    tree->numChildren = 0;
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

//void barnesHut::treePrune(node* tree){
//    //cout << "tree->children: " << tree->children << endl;
//    //cout << "tree parent: " << tree->parent << endl;
//    if (tree->numChildren != 0){
//		for(int i=0; i<8; i++){
//			treePrune(tree->children[i]);
//		}
//	} else if (tree->num == 0 && tree->parent != nullptr){
////        cout << "liveChildren: " << endl;
//        //cout << "not null tree parent: " << tree->parent << endl;
//        //for(int i : tree->parent->liveChildren){
//        //    cout << i << ", ";
//        //}
//        //cout << endl;
//        tree->parent->numChildren -= 1;
//        //cout << "Erasing vector: " << tree->childIndx << endl;
//        auto indx = find(tree->parent->liveChildren.begin(), tree->parent->liveChildren.end(), tree->childIndx);
//        //cout << "index: " << *indx << endl;
//        tree->parent->liveChildren.erase(indx);
////        cout << "number of root children: " << root->numChildren << endl;
//		delete tree;
//        //cout << "ok" << endl;
//	}
//}

void barnesHut::treeChop(node* tree){
	if (tree->numChildren != 0){
	    vector<int> children = tree->liveChildren;
	    /// Deallocate memory
		for(int i : children){
			treeChop(tree->children[i]);
		}
	}
	if(tree->parent != nullptr){
		delete tree;
	} else{
        tree->numChildren = 0;
        tree->liveChildren = {};
        tree->num = 0;
	}
}

vector<vector<int>> barnesHut::segment(node* nod, vector<body> bods){
    vector<vector<int>> out = {{}, {}, {}, {}, {}, {}, {}, {}};
    if (!bods.empty()) {
        if (nod->liveChildren.empty()) {
            addChildren(nod);
        }
        for (int i = 0; i < bods.size(); i++) {
            for (int j=0; j<nod->children.size(); j++){
                if (inNode(bods[i].pos.back(), nod->children[j]) && bods[i].active[0])
                    out[j].emplace_back(i);
            }
        }
    }
//    cout << "segments: {";
//    for (int i=0; i<nod->children.size(); i++){
//        cout << "{";
//        for (int j=0; j<out[i].size(); j++){
//            cout << out[i][j];
//            if (j<out[i].size()-1)
//                cout << ", ";
//        }
//        cout << "}";
//        if (i<nod->children.size()-1)
//            cout << ", ";
//    }
//    cout << "}" << endl;
    return out;
}

void barnesHut::updateRoot(){
    vector<double> cm = {0.,0.,0.};
    for (auto i : root->liveChildren){
        root->num += root->children[i]->num;
        root->mass += root->children[i]->mass;
        cm = vecAdd(cm, scalMult(root->children[i]->mass, root->children[i]->pos));
    }
    root->pos = scalMult(1/root->mass, cm);
    root->numChildren = int(root->liveChildren.size());
}

// Tree building functions
void barnesHut::treeBuild(){ // double width, double centre){
//    cout << "building..." << endl;

    if(!root){
        root = new node(width, centre);
    }
    vector<vector<int>> segments = segment(root, (*bodies));
//#pragma omp parallel for
    vector<thread> threads(segments.size());
    for (int i=0; i<segments.size(); i++) {
        threads[i] = thread([&, i]() {
//            cout << "i: " << i << endl;
            for (int j : segments[i]) {
//                cout << "j: " << j << endl;
                if (inNode((*bodies)[j].pos.back(), root->children[i])) {
                    treeInsert(root->children[i], j);
                }
            }
        });
        //printTree(root->children[i], 0);
    }
    for (auto&& thrd : threads) {
        thrd.join();
    }
    updateRoot();
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
        /// Update current pos and COM
        current->num += 1;
        current->pos = vecAdd(scalMult(current->mass, current->pos),
                              scalMult((*bodies)[i].mass.back(), (*bodies)[i].pos.back()));
        current->mass += (*bodies)[i].mass.back();
        current->pos = scalMult(1 / current->mass, current->pos);
        if (current->num == 2){
//            cout << "Add children to saturated node" << endl;
            /// Add children to saturated node
            addChildren(current);

            /// Move original body
            node* child = whichChild(current, current->bodyindx);
            current->liveChildren.emplace_back(child->childIndx);
            current->numChildren += 1;
            child->num += 1;
            child->bodyindx = current->bodyindx;
            child->pos = (*bodies)[current->bodyindx].pos.back();
    		child->mass = (*bodies)[current->bodyindx].mass.back();
//    		cout << "preliminary body: " << current->bodyindx << endl;
//    		cout << "new body: " << i << endl;

    		/// Update current
            current = whichChild(current, i);
        } else if (current->num > 2) {
            /// Children already exist so update current
//            cout << "Children already exist so update current" << endl;
            current = whichChild(current, i);
        }
    }

    /// Fill in leaf data
    current->num += 1;
    current->bodyindx = i;
    current->pos = (*bodies)[i].pos.back();
    current->mass = (*bodies)[i].mass.back();
    if (current->parent != nullptr) {
        current->parent->liveChildren.emplace_back(current->childIndx);
        current->numChildren += 1;
    }

//    printTree(root, 0);
}

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
//#pragma omp parallel for
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
		double distance = m_modulus(vecAdd(scalMult(-1, (*bodies)[i].pos.back()), tree->pos), false);
//		cout << "distance: " << distance << endl;
		vector<double> w = tree->width;
		double minWidth = *min_element(w.begin(), w.end());
//		cout << "checking approx..." << endl;
		if (minWidth/distance > theta){
//		    cout << "computing force on " << i << " from bulk node" << endl;
//		    cout << "tree->pos: ";
//		    printVec(tree->pos);
			f = ngl((*bodies)[i].pos.back(), tree->pos, tree->mass);
		} else{
//#pragma omp parallel for
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
	out = scalMult(mass*G/pow(m_modulus(delta, false), 3), delta);
//    cout << "ngl out: " << out[0] << endl;
	return out;
}