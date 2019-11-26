#include <vector>
#include <tuple>
#include "bodies.h"

#ifndef TREES
#define TREES

using namespace std;

// General variables
#define G 1.6e-11;
//#define theta 1;

// node could be a group of particles or a single particle
struct node{
	double mass{};
	vector<double> pos;

	double width{};
	vector<double> centre;

	int num{}; // number of particles in tree
	int numChildren{};
	node* children{};
	node* parent{};

	int bodyindx{};

	node(double, vector<double>&);
	explicit node(node*);
	node();
};

node::node() = default;

// node constructors
node::node(double w, vector<double> &c){
	num = 0;
	numChildren = 0;
	width = w;
	centre = c;
};

node::node(node* tree){
	num = 0;
	numChildren = 0;
	parent = tree;
	children = nullptr;
	width = tree->width/2;
};

// Macros to retrieve node data; x is a pointer
#define Mass(x) (((node*) (x))->mass)
#define pos(x) (((node*) (x))->pos)

class barnesHut{
public:

}

// Function definitions

// Node functions
void addChildren(node*);
bool inNode(vector<double>, node*);

// Tree building functions
void treePrune(node*);
void treeBuild();
void treeInsert(node*, int);
void treeChop(node*);

// Kinematic functions
void acceleration(node*);
vector<double> treeAcc(node*, int);
vector<double> ngl(vector<double>&, vector<double>&, double);

#endif //TREES