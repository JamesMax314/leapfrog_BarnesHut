#include <vector>
#include <tuple>
#include "bodies.h"

#ifndef TREES
#define TREES

using namespace std;

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

// Macros to retrieve node data; x is a pointer
#define Mass(x) (((node*) (x))->mass)
#define pos(x) (((node*) (x))->pos)

class barnesHut{
    void addChildren(node*);
    bool inNode(vector<double>, node*);
    void treePrune(node*);
    void treeInsert(node*, int);
public:
    vector<body> bodies;
    node* root;
    double theta = 0.9;
    double G = 1.6e-11;
    double width;
    vector<double> centre;

    // Tree building functions
    void treeBuild();
    void treeChop(node*);

    // Kinematic functions
    void acceleration(node*);
    vector<double> treeAcc(node*, int);
    vector<double> ngl(vector<double>&, vector<double>&, double);

    explicit barnesHut(vector<body>&, vector<double>);
};

// Function definitions



#endif //TREES