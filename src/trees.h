#ifndef TREES
#define TREES

#include <vector>
#include <tuple>
#include "bodies.h"

using namespace std;

// node could be a group of particles or a single particle
struct node{
	double mass{};
	vector<double> pos;

	vector<double> width;
	vector<double> centre;

	int num{}; // number of particles in tree
	int childIndx{};
	vector<int> liveChildren;
	vector<node*> children;
	node* parent{};

	int bodyindx{};

	node(vector<double> width, vector<double>& centre);
	node(node* tree, int chldIndx);
	explicit node(node* root);
	node();
//	~node();
};

// Macros to retrieve node data; x is a pointer
#define Mass(x) (((node*) (x))->mass)
#define pos(x) (((node*) (x))->pos)

class barnesHut{
    node* whichChild(node* tree, int i);
    void addChildren(node*);
    bool inNode(const vector<double>&, node*);
    void treeInsert(node*, int);
public:
    vector<body>* bodies;
    node* root;
    double theta = 0.9;
    double G = 6.674e-11;
    vector<double> width;
    vector<double> centre;

    // Tree building functions
    void treeBuild();
    void treeChop(node*);
    vector<vector<int>> segment(node* root, vector<body> bodies);
    void updateRoot();

    // Kinematic functions
    void acceleration(node*);
    vector<double> treeAcc(node*, int);
    vector<double> ngl(vector<double>& r1, vector<double>& r2, double mass, double softening);

    explicit barnesHut(vector<body>& bods, vector<double> dim);
    explicit barnesHut(vector<body>& bods, vector<double>& dim, vector<double>& cent);
};

// Function definitions



#endif //TREES