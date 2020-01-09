#ifndef LEAPFROG
#define LEAPFROG

#include "bodies.h"
#include "trees.h"

void treeMake(barnesHut& hut);
void interaction(barnesHut& hut);
void treeBreak(barnesHut& hut);
void bodiesUpdate(vector<body>* bodies, double dt);
void partialUpdate(node* tree, vector<body>* bodies, double dt);
void boundaryInteract(barnesHut& bh, vector<body>& boundary);

#endif //LEAPFROG