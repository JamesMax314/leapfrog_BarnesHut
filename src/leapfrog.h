#ifndef LEAPFROG
#define LEAPFROG

#include "bodies.h"
#include "trees.h"

void treeMake(barnesHut&);
void treeBreak(barnesHut&);
void bodiesUpdate(barnesHut&, double);
void partialUpdate(node*, vector<body>&, double);

#endif //LEAPFROG