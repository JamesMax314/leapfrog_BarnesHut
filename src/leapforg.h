#ifndef LEAPFROG
#define LEAPFROG

#include "bodies.h"
#include "trees.h"

node* treeMake(body*, double, vector<double>);
void treeBreak(node*);
void bodiesUpdate(node*, double);

#endif //LEAPFROG