#ifndef LEAPFROG
#define LEAPFROG

#include "bodies.h"
#include "trees.h"

node* treeMake(vector<body>&, vector<double>&);
void treeBreak(node*);
void bodiesUpdate(node*,vector<body>& ,double);

#endif //LEAPFROG