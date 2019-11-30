#include <iostream>
#include <algorithm>

#include "treeShow.h"
#include "trees.h"
#include "vecMaths.h"

using namespace std;

void addSpace(int space){
    for (int j=COUNT; j<(space); j++)
        cout << " ";
    cout << "|";
}

void printTree(node* nod, int space){
    // Increase distance between levels
    space += COUNT;

    if (nod->parent == nullptr){
        cout << INDENT << "root:" << endl;
        cout << INDENTCLEAR << "~ num: " << nod->num << endl;
        cout << INDENTCLEAR << "~ numChildren: " << nod->liveChildren.size() << endl;
        cout << INDENTCLEAR << "~ width: "; printVec(nod->width); cout << endl;
        cout << INDENTCLEAR << "~ centre: "; printVec(nod->centre); cout << endl;
        cout << INDENTCLEAR << "~ COM: "; printVec(nod->pos); cout << endl;
    } else{
        vector<int> v = nod->parent->liveChildren;
        if (find(v.begin(), v.end(), nod->childIndx) != v.end()) {
            addSpace(space); cout << INDENT << nod->childIndx << ": " << nod->num << endl;
            addSpace(space); cout << INDENTCLEAR << "~ width: "; printVec(nod->width); cout << endl;
            addSpace(space); cout << INDENTCLEAR << "~ centre: "; printVec(nod->centre); cout << endl;
            addSpace(space); cout << INDENTCLEAR << "~ COM: "; printVec(nod->pos); cout << endl;
        }
    }

    // Process right child first
    for(int i=0; i<nod->children.size(); i++) {
        vector<int> v = nod->liveChildren;
        if (find(v.begin(), v.end(), i) != v.end()) {
            printTree(nod->children[i], space);
        } else{
            addSpace(space+COUNT); cout << endl;
            addSpace(space+COUNT); cout << INDENT << i << ": #" << endl;
        }
    }
}