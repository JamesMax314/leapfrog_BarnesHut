#include <iostream>
#include <algorithm>
#include <windows.h>

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
        if (nod->num==1)
            SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 11);
        if (nod->num!=1)
            SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 13);
        addSpace(0); cout << INDENT << "root:" << endl;
        addSpace(0); cout << INDENTCLEAR << "~ num: " << nod->num << endl;
        addSpace(0); cout << INDENTCLEAR << "~ numChildren: " << nod->liveChildren.size() << endl;
        addSpace(0); cout << INDENTCLEAR << "~ width: "; printVec(nod->width); cout << endl;
        addSpace(0); cout << INDENTCLEAR << "~ centre: "; printVec(nod->centre); cout << endl;
        addSpace(0); cout << INDENTCLEAR << "~ COM: "; printVec(nod->pos); cout << endl;
    } else{
        vector<int> v = nod->parent->liveChildren;
        if (find(v.begin(), v.end(), nod->childIndx) != v.end()) {
            if (nod->num==1)
                SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 11);
            if (nod->num!=1)
                SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 13);
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
            SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 12);
            addSpace(space+COUNT); cout << endl;
            addSpace(space+COUNT); cout << INDENT << i << ": #" << endl;
        }
    }
}