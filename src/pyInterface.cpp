#include <iostream>
#include <fstream>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "bodies.h"
#include "trees.h"
#include "leapfrog.h"
#include "pyInterface.h"
#include "treeShow.h"

namespace py = pybind11;

vector<body> basicRun(vector<body>& bodies, vector<double> centre, vector<double> width, int numIter, double dt){
    barnesHut bh = barnesHut(bodies, width, centre);
    for(int j=0; j<numIter; j++) {
//        cout << "Making tree" << endl;
        treeMake(bh);
//        cout << "tree made" << endl;
        if (j == 0)
            printTree(bh.root, 0);
//        cout << "computing interactions" << endl;
        interaction(bh);
//        cout << "updating bodies" << endl;
        bodiesUpdate(bh, dt);
//        cout << "breaking tree" << endl;
        treeBreak(bh);
    }
    return *bh.bodies;
}

//void test(body& b){
//    body* d = &b;
//    cout << (*d).mass[0] << endl;
//}

PYBIND11_MODULE(treecode, m) {
    //m.def("test", &test);
    m.def("basicRun", &basicRun);
    py::class_<barnesHut>(m, "barnesHut")
            .def(py::init<vector<body>&, vector<double>&, vector<double>&>());
    py::class_<body>(m, "body")
            .def(py::init<>())
            .def(py::init<double&, vector<double>&,
                    vector<double>&, vector<double>&>())
            .def("setPos", &body::setPos)
            .def("setAcc", &body::setAcc)
            .def("setVel", &body::setVel)
            .def("setMass", &body::setMass)
            .def("getPos", &body::getPos)
            .def("getAcc", &body::getAcc)
            .def("getVel", &body::getVel)
            .def("getMass", &body::getMass);
}