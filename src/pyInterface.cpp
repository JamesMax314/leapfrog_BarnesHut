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

progbar::progbar(double maximum, double barLength): max(maximum), length(barLength){}

void progbar::update(int current){
//    cout << "\r";
//    for (int i=0; i<length+20; i++)
//        cout << " ";
    cout << "\rProgress: ";
    int bars = int(length*((current+1) / max));
    for (int i=0; i<bars; i++)
        cout << "#";
    for (int i=0; i<length-bars; i++)
        cout << "_";
    cout << " | (" << int(100*((current+1) / max)) << "%)" << "\r";
}

vector<body> basicRun(vector<body>& bodies, vector<double> centre, vector<double> width, int numIter, double dt){
    barnesHut bh = barnesHut(bodies, width, centre);
//    bh.theta = 0;
//    bh.G = bh.G*5;
    progbar prog = progbar(numIter, 20);
    for(int j=0; j<numIter; j++) {
//        cout << "Making tree" << endl;
        treeMake(bh);
//        cout << "tree made" << endl;
//        if (j == 1)
//            printTree(bh.root, 0);
//            cout << endl;
//        cout << "computing interactions" << endl;
        interaction(bh);
//        cout << "updating bodies" << endl;
        bodiesUpdate(bh, dt);
//        cout << "breaking tree" << endl;
        treeBreak(bh);
//        printTree(bh.root, 0);
        prog.update(j);
    }
    cout << endl;
    return *bh.bodies;
}

vector<body> fixedBoundary(vector<body>& bodies, vector<body>& boundary, vector<double> centre, vector<double> width, int numIter, double dt){
    barnesHut bh = barnesHut(bodies, width, centre);
    cout << "bodLen: " << boundary.size() << endl;
    progbar prog = progbar(numIter, 20);
    for(int j=0; j<numIter; j++) {
        treeMake(bh);
        interaction(bh);
        boundaryInteract(bh, boundary);
        bodiesUpdate(bh, dt);
        treeBreak(bh);
        prog.update(j);
    }
    cout << endl;
    return *bh.bodies;
}

PYBIND11_MODULE(treecode, m) {
    //m.def("test", &test);
    m.def("basicRun", &basicRun);
    m.def("fixedBoundary", &fixedBoundary);
    py::class_<barnesHut>(m, "barnesHut")
            .def(py::init<vector<body>&, vector<double>&, vector<double>&>());
    py::class_<body>(m, "body")
            .def(py::init<>())
            .def(py::init<double&, vector<double>&,
                    vector<double>&, vector<double>&>())

            .def_property("pos", &body::getAcc, &body::setAcc)
            .def_property("pos", &body::getVel, &body::setVel)
            .def_property("pos", &body::getMass, &body::setMass)
            .def_property("pos", &body::getSoftening, &body::setSoftening)
            .def_property("pos", &body::getPos, &body::setPos);
}