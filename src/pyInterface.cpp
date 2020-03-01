#include <iostream>
#include <fstream>
#include <chrono>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "bodies.h"
#include "trees.h"
#include "leapfrog.h"
#include "pyInterface.h"
#include "treeShow.h"
#include "poisson.h"
#include "tpm.h"

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
    bh.initActiveBods();
    progbar prog = progbar(numIter, 20);
    for(int j=0; j<numIter; j++) {
        treeMake(bh);
        interaction(bh);
        bodiesUpdate(bh.bodies, bh.activeBods, dt);
        treeBreak(bh);
        prog.update(j);
    }
    cout << endl;
    return *bh.bodies;
}

vector<body> fixedBoundary(vector<body>& bodies, vector<body>& boundary, vector<double> centre,
        vector<double> width, int numIter, double dt){
    barnesHut bh = barnesHut(bodies, width, centre);
    bh.initActiveBods();
    cout << "bodLen: " << boundary.size() << endl;
    progbar prog = progbar(numIter, 20);
    for(int j=0; j<numIter; j++) {
        treeMake(bh);
        interaction(bh);
        boundaryInteract(bh, boundary);
        bodiesUpdate(bh.bodies, bh.activeBods, dt);
        treeBreak(bh);
        prog.update(j);
    }
    cout << endl;
    return *bh.bodies;
}

vector<body> particleMesh(vector<body>& bodies, double spacing, double width, int numIter, double dt){
    vector<body>* bods = &bodies;
    grid g = grid(spacing, width);
    g.initActiveBods(bods);
    progbar prog = progbar(numIter, 20);
    for(int j=0; j<numIter; j++) {
        g.updateGrid(bods);
        g.solveField();
//        g.interp(bods);
        g.interpW(bods, true);
        bodiesUpdate(bods, g.activeBods, dt, g.dim);
//        g.spacing = g.spacing*1.01;
//        for (int i=0; i<3; i++){
//            g.dim[i] = g.dim[i]*1.01;
//        }
        prog.update(j);
    }
    cout << endl;
    return bodies;
}

vector<vector<vector<double>>> PMTest(vector<body>& bodies, double spacing, double width, int axis){
    vector<body>* bods = &bodies;
    grid g = grid(spacing, width);
    g.initActiveBods(bods);
    g.updateGrid(bods);
    g.solveField();
//    g.getF1(axis);
//    cout << g.realField[0][int(100)] << endl;
    auto out = g.getF(axis);
//    g.ctor(g.realPot);
//    g.magF();
//    g.interp(bods);
//    bodiesUpdate(bods, dt);
    return out;
}

vector<body> PMTest1(vector<body>& bodies, double spacing, double width, double dt){
    vector<body>* bods = &bodies;
    grid g = grid(spacing, width);
    g.initActiveBods(bods);
    g.updateGrid(bods);
    g.solveField();
//    g.ctor(g.realPot);
//    g.magF();
    g.interp(bods);
    bodiesUpdate(bods, g.activeBods, dt);
    return bodies;
}

vector<vector<vector<double>>> PMTestPot(vector<body>& bodies, double spacing, double width){
    vector<body>* bods = &bodies;
    grid g = grid(spacing, width);
    g.initActiveBods(bods);
    g.updateGrid(bods);
//    g.solveField();
//    g.ctor(g.realPot);
//    g.magF();
//    g.interp(bods);
//    bodiesUpdate(bods, g.activeBods, dt);
    double* out = new double[int(g.numPts[0]*g.numPts[1]*g.numPts[2])];
    g.ctor(g.realPot, out);
    vector<vector<vector<double>>> vOut;
    for (int i=0; i<g.numPts[0]; i++){
        vector<vector<double>> tmp1;
        for (int j=0; j<g.numPts[0]; j++) {
            vector<double> tmp2;
            for (int k = 0; k < g.numPts[0]; k++) {
                tmp2.emplace_back(out[int(i * g.numPts[2] * g.numPts[1] + j * g.numPts[2] + k)]);
            }
            tmp1.emplace_back(tmp2);
        }
        vOut.emplace_back(tmp1);
    }
    return vOut;
}

vector<vector<vector<double>>> PMTestForce(vector<body>& bodies, double spacing, double width, int axis){
    vector<body>* bods = &bodies;
    grid g = grid(spacing, width);
    g.initActiveBods(bods);
    g.updateGrid(bods);
    g.solveField();
//    g.ctor(g.realPot);
//    g.magF();
//    g.interp(bods);
//    bodiesUpdate(bods, g.activeBods, dt);
    double* out = g.realField[axis];
    vector<vector<vector<double>>> vOut;
    for (int i=0; i<g.numPts[0]; i++){
        vector<vector<double>> tmp1;
        for (int j=0; j<g.numPts[0]; j++) {
            vector<double> tmp2;
            for (int k = 0; k < g.numPts[0]; k++) {
                tmp2.emplace_back(out[int(i * g.numPts[2] * g.numPts[1] + j * g.numPts[2] + k)]);
            }
            tmp1.emplace_back(tmp2);
        }
        vOut.emplace_back(tmp1);
    }
    return vOut;
}

vector<body> TreePareticleMesh(vector<body>& bodies, double spacing, double width,
        double density, int numIter, double dt, double t){
    tree_PM tpm = tree_PM(bodies, spacing, width, density, dt, t);
    progbar prog = progbar(numIter, 20);
    auto t1 = chrono::high_resolution_clock::now();
    auto tT = chrono::high_resolution_clock::now();
    auto genSeeds = 0.;
    auto genSubGrids = 0.;
    auto classiftBods = 0.;
    auto runTrees = 0.;
    auto update = 0.;
    auto tot = 0.;
    for(int j=0; j<numIter; j++) {
        tpm.genSeeds();
        genSeeds = genSeeds + chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - t1).count();
        t1 = chrono::high_resolution_clock::now();
        tpm.genSubGrids();
        genSubGrids = genSubGrids + chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - t1).count();
        t1 = chrono::high_resolution_clock::now();
        tpm.classiftBods();
        classiftBods = classiftBods + chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - t1).count();
        t1 = chrono::high_resolution_clock::now();
        tpm.runTrees();
        runTrees = runTrees + chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - t1).count();
        t1 = chrono::high_resolution_clock::now();
        prog.update(j);
        update = update + chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - t1).count();
        t1 = chrono::high_resolution_clock::now();
    }
    tot = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - tT).count();
    cout << endl;
    cout << "genSeeds: " << genSeeds/tot << endl;
    cout << "genSubGrids: " << genSubGrids/tot << endl;
    cout << "classiftBods: " << classiftBods/tot << endl;
    cout << "runTrees: " << runTrees/tot << endl;
    cout << "update: " << update/tot << endl;
    return bodies;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

vector<body> runTPM(vector<body>& bodies, double spacing, double width,
                    double density, int numIter, double dt){
    tree_PM tpm = tree_PM(bodies, spacing, width, density, dt);

    progbar prog = progbar(numIter, 20);
    for(int j=0; j<numIter; j++) {
        tpm.genSeeds();
        tpm.genSubGrids();
        tpm.classiftBods();
        tpm.runTrees();
        prog.update(j);
    }
    return *tpm.bodies;
}

void resetTPM(vector<body>& bodies){
    for (auto & i : bodies){
        vector<double> tmpPos = i.pos.back();
        i.pos.clear();
        i.pos.emplace_back(tmpPos);

        vector<double> tmpVel = i.vel.back();
        i.vel.clear();
        i.vel.emplace_back(tmpVel);

        vector<double> tmpAcc = i.acc.back();
        i.acc.clear();
        i.acc.emplace_back(tmpAcc);
    }
}


PYBIND11_MODULE(treecode, m) {
    //m.def("test", &test);
    m.def("basicRun", &basicRun);
    m.def("fixedBoundary", &fixedBoundary);
    m.def("particleMesh", &particleMesh);
    m.def("PMTest", &PMTest);
    m.def("PMTest1", &PMTest1);
    m.def("TreePareticleMesh", &TreePareticleMesh);
    m.def("PMTestPot", &PMTestPot);
    m.def("PMTestForce", &PMTestForce);

    m.def("runTPM", &runTPM);
    m.def("resetTPM", &resetTPM);

    py::class_<barnesHut>(m, "barnesHut")
            .def(py::init<vector<body>&, vector<double>&, vector<double>&>());
    py::class_<body>(m, "body")
            .def(py::init<>())
            .def(py::init<double&, vector<double>&,
                    vector<double>&, vector<double>&>())

            .def_property("acc", &body::getAcc, &body::setAcc)
            .def_property("vel", &body::getVel, &body::setVel)
            .def_property("mass", &body::getMass, &body::setMass)
            .def_property("soft", &body::getSoftening, &body::setSoftening)
            .def_property("pos", &body::getPos, &body::setPos);
    py::class_<tree_PM>(m, "tree_PM")
            .def(py::init<>());

    py::class_<grid>(m, "grid", py::buffer_protocol())
            .def("getF", &grid::getF)
            .def_buffer([](grid &m) -> py::buffer_info {
                return py::buffer_info(
                        m.data(),                               /* Pointer to buffer */
                        sizeof(double),                          /* Size of one scalar */
                        py::format_descriptor<double>::format(), /* Python struct-style format descriptor */
                        3,                                      /* Number of dimensions */
                        { m.size()[0], m.size()[1], m.size()[2] },                 /* Buffer dimensions */
                        { sizeof(double) * m.size()[2]*m.size()[1], sizeof(double) * m.size()[2],             /* Strides (in bytes) for each index */
                          sizeof(double) });
            });
}