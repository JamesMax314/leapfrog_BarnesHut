#include <iostream>
#include <fstream>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "bodies.h"
#include "trees.h"
#include "leapfrog.h"
#include "pyInterface.h"

namespace py = pybind11;

vector<vector<body>> add(vector<body> bodies) {
    return {bodies, bodies};
}

PYBIND11_MODULE(treecode, m) {
    m.def("add", &add);
    py::class_<body>(m, "body")
            .def(py::init<>())
            .def(py::init<double&, vector<double>&,
                    vector<double>&, vector<double>&,
                    double&, double&>())
            .def("setPos", &body::setPos)
            .def("setAcc", &body::setAcc)
            .def("setVel", &body::setVel)
            .def("setMass", &body::setMass)
            .def("getPos", &body::getPos)
            .def("getAcc", &body::getAcc)
            .def("getVel", &body::getVel)
            .def("getMass", &body::getMass);
}