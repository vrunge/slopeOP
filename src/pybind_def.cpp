// Copyright (c) 2019 Marco Pascucci

#ifdef BUILD_PYTHON_MODULE

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "op2d.h"
#include "Omega.h"
#include <string>

namespace py = pybind11;

PeltResult<float, float> op2D_float(vector<float> &x, vector<float> &y, float beta)
{
    return op2D(x, y, beta);
}

Omega *slopeOP(std::vector<double> data, std::vector<double> states, double penalty,
               std::string constraint = "null", double minAngle = 0, std::string type = "channel")
{
    Omega *omega = new Omega(states, penalty, data.size());
    //DIFFERENT PRUNING
    if (type == "null" && constraint == "null")
    {
        omega->algo(data);
    }
    if (type == "channel" && constraint == "null")
    {
        omega->algoChannel(data);
    }
    if (type == "pruning" && constraint == "null")
    {
        omega->algoPruning(data);
    }
    if (type == "pruningMax" && constraint == "null")
    {
        omega->algoPruningMyList(data);
    }

    //DIFFERENT CONSTRAINTS
    if (constraint == "isotonic")
    {
        omega->algoISOTONIC(data);
    }
    if (constraint == "unimodal")
    {
        omega->algoUNIMODAL(data);
    }
    if (constraint == "smoothing")
    {
        omega->algoSMOOTHING(data, minAngle);
    }

    omega->backtracking(data.size());

    return omega;
}

PYBIND11_MODULE(slopeOP, m)
{
    py::class_<PeltResult<float, float>>(m, "PeltResult")
        .def(py::init([](vector<unsigned int> cp, vector<float> x, vector<float> y, double cost) {
            return new PeltResult<float, float>(cp, x, y, cost);
        }))
        .def_readwrite("cp", &PeltResult<float, float>::cp)
        .def_readwrite("x", &PeltResult<float, float>::x)
        .def_readwrite("y", &PeltResult<float, float>::y)
        .def_readwrite("cost", &PeltResult<float, float>::cost);

    py::class_<Omega>(m, "Omega")
        .def("GetChangepoints", &Omega::GetChangepoints, pybind11::return_value_policy::copy)
        .def("GetParameters", &Omega::GetParameters)
        .def("GetGlobalCost", &Omega::GetGlobalCost)
        .def("GetPruning", &Omega::GetPruning);

    m.doc() = "python interface to segmentaition algorithms"; // optional module docstring
    m.def("op2D", &op2D_float, "OP algorithm 2D (piece-wise linear fit)", py::arg("x"), py::arg("y"), py::arg("penality"));
    m.def("slopeOP", &slopeOP, "slopeOP algorithm", py::arg("data"), py::arg("states"), py::arg("penality"),
          py::arg("constraint") = "null", py::arg("minAngle") = 0, py::arg("type") = "channel");
}
#endif
