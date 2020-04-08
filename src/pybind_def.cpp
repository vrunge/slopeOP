// Copyright (c) 2019 Marco Pascucci

#ifdef BUILD_PYTHON_MODULE

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "op2d.h"
#include "OmegaOP.h"
#include "OmegaSN.h"
#include <string>
#include <vector>

namespace py = pybind11;

PeltResult<float, float> op2D_float(vector<float> &x, vector<float> &y, float beta)
{
    return op2D(x, y, beta);
}

OmegaOP *slopeOP(std::vector<double> data, std::vector<double> states, double penalty,
               std::string constraint = "null", double minAngle = 0, std::string type = "channel")
{
    OmegaOP *omega = new OmegaOP(states, data[0], penalty, data.size());

    //DIFFERENT PRUNING + NO CONSTRAINT
    if(type == "null" && constraint == "null"){omega->algo(data);}
    if(type == "channel" && constraint == "null"){omega->algoChannel(data);}
    if(type == "pruning" && constraint == "null"){omega->algoPruning(data);}

    //DIFFERENT CONSTRAINTS
    if(constraint == "isotonic"){omega->algoISOTONIC(data);}
    if(constraint == "unimodal"){omega->algoUNIMODAL(data);}
    if(constraint == "smoothing"){omega->algoSMOOTHING(data, minAngle);}

    omega->backtracking(data.size());

    return omega;
}


OmegaSN *slopeSN(std::vector<double> data, std::vector<double> states,
                    unsigned int nbSegments, std::string constraint = "null")
{
  OmegaSN *omega = new OmegaSN(states, data[0], nbSegments, data.size());
  //DIFFERENT CONSTRAINTS
  if(constraint == "null"){omega->algoNULL(data);}
  if(constraint == "isotonic"){omega->algoISOTONIC(data);}

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

    py::class_<OmegaOP>(m, "OmegaOP")
        .def("GetChangepoints", [](const OmegaOP &a) {
                //
                std::vector<int> cp = a.GetChangepoints();
                for (size_t i=0; i<cp.size(); i++) {
                    cp[i] --;
                }
                return cp;
            }, pybind11::return_value_policy::copy)
        .def("GetParameters", &OmegaOP::GetParameters)
        .def("GetGlobalCost", &OmegaOP::GetGlobalCost)
        .def("GetPruning", &OmegaOP::GetPruning);

    py::class_<OmegaSN>(m, "OmegaSN")
        .def("GetChangepoints", [](const OmegaSN &a) {
                //
                std::vector<int> cp = a.GetChangepoints();
                for (size_t i=0; i<cp.size(); i++) {
                    cp[i] --;
                }
                return cp;
            }, pybind11::return_value_policy::copy)
        .def("GetParameters", &OmegaSN::GetParameters)
        .def("GetGlobalCost", &OmegaSN::GetGlobalCost)
        .def("GetPruning", &OmegaSN::GetPruning);

    m.doc() = "python interface to segmentaition algorithms"; // optional module docstring
    m.def("op2D", &op2D_float, "OP algorithm 2D (piece-wise linear fit)", py::arg("x"), py::arg("y"), py::arg("penality"));
    m.def("slopeOP", &slopeOP, "slopeOP algorithm", py::arg("data"), py::arg("states"), py::arg("penality"),
          py::arg("constraint") = "null", py::arg("minAngle") = 0, py::arg("type") = "channel");
    m.def("slopeSN", &slopeSN, "slopeSN algorithm", py::arg("data"), py::arg("states"), py::arg("nb_segments"),
          py::arg("constraint") = "null");
}
#endif
