#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>
#include "homolig.hpp"

PYBIND11_MAKE_OPAQUE(std::vector<double>);
PYBIND11_MAKE_OPAQUE(std::vector<std::string>);

namespace py = pybind11;

PYBIND11_MODULE(homoligcpp, m){ //for some reason this works. 
    m.doc() = "homolig alignment in cpp"; // optional module docstring
    py::bind_vector<std::vector<double>>(m, "VectorDouble");
    py::bind_vector<std::vector<std::string>>(m, "VectorString"); // need to build this as an inherited class so can construct vector from strings
    m.def("homolig", &homolig, "A function that creates the amino-acid matrix and gets alignment scores of amino acids");
};
