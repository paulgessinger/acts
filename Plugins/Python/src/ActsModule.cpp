#include <pybind11/pybind11.h>

namespace py = pybind11;

int add(int i, int j) {
  return i + j + 1;
}

PYBIND11_MODULE(_acts, m) {
  m.doc() = "Acts";
  m.def("add", &add, "Add stuff");
}