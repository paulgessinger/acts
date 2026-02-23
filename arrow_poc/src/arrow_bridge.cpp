#include <random>

#include <arrow/api.h>
#include <arrow/c/bridge.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

// Build a small RecordBatch in C++ and export it via the Arrow C Data
// Interface. Returns (array_ptr, schema_ptr) as Python integers so Python can
// call:
//
//   batch = pa.RecordBatch._import_from_c(array_ptr, schema_ptr)
//
// Memory ownership note
// ─────────────────────
// ArrowArray / ArrowSchema are stack-allocated. Arrow fills in a `release`
// callback that keeps the underlying *buffer* memory alive independently of
// the local `batch` shared_ptr; PyArrow calls that callback when the returned
// RecordBatch is garbage-collected. Because _import_from_c is called
// synchronously here and sets release=nullptr once it has taken ownership, the
// structs are already consumed before the stack frame unwinds — no leak.
py::object make_record_batch(int num_tracks) {
  auto throw_if_not_ok = [](const arrow::Status& s) {
    if (!s.ok())
      throw std::runtime_error(s.ToString());
  };

  // ── Build columns ──────────────────────────────────────────────────────
  // Refactored: generate configurable number of tracks with random values
  // Use C++ <random> library for RNG

  std::vector<int64_t> id_vals(num_tracks);
  std::vector<double> pt_vals(num_tracks);
  std::vector<double> eta_vals(num_tracks);

  std::mt19937_64 rng{std::random_device{}()};
  std::uniform_real_distribution<double> pt_dist(
      0.1, 5.0);  // "pt" between 0.1 and 5.0
  std::uniform_real_distribution<double> eta_dist(
      -3.0, 3.0);  // "eta" between -3.0 and 3.0

  for (int i = 0; i < num_tracks; ++i) {
    id_vals[i] = i + 1;
    pt_vals[i] = pt_dist(rng);
    eta_vals[i] = eta_dist(rng);
  }

  arrow::Int64Builder id_builder;
  throw_if_not_ok(id_builder.AppendValues(id_vals));
  std::shared_ptr<arrow::Array> ids;
  throw_if_not_ok(id_builder.Finish(&ids));

  arrow::DoubleBuilder pt_builder;
  throw_if_not_ok(pt_builder.AppendValues(pt_vals));
  std::shared_ptr<arrow::Array> pts;
  throw_if_not_ok(pt_builder.Finish(&pts));

  arrow::DoubleBuilder eta_builder;
  throw_if_not_ok(eta_builder.AppendValues(eta_vals));
  std::shared_ptr<arrow::Array> etas;
  throw_if_not_ok(eta_builder.Finish(&etas));

  // ── Assemble RecordBatch ───────────────────────────────────────────────
  auto schema = arrow::schema({
      arrow::field("track_id", arrow::int64()),
      arrow::field("pt", arrow::float64()),
      arrow::field("eta", arrow::float64()),
  });
  auto batch = arrow::RecordBatch::Make(schema, num_tracks, {ids, pts, etas});

  // ── Export via C Data Interface ────────────────────────────────────────
  // Stack-allocate: _import_from_c is called synchronously below, so the
  // structs only need to survive until it returns. Arrow marks them as
  // consumed by setting release=nullptr, so the stack frame cleanup is safe.
  ArrowArray c_array{};
  ArrowSchema c_schema{};

  auto status = arrow::ExportRecordBatch(*batch, &c_array, &c_schema);
  if (!status.ok())
    throw std::runtime_error("ExportRecordBatch failed: " + status.ToString());

  auto py_module = py::module_::import("pyarrow");
  auto rtbatch =
      py_module.attr("RecordBatch")
          .attr("_import_from_c")(reinterpret_cast<uintptr_t>(&c_array),
                                  reinterpret_cast<uintptr_t>(&c_schema));

  return rtbatch;
}

PYBIND11_MODULE(arrow_bridge, m) {
  m.doc() = "Arrow C Data Interface zero-copy proof-of-concept";
  m.def("make_record_batch", &make_record_batch, py::arg("num_tracks") = 10,
        "Build a 5-row RecordBatch (track_id, pt, eta) in C++ and return\n"
        "(array_ptr, schema_ptr) for zero-copy import via\n"
        "pyarrow.RecordBatch._import_from_c().");
}
