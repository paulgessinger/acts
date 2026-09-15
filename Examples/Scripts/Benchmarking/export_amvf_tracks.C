// Export one event of fitted track parameters from a track-summary ROOT file.

#include <array>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>

void export_amvf_tracks(const char* input, const char* output,
                        long long event = 0) {
  TFile file(input);
  auto* tree = file.Get<TTree>("tracksummary");
  if (tree == nullptr || event < 0 || event >= tree->GetEntries()) {
    throw std::runtime_error("Requested track-summary event does not exist");
  }

  const std::array<std::string, 12> names = {
      "eLOC0_fit",    "eLOC1_fit",      "ePHI_fit",      "eTHETA_fit",
      "eQOP_fit",     "eT_fit",         "err_eLOC0_fit", "err_eLOC1_fit",
      "err_ePHI_fit", "err_eTHETA_fit", "err_eQOP_fit",  "err_eT_fit"};
  std::array<std::vector<float>*, 12> values{};
  for (std::size_t i = 0; i < names.size(); ++i) {
    tree->SetBranchAddress(names[i].c_str(), &values[i]);
  }
  tree->GetEntry(event);

  std::ofstream stream(output);
  if (!stream) {
    throw std::runtime_error("Could not create output CSV");
  }
  stream << "loc0,loc1,phi,theta,qop,time,sigma_loc0,sigma_loc1,sigma_phi,"
            "sigma_theta,sigma_qop,sigma_time\n";
  for (std::size_t row = 0; row < values[0]->size(); ++row) {
    for (std::size_t column = 0; column < values.size(); ++column) {
      stream << (column == 0 ? "" : ",") << values[column]->at(row);
    }
    stream << '\n';
  }
}
