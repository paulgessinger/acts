#include <iostream>

#include <podio/ROOTReader.h>

int main(int argc, char** argv) {
  podio::ROOTReader reader;
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <input file>" << std::endl;
    return 1;
  }
  std::filesystem::path file = argv[1];
  reader.openFile(file);

  auto frame = reader.readNextEntry("events");
  while (frame != nullptr) {
    std::cout << "EVENT" << std::endl;
    frame = reader.readNextEntry("events");
  }
}
