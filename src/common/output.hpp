#pragma once
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

namespace etch {
template <class T> void write_binary(const std::string& path, const std::vector<T>& data) {
  std::ofstream file(path, std::ios::binary);
  file.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(T));
}

class Summary {
  std::ofstream file;
  bool first = true;

  void key(const std::string& name) {
    file << (first ? "\n  " : ",\n  ") << '"' << name << "\": ";
    first = false;
  }

public:
  explicit Summary(const std::string& path) : file(path) {
    file << '{' << std::setprecision(12) << std::boolalpha;
  }

  ~Summary() {
    file << "\n}\n";
  }

  template <class T> void field(const std::string& name, const T& value) {
    key(name);
    file << value;
  }

  void field(const std::string& name, const std::string& value) {
    key(name);
    file << '"' << value << '"';
  }
};
}
