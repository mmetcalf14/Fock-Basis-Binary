#include <iostream>
#include <vector>
#include "search.h"

int main() {
  std::vector<double> values = {0, 0.5, 1, 5, 7.5, 10, 12.5};
  std::vector<double> tests = {0, 0.4, 0.5, 3, 7.5, 11.5, 12.5, 13};
  for(double d : tests) {
    auto it = binary_locate(values.begin(), values.end(), d);
    std::cout << "found " << d << " right after index " << std::distance(values.begin(), it) << " which has value " << *it << std::endl;
  }
  return 0;
}
