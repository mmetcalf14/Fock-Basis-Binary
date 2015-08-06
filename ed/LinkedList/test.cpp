#include <iostream>
#include "LinkedList/Site.hpp"

int main(int argc, char const *argv[]) {
  Site<double, char> A('a');
  Site<double, char> B('b', &A);
  // A.LinkTo(&B, 1.0);
  // B.LinkTo(&A, 1.0);

  std::cout << A.TotalSites() << " " << A.NumNeighbors() << " " <<
    A.VerifySite() << std::endl;
  for (auto &j : A.getNeighbors() ){
    std::cout << j->data << std::endl;
  }
  std::cout << B.TotalSites() << " " << B.NumNeighbors() << " " <<
    B.VerifySite() << std::endl;
  for (auto &j : B.getNeighbors() ){
    std::cout << j->data << std::endl;
  }
  return 0;
}
