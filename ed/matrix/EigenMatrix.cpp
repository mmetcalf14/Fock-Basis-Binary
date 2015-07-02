#include <iostream>
#include <vector>
#include "EDType.h"
#include "matrix/EigenMatrix.h"

/*NOTE: This is an example to initialize sparse matrix.*/

int main(int argc, char const *argv[]) {
  size_t rows = 10, cols = 10;
  std::vector<RealTriplet> tripletList;
  tripletList.reserve(20);
  for (size_t cnt = 0; cnt < 9; cnt++) {
    tripletList.push_back(RealTriplet(cnt, cnt+1, -1.0));
    tripletList.push_back(RealTriplet(cnt+1, cnt, -1.0));
    // tripletList.push_back(RealTriplet(cnt+1, cnt, -2.0));
  }
  RealSparseMatrixType mat(rows,cols);
  mat.setFromTriplets(tripletList.begin(), tripletList.end());
  for (size_t k=0; k < mat.outerSize(); ++k){
    for (RealSparseMatrixType::InnerIterator it(mat,k); it; ++it)
    {
      INFO(it.index()
          << " - ( " << it.row() << ", " << it.col() << " ) - "
          << it.value() );
    }
  }

  // RealMatrixType dense = RealMatrixType(mat);
  INFO(mat);
  return 0;
}
