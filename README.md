# MEX Without Tears!
mexcpp provides C++ classes to wrap MATLAB objects. The goal is to
make writing MEX files less painful.

This is the only library I know of that wraps cell and structure
arrays, N-dimensional arrays, and sparse matrices (forthcoming),
which are significantly more painful to deal with than plain numeric
vectors.  If you know about [Rcpp] [1] or [nr3matlab.h] [2], this
library is similar, but covers MATLAB's API more completely.

The classes are fully templatized, do not copy memory, and do not
check bounds unless requested (forthcoming). Therefore, with
optimization enabled, the wrapper should add no runtime penalty.

Typechecking (mxGetClassID) is enabled by default, but can be
disabled at your choice/risk if you need the performance.

# Requirements
This project is in its infancy, but since the basic functionality
works, I wanted to share it as soon as possible. You will need a
compiler which supports a reasonable amount of C++11.  I developed
and tested this package on Mac OS X 10.8 with Apple clang version
4.1 and MATALB R2011b. Let me know if you run into trouble on your
platform.

# Usage
To use, just `#include "mexcpp.h"`.  For a full list of examples,
see mexcpp_test.cpp. Some examples are:

Let `pRhs` be the input pointer to mexFunction. Then:

 - To wrap a scalar, `double d = scalar<double>(pRhs[i])`
 - To wrap a matrix, `Mat<double> m(pRhs[i])`
    - To address by linear index, `double x = m[3]`
    - To address by subscript, `double y = m(1, 2)`
 - Homogeneous cell arrays are currently supported. If we have a
   cell array of double matrices (all of which could have different
   dimensions), `CellMat<Mat<double> > cm(pRhs[i])`
    - To get the linear entry i, `auto entry = cm.get(i)`
      and entry has type `Mat<double>`
    - The cell array template can take any type, so you can have
      cell arrays of cell arrays of structure arrays with
      `CellMat<CellMat<StructMat> > > nestedCM`
 - Structure arrays. Interface to be finalized; see mexcpp_test.cpp.

[1]: http://dirk.eddelbuettel.com/code/rcpp.html
[2]: http://www.nr.com/nr3_matlab.html

