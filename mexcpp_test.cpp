// Run this line in MATLAB:
// mexcpp_test(3, 'c', [1 2 3], single([1 2; 3 4]), {[1 2], [1 2; 3 4]}, struct('f0', [10 20; 30 40], 'f1', 52, 'f2', int32(99), 'f3', 'blahblah'))
#include "mexcpp.h"

enum {
  iDouble,
  iChar,
  iDoubleVector,
  iSingleMatrix,
  iCellMat,
  iStructMat,
  nI
};

using namespace mexcpp;

void mexFunction(int nOut, mxArray *pOut[], int nIn, const mxArray *pIn[]) {
  if (nIn != nI) {
    mexErrMsgIdAndTxt("mexcpp_test:nIn", "Usage: mexcpp_test(double, char, doubleVec, singleMat, cellMat, structMat)", nI);
  }

  double d = scalar<double>(pIn[iDouble]);
  mexPrintf("Double scalar read; d = %g\n", d);

  char c   = scalar<char>(pIn[iChar]);
  mexPrintf("Char scalar read; c = %c\n", c);

  Mat<double> dv(pIn[iDoubleVector]);
  mexPrintf("Double vector read; N = %d, M = %d, end = %g\n", dv.N, dv.M, dv[dv.numEl - 1]);

  Mat<float> m(pIn[iSingleMatrix]);
  mexPrintf("Single matrix read; N = %d, M = %d, end = %g\n", m.N, m.M, m(m.N - 1, m.M - 1));

  CellMat<Mat<double> > cm(pIn[iCellMat]);
  mexPrintf("Cell matrix read; N = %d, M = %d\n", cm.N, cm.M);

  for (int i = 0; i < cm.numEl; i++) {
    mexPrintf("Entry %d was %d x %d matrix\n", i, cm[i].N, cm[i].M);
  }

  StructMat sm(pIn[iStructMat]);
  mexPrintf("Structure matrix read; N = %d, M = %d, nFields = %d\n", sm.N, sm.M, sm.nFields);
  for (int i = 0; i < sm.numEl; i++) {
    auto f0 = sm[i].f<Mat<double> >("f0");
    auto f1 = sm[i].sf<double>("f1");
    auto f2 = sm[i].sf<int32_t>("f2");
    auto f3 = str(sm[i].field("f3"));
    mexPrintf("Entry %d had f0 = [%d x %d] f1 = %g f2 = %d f3 = %s\n", i, f0.N, f0.M, f1, f2, f3.c_str());
  }

  mexPrintf("All done!\n");

}
