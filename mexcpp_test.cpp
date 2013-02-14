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
  mexPrintf("Double vector read; N = %d, M = %d\n", dv.N, dv.M);

  Mat<float> m(pIn[iSingleMatrix]);
  mexPrintf("Single matrix read; N = %d, M = %d\n", m.N, m.M);

  CellMat<Mat<double> > cm(pIn[iCellMat]);
  mexPrintf("Cell matrix read; N = %d, M = %d\n", cm.N, cm.M);

  for (int i = 0; i < cm.numEl; i++) {
    auto entry = cm.get(i);
    mexPrintf("Entry %d was %d x %d matrix\n", i, entry.N, entry.M);
  }

  StructMat sm(pIn[iStructMat]);
  mexPrintf("Structure matrix read; N = %d, M = %d, nFields = %d\n", sm.N, sm.M, sm.nFields);
  for (int i = 0; i < sm.numEl; i++) {
    auto e = sm.entry(i);
    Mat<double> f0 = e.getField<Mat<double> >("f0");
    double  f1 = e.getScalarField<double>("f1");
    int32_t f2 = e.getScalarField<int32_t>("f2");
    mexPrintf("Entry %d had f0 = [%d x %d] f1 = %g f2 = %d\n", i, f0.N, f0.M, f1, f2);
  }

  mexPrintf("All done!");

}
