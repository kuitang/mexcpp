/* mexcpp_test.cpp -- Test mexcpp library
 *
 * Copyright (c) 2013-5 Kui Tang <kuitang@gmail.com>
 * License: MIT
 *
 * Inspired by Rcpp and nr3matlab.h
 *
 * Run this line in MATLAB:
 * [od, os, ocm, osm] = mexcpp_test(3, 'c', 'string', [1 2 3], single([1 2; 3 4]), {[1 2], [1 2; 3 4]}, struct('f0', [10 20; 30 40], 'f1', 52, 'f2', int32(99), 'f3', 'blahblah'))
 */

#include "mexcpp.h"
#include <string>

#ifdef HAVE_EIGEN
#include <Eigen/Dense>
#include <sstream>
#endif

enum {
  iDouble,
  iChar,
  iStr,
  iDoubleVector,
  iSingleMatrix,
  iCellMat,
  iStructMat,
  nI
};

enum {
  oDoubleMatrix,
  oString,
  oCellMat,
  oStructMat,
  nO
};

using namespace mexcpp;

// test pointer overloading
void printPointerArray(int n, double *p) {
  for (int i = 0; i < n; i++) {
    mexPrintf("printPointerArray: p[%d] = %g\n", i, p[i]);
  }
}

void mexFunction(int nOut, mxArray *pOut[], int nIn, const mxArray *pIn[]) {
  if (nIn != nI || nOut != nO) {
    mexErrMsgIdAndTxt("mexcpp_test:nIn", "Usage: mexcpp_test(double, char, doubleVec, singleMat, cellMat, structMat)", nI);
  }

  double d = scalar<double>(pIn[iDouble]);
  mexPrintf("Double scalar read; d = %g\n", d);

  char c = scalar<char>(pIn[iChar]);
  mexPrintf("Char scalar read; c = %c\n", c);

  std::string s = scalar<std::string>(pIn[iStr]);
  mexPrintf("String scalar read; s = %s\n", s.c_str());

  Mat<double> dv(pIn[iDoubleVector]);
  mexPrintf("Double vector read; N = %d, M = %d, end = %g\n", dv.N, dv.M, dv[dv.length - 1]);

  printPointerArray(dv.length, dv);

  Mat<float> m(pIn[iSingleMatrix]);
  mexPrintf("Single matrix read; N = %d, M = %d, end = %g\n", m.N, m.M, m(m.N - 1, m.M - 1));

#ifdef HAVE_EIGEN
  std::stringstream ss;
  const Mat<float>::EigenMap mEigenMap(m.asEigenMap());
  ss << "mEigenMap contents: \n" << mEigenMap << "\n";
  mexPrintf(ss.str().c_str());
  ss.str("");

  Mat<float>::EigenMatrix mEigenMat(m.asEigenMatrix());
  ss << "mEigenMatrix contents: \n" << mEigenMap << "\n";
  mexPrintf(ss.str().c_str());
  ss.str() = "";

#endif

  CellMat<Mat<double> > cm(pIn[iCellMat]);
  mexPrintf("Cell matrix read; N = %d, M = %d\n", cm.N, cm.M);

  for (int i = 0; i < cm.length; i++) {
    mexPrintf("Entry %d was %d x %d matrix\n", i, cm[i].N, cm[i].M);
  }

  StructMat sm(pIn[iStructMat]);
  mexPrintf("Structure matrix read; N = %d, M = %d, nFields = %d\n", sm.N, sm.M, sm.nFields);
  for (int i = 0; i < sm.length; i++) {
    Mat<double> f0 = sm[i].get<Mat<double> >("f0");
    double      f1 = sm[i].getS<double>("f1");
    int32_t     f2 = sm[i].getS<int32_t>("f2");
    std::string f3 = sm[i].getS<std::string>("f3");
    mexPrintf("Entry %d had f0 = [%d x %d] f1 = %g f2 = %d f3 = %s\n", i, f0.N, f0.M, f1, f2, f3.c_str());
  }

  pOut[oString] = scalar<std::string>("String output");

  // Making stuff
  Mat<double> om(2,2);
  om(0,0) = 1;
  om(0,1) = 2;
  om(1,0) = 3;
  om(1,1) = 4;
  pOut[oDoubleMatrix] = om;
  mexPrintf("Double matrix output created.\n");

#ifdef HAVE_EIGEN
  Mat<double>::EigenMap omEigenMap(om.asEigenMap());
  omEigenMap(0,0) += 0.1;
  mexPrintf("Created (nonconst) EigenMap on output matrix and manipulated it.");
#endif

  CellMat<Mat<double> > ocm(1,2);
  ocm.setS(0, 5.1);
  ocm.set(1, Mat<double>(3, 3));
  pOut[oCellMat] = ocm;
  mexPrintf("Cell matrix output created.\n");
  pOut[oCellMat] = ocm;

  // No C++11 features
  std::vector<std::string> fns;
  fns.push_back("bar");
  fns.push_back("quuz");
  StructMat osm(2, 1, fns);
  osm[0].set("bar", Mat<double>(3,3));
  osm[0].setS("quuz", 1);
  osm[1].set(0, Mat<double>(2,2));
  osm[1].setS(1, 2);
  pOut[oStructMat] = osm;
  mexPrintf("Struct matrix output created.\n");

  mexPrintf("All done!\n");
}
