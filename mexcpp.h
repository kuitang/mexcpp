/* mexcpp.h -- Seamless MATLAB and C++ Integration
 *
 * This file wraps MATLAB datatypes (arrays, structures, cell arrays) in C++
 * template classes. Numeric arrays can be further wrapped into Eigen.
 *
 * This version does not support ownership: the classes are intended to be used
 * for arguments to mexFunction only. Ownership support (so we can use Engine)
 * is planned.
 *
 * TODO: Support complex!
 *
 * Version 0.1
 * Copyright (c) 2013 Kui Tang <kuitang@gmail.com>
 *
 * Inspired by Rcpp nr3matlab.h
 */

#pragma once
#include <algorithm>
#include <cmath>
#include <stdint.h>
#include <complex>
#include <exception>

#include "mex.h"
#include "matrix.h"

namespace mexcpp {

  // Template nonsense for MATLAB types
  template<class T> inline mxClassID ty() {return mxUNKNOWN_CLASS;}

  // Numeric types
  template<> inline mxClassID ty<double>() {return mxDOUBLE_CLASS;}
  template<> inline mxClassID ty<float>() {return mxSINGLE_CLASS;}
  template<> inline mxClassID ty<int32_t>() {return mxINT32_CLASS;}
  template<> inline mxClassID ty<uint32_t>() {return mxUINT32_CLASS;}
  template<> inline mxClassID ty<char>() {return mxCHAR_CLASS;}
  template<> inline mxClassID ty<uint8_t>() {return mxUINT8_CLASS;}
  template<> inline mxClassID ty<int64_t>() {return mxINT64_CLASS;}
  template<> inline mxClassID ty<uint64_t>() {return mxUINT64_CLASS;}
  template<> inline mxClassID ty<bool>() {
    if (sizeof(bool)==1) return mxLOGICAL_CLASS;
    else throw("bool and mxLOGICAL_CLASS have incompatible sizes");
  }

  // Symbolic types
  struct mxCell;
  struct mxStruct;
  template<> inline mxClassID ty<mxCell>() {return mxCELL_CLASS;}
  template<> inline mxClassID ty<mxStruct>() {return mxSTRUCT_CLASS;}

  inline mxClassID ty(const mxArray *p) {return mxGetClassID(p);}

#ifdef NTYPECHECK
  template<class T>
  inline void checkTypeOrErr(const mxArray *m) { }
#else
  template<class T>
  inline void checkTypeOrErr(const mxArray *m) {
    if (mxGetClassID(m) != ty<T>()) {
      mexErrMsgTxt("Wrong type");
    }
  }
#endif

  // Scalar library
  // just a placeholder
  template<class T>
  struct Scalar { };

  template<class T>
  const T &scalar(const mxArray *prhs) {
    checkTypeOrErr<T>(prhs);
    return *(static_cast<T*>(mxGetData(prhs)));
  }

  // This function copies memory. It is expected to be used for small
  // user-supplied strings.
  std::string str(const mxArray *prhs) {
    checkTypeOrErr<char>(prhs);
    char *s = mxArrayToString(prhs);
    std::string ret(s);
    mxFree(s);
    return ret;
  }

  // TODO: Check for memory allocation error (when running outside of
  // mexFunction)
  template<class T> mxArray *createScalar(T x) {
    mxArray *ret = mxCreateNumericMatrix(1, 1, ty<T>());
    *(static_cast<T *>(mxGetData(ret))) = x;
    return ret;
  }

  template<> mxArray *createScalar<double>(double x) {
    return mxCreateDoubleScalar(x);
  }

  template<> mxArray *createScalar<bool>(bool x) {
    return mxCreateLogicalScalar(x);
  }

  // TODO: Check the generated code
  struct BaseMat {
    const mxArray *pm;
    size_t N, M, numEl;

    BaseMat(const mxArray *p) : pm(p) {
      //mexPrintf("BaseMat: p = %lx\n", p);
      N = mxGetN(p);
      M = mxGetM(p);
      numEl = N * M;
    }

    size_t sub2ind(size_t r, size_t c) { return c*N + r; }
  };

  // TODO: Specialize for a complex!
  template<class T>
  struct Mat : public BaseMat {
    T *re;
    int complexity;

    Mat(const mxArray *p) : BaseMat(p) {
      checkTypeOrErr<T>(p);
      re = static_cast<T *>(mxGetData(p));
    }

    // NO RANGE CHECKING!
    T *col(size_t c) { return re + sub2ind(0, c); }
    T &operator[](size_t ind) { return re[ind]; }
    T &operator()(size_t r, size_t c) { return re[sub2ind(r,c)]; }
  };

  // Homogeneous cell array (Mat, CellMat, StructMat, etc.)
  template<class CellT>
  struct CellMat : public BaseMat {
    CellMat(const mxArray *p) : BaseMat(p) { checkTypeOrErr<mxCell>(p); }

    mxArray *getPtr(size_t ind) { return mxGetCell(pm, ind); }
    CellT operator[](size_t ind) { return CellT(getPtr(ind)); }
    CellT operator()(size_t r, size_t c) { return CellT(getPtr(sub2ind(r, c))); }
  };

  // For use in StructMat and StructNDArray
  struct Entry {
    const mxArray *pmm;
    size_t ind;

    Entry(const mxArray *pmm_, size_t ind_) : pmm(pmm_), ind(ind_) {
      //mexPrintf("Entry constructed; pmm = %lx ,ind = %d\n", pmm, ind);
    }

    mxArray *field(size_t fnum) {
      //mexPrintf("field: pmm = %x, ind = %d, fnum = %d\n", pmm, ind, fnum);
      mxArray *f = mxGetFieldByNumber(pmm, ind, fnum);
      if (f == 0) {
        mexErrMsgIdAndTxt("mexcpp:field", "Field %d of index %d of struct matrix %lx was invalid.\n", fnum, ind, pmm);
      }
      return f;
    }

    mxArray *field(const char *fname) {
      //mexPrintf("field: pmm = %x, ind = %d, fname = %s\n", pmm, ind, fname);
      mxArray *f = mxGetField(pmm, ind, fname);
      if (f == 0) {
        mexErrMsgIdAndTxt("mexcpp:field", "Field %s of index %d of struct matrix %lx was invalid.\n", fname, ind, pmm);
      }
      return f;
    }

    // Get an instantiated field
    template <class T, class F> T f(F fn) { return T(field(fn)); }
    template <class S, class F> S sf(F fn) { return scalar<S>(field(fn)); }
  };

  struct StructMat : public BaseMat {
    size_t nFields;

    StructMat(const mxArray *p) : BaseMat(p) {
      checkTypeOrErr<mxStruct>(p);
      nFields = mxGetNumberOfFields(p);
    }

    Entry operator[](size_t ind) { return Entry(pm, ind); }
    Entry operator()(size_t r, size_t c) { return (*this)[(sub2ind(r,c))]; }
  };

  /*
  template<class T>
  class NDArray {
    mxArray *pm;
    int complexity;
    size_t nDims, numEl;
    size_t *dims;
    T *re;

    NDArray(const mxArray *p) : pm(p) {
      re = mxGetData(p);
      nDims = mxGetNumberOfDimensions(pm);
      dims  = mxGetDimensions(p);
      numEl = mxGetNumberOfElements(pm);
    }

    T &operator[](size_t ind) {
      return re[ind];
    }

    T &operator()(size_t *subs) {
      size_t ind = mxCalcSingleSubscript(pm, nDims, subs);
      return re[ind];
    }
  };
  */

  // NDArray: Need variadic!

}
