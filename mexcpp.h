/* mexcpp.h -- Seamless MATLAB and C++ Integration
 *
 * This file wraps MATLAB datatypes (arrays, classures, cell arrays) in C++
 * template classes. Numeric arrays can be further wrapped into Eigen.
 *
 * This version does not support ownership: the classes are intended to be used
 * for arguments to mexFunction only. Ownership support (so we can use Engine)
 * is planned.
 *
 * TODO: Support complex!
 *
 * Version 0.2 -- compiles on 4.6.3, no C++11 needed. scalar supports string.
 * Copyright (c) 2013 Kui Tang <kuitang@gmail.com>
 *
 * Inspired by Rcpp and nr3matlab.h
 */

#pragma once
#include <algorithm>
#include <cmath>
#include <stdint.h>
#include <complex>
#include <exception>
#include <vector>
#include <string>

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
  class mxCell;
  class mxStruct;
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

  // Extract a scalar from a MATLAB pointer
  template<class T>
  T scalar(const mxArray *prhs) {
    checkTypeOrErr<T>(prhs);
    return *(static_cast<T*>(mxGetData(prhs)));
  }

  template<>
  std::string scalar(const mxArray *prhs) {
    return std::string(mxArrayToString(prhs));
  }

  // Set a MATLAB pointer from a scalar
  template<class T>
  mxArray *scalar(const T &x) {
    mxArray *pm = mxCreateNumericMatrix(1, 1, ty<T>(), mxREAL);
    *(static_cast<T*>(mxGetData(pm))) = x;
    return pm;
  }

  template<>
  mxArray *scalar(const std::string &s) {
    mxArray *pm = mxCreateString(s.c_str());
    return pm;
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
    // M is rows, N is cols
    size_t N, M, length;

    BaseMat(const mxArray *p) : pm(const_cast<mxArray *>(p)) {
      //mexPrintf("BaseMat: p = %lx\n", p);
      M = mxGetM(p);
      N = mxGetN(p);
      length = N * M;
    }

    BaseMat(size_t M_, size_t N_) : pm(0), M(M_), N(N_), length(N_ * M_) { }

    // Conversion cast operator
    operator mxArray *() { return pm; }

    inline size_t sub2ind(size_t r, size_t c) { return c*M + r; }

  protected:
    mxArray *pm;
  };

  // TODO: Specialize for a complex!
  template<class T>
  struct Mat : public BaseMat {
    int complexity;
    T *re;

    // Construct from MATLAB pointer
    Mat(const mxArray *p) : BaseMat(p) {
      checkTypeOrErr<T>(p);
      re = static_cast<T *>(mxGetData(p));
    }

    // Construct our own (in MATLAB's memory)
    Mat(size_t rows, size_t cols) : BaseMat(rows, cols) {
      pm = mxCreateNumericMatrix(rows, cols, ty<T>(), mxREAL);
      re = static_cast<T *>(mxGetData(pm));
    }

    // Construct from existing array (copy to MATLAB's memory)
    Mat(size_t rows, size_t cols, T *p) : BaseMat(rows, cols) {
      pm = mxCreateNumericMatrix(rows, cols, ty<T>(), mxREAL);
      re = static_cast<T *>(mxGetData(pm));
      std::copy(p, p + rows*cols, re);
    }

    // Construct from existing std::vector (copy to MATLAB's memory)
    Mat(const std::vector<T> v, bool colMat=true) : BaseMat(0, 0) {
      length = v.size();
      if (colMat) {
        M = length;
        N = 1;
      } else {
        M = 1;
        N = length;
      }

      pm = mxCreateNumericMatrix(M, N, ty<T>(), mxREAL);
      re = static_cast<T *>(mxGetData(pm));
      std::copy(v.data(), v.data() + v.size(), re);
    }


    // Pointer cast operator. Only defined for numeric arrays.
    // Allows usage in cases where e.g. double *vec is expected
    operator T*() { return re; }

    // NO RANGE CHECKING!
    T *col(size_t c) { return re + sub2ind(0, c); }
    T &operator[](size_t ind) { return re[ind]; }
    T &operator()(size_t r, size_t c) { return re[sub2ind(r,c)]; }
  };

  // Homogeneous cell array (Mat, CellMat, classMat, etc.)
  template<class CellT>
  struct CellMat : public BaseMat {
    CellMat(const mxArray *p) : BaseMat(p) { checkTypeOrErr<mxCell>(p); }

    // Construct our own (in MATLAB's memory)
    CellMat(size_t rows, size_t cols) : BaseMat(rows, cols) {
      pm = mxCreateCellMatrix(rows, cols);
    }

    template<class T>
    void set(size_t ind, const T &x) { mxSetCell(pm, ind, const_cast<T &>(x)); }
    template<class S>
    void setS(size_t ind, const S &x) { mxSetCell(pm, ind, scalar(x)); }

    template<class T>
    void set(size_t r, size_t c, const T &x) { set(sub2ind(r, c), x); }
    template<class S>
    void setS(size_t r, size_t c, const S &x) { setS(sub2ind(r, c), x); }

    mxArray *ptr(size_t ind) { return mxGetCell(pm, ind); }
    mxArray *ptr(size_t r, size_t c) { return ptr(sub2ind(r, c)); }

    CellT operator[](size_t ind) { return CellT(ptr(ind)); }
    CellT operator()(size_t r, size_t c) { return CellT(ptr(sub2ind(r, c))); }
  };

  // For use in classMat and classNDArray. Just stores a pointer to
  // mxArray (a struct array) and in index.
  struct Entry {
    mxArray *pmm;
    mwIndex ind;

    Entry(mxArray *pmm_, mwIndex ind_) : pmm(pmm_), ind(ind_) {
      //mexPrintf("Entry constructed; pmm = %lx ,ind = %d\n", pmm, ind);
    }

    mxArray *field(mwIndex fnum) {
      //mexPrintf("field: pmm = %x, ind = %d, fnum = %d\n", pmm, ind, fnum);
      mxArray *f = mxGetFieldByNumber(pmm, ind, fnum);
      if (f == 0) {
        mexErrMsgIdAndTxt("mexcpp:field", "Field %d of index %d of matrix %lx was invalid.\n", fnum, ind, pmm);
      }
      return f;
    }

    mxArray *field(const char *fname) {
      //mexPrintf("field: pmm = %x, ind = %d, fname = %s\n", pmm, ind, fname);
      mxArray *f = mxGetField(pmm, ind, fname);
      if (f == 0) {
        mexErrMsgIdAndTxt("mexcpp:field", "Field %s of index %d of matrix %lx was invalid.\n", fname, ind, pmm);
      }
      return f;
    }

    // Get an instantiated field
    template <class T, class F> T get(F fn) { return T(field(fn)); }
    template <class S, class F> S getS(F fn) { return scalar<S>(field(fn)); }

    template <class T> void set(const char *fn, T x) {
      mxSetField(pmm, ind, fn, x);
    }

    template <class S> void setS(const char *fn, S x) {
      mxSetField(pmm, ind, fn, scalar<S>(x));
    }

    template <class T> void set(mwIndex fi, T x) {
      mxSetFieldByNumber(pmm, ind, fi, x);
    }

    template <class S> void setS(mwIndex fi, S x) {
      mxSetFieldByNumber(pmm, ind, fi, scalar<S>(x));
    }

  };

  struct StructMat : public BaseMat {
    size_t nFields;

    StructMat(const mxArray *p) : BaseMat(p) {
      checkTypeOrErr<mxStruct>(p);
      nFields = mxGetNumberOfFields(p);
    }

    StructMat(size_t rows, size_t cols, std::vector<std::string> fieldNames) :
      BaseMat(rows, cols) {
      nFields = fieldNames.size();
      const char *fns[nFields];

      for (int i = 0; i < fieldNames.size(); i++) {
        fns[i] = fieldNames[i].c_str();
      }

      pm = mxCreateStructMatrix(rows, cols, nFields, fns);
    }

    Entry operator[](size_t ind) { return Entry(pm, ind); }
    Entry operator()(size_t r, size_t c) { return (*this)[(sub2ind(r,c))]; }

    // Syntactic sugar; access the first element.
    // Deals with the common use case with a 1x1 "matrix".
    template <class T> void set(const char *fn, T x) {
      (*this)[0].set(fn, x);
    }

    template <class S> void setS(const char *fn, S x) {
      (*this)[0].setS(fn, x);
    }

    template <class T> void set(mwIndex fi, T x) {
      (*this)[0].set(fi, x);
    }

    template <class S> void setS(mwIndex fi, S x) {
      (*this)[0].set(fi, x);
    }

    template <class T, class F> T get(F fn) { return (*this)[0].get<T>(fn); }
    template <class S, class F> S getS(F fn) { return (*this)[0].getS<S>(fn); }
  };

  // Compressed sparse column sparse matrix (MATLAB's format)
  struct SparseMat : public BaseMat {
    void construct(const mxArray *pm) {
      if (!mxIsDouble(pm) || !mxIsSparse(pm)) {
        mexErrMsgIdAndTxt("SparseMat:type", "matrix pm must be double and sparse");
      }

      N = mxGetN(pm);
      M = mxGetM(pm);
      nzMax = mxGetNzmax(pm);
      pr = mxGetPr(pm);
      ir = mxGetIr(pm);
      jc = mxGetJc(pm);
    }

    size_t nzMax;
    double *pr;
    mwIndex *ir, *jc;

    // Grab from MATLAB
    SparseMat(mxArray *pm) : BaseMat(pm) { construct(pm); }

    // Construct our own
    SparseMat(size_t rows, size_t cols, size_t nnz) :
      BaseMat(rows, cols) {
      pm = mxCreateSparse(rows, cols, nnz, mxREAL);
      construct(pm);
    }
  };

  /*
  template<class T>
  class NDArray {
    mxArray *pm;
    int complexity;
    size_t nDims, length;
    size_t *dims;
    T *re;

    NDArray(const mxArray *p) : pm(p) {
      re = mxGetData(p);
      nDims = mxGetNumberOfDimensions(pm);
      dims  = mxGetDimensions(p);
      length = mxGetNumberOfElements(pm);
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
