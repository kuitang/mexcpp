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
 *     : Fix size_t and mwIndex. (mxGetN is size_t but everything else mwIndx)
 *
 * Version 0.4 -- added copy constructor to Mat and CTRL-C catcher
 *             -- added const-correct &operator[], &operator(), and *col methods.
 * Version 0.3 -- added sparse matrices to ordinary interface
 * Version 0.2 -- compiles on 4.6.3, no C++11 needed. scalar supports string.
 *
 * Copyright (c) 2013-4 Kui Tang <kuitang@gmail.com>
 *
 * Inspired by Rcpp and nr3matlab.h
 */

#pragma once
#include <algorithm>
//#include <cmath>
#include <cstdint>
#include <complex>
#include <exception>
#include <vector>
#include <string>

#include <signal.h>

#include "mex.h"
#include "matrix.h"


namespace mexcpp {
  /* * Ctrl-C signal catching infrastructure * */
  // RAII to change/restore signal handlers.
  // To use, instatiate before the main computation loop and check the
  // interrupted field at each loop iteration.
  // Based on http://linuxtoosx.blogspot.com/2010/10/ctrl-c-signal-catching-from-c-program.html
  struct SigintHandler {
    struct sigaction sa, osa;
    static short interrupted;
    // protect against multiple instantiations; cheaper singleton
    // static int allocs; = 0;

    static void sigprocCtrlC(int sig) {
      interrupted = 1;
      mexErrMsgTxt("Ctrl-C caught in mex file; interrupted");
    }

    SigintHandler() {
      sa.sa_handler = sigprocCtrlC;
      sigaction(SIGINT, &sa, &osa);
    }

    ~SigintHandler() {
      sigaction(SIGINT, &osa, &sa);
    }

  };

  short SigintHandler::interrupted = 0;

  /* Global variables for storing the signal handlers */
  // Call before the main computation loop
  void saveSigact() {
  }

  // Call before returning to MATLAB
  void restoreSigact() {
  }

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
  // These are NOT memory managed, because they are not classes!
  // If you are not directly putting them into pRhs, make sure to call
  // mxDestroyArray.

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

  // TODO: Check the generated code
  struct BaseMat {
    // M is rows, N is cols
    mxArray *pm;
    size_t N, M, length;
    // Do *I* own the memory, or does MATLAB?
    bool owned;

    BaseMat(const mxArray *p, bool owned_=false) : pm(const_cast<mxArray *>(p)), owned(owned_) {
      //mexPrintf("BaseMat: p = %lx\n", p);
      M = mxGetM(p);
      N = mxGetN(p);
      length = N * M;
    }

    BaseMat(size_t M_, size_t N_, bool owned_=false) :
      pm(0), M(M_), N(N_), length(N_ * M_), owned(owned_) { }

    ~BaseMat() {
      if (owned) {
        mxDestroyArray(pm);
      }
    }

    // Conversion cast operator
    operator mxArray *() { return pm; }

    inline size_t sub2ind(size_t r, size_t c) const { return c*M + r; }
  };

  // TODO: Specialize for a complex!
  template<class T>
  struct Mat : public BaseMat {
    int complexity;
    T *re;

    // Construct from MATLAB pointer
    Mat(const mxArray *p, bool owned_=false) : BaseMat(p, owned_) {
      checkTypeOrErr<T>(p);
      re = static_cast<T *>(mxGetData(p));
    }

    // Construct our own (in MATLAB's memory)
    Mat(size_t rows, size_t cols, bool owned_=false) : BaseMat(rows, cols, owned_) {
      pm = mxCreateNumericMatrix(rows, cols, ty<T>(), mxREAL);
      re = static_cast<T *>(mxGetData(pm));
    }

    // Construct from existing array (copy to MATLAB's memory)
    Mat(size_t rows, size_t cols, T *p, bool owned_=false) : BaseMat(rows, cols, owned_) {
      pm = mxCreateNumericMatrix(rows, cols, ty<T>(), mxREAL);
      re = static_cast<T *>(mxGetData(pm));
      std::copy(p, p + rows*cols, re);
    }

    // Copy constructor (matlab owns)
    Mat(const Mat<T> &m, bool owned_=false) : BaseMat(m.M, m.N, owned) {
      pm = mxCreateNumericMatrix(m.M, m.N, ty<T>(), mxREAL);
      re = static_cast<T *>(mxGetData(pm));
      std::copy(m.re, m.re + m.length, re);
    }

    // TODO: Figure out correct ownership semantics.
    // Construct from existing std::vector (copy to MATLAB's memory)
    Mat(const std::vector<T> v, bool colMat=true, bool owned_=false) : BaseMat(0, 0, owned_) {
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

    const T *col(size_t c) const { return re + sub2ind(0, c); }
    const T &operator[](size_t ind) const { return re[ind]; }
    const T &operator()(size_t r, size_t c) const { return re[sub2ind(r,c)]; }
  };

  // Homogeneous cell array (Mat, CellMat, classMat, etc.)
  template<class CellT>
  struct CellMat : public BaseMat {
    CellMat(const mxArray *p, bool owned_=false) : BaseMat(p, owned_) { checkTypeOrErr<mxCell>(p); }

    // Construct our own (in MATLAB's memory)
    CellMat(size_t rows, size_t cols, bool owned_) : BaseMat(rows, cols, owned_) {
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

    mxArray *ptr(size_t ind) const { return mxGetCell(pm, ind); }
    mxArray *ptr(size_t r, size_t c) const { return ptr(sub2ind(r, c)); }

    CellT operator[](size_t ind) const { return CellT(ptr(ind)); }
    CellT operator()(size_t r, size_t c) const { return CellT(ptr(sub2ind(r, c))); }
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

    StructMat(const mxArray *p, bool owned_) : BaseMat(p, owned_) {
      checkTypeOrErr<mxStruct>(p);
      nFields = mxGetNumberOfFields(p);
    }

    StructMat(size_t rows,
              size_t cols,
              std::vector<std::string> fieldNames,
              bool owned_=false) :
      BaseMat(rows, cols, owned_) {
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
    size_t nzMax;
    double *pr;
    mwIndex *ir, *jc;

    // Construct from MATLAB pointer
    SparseMat(mxArray *pm) : BaseMat(pm) {
      if (!mxIsDouble(pm) || !mxIsSparse(pm)) {
        mexErrMsgIdAndTxt("SparseMat:type", "matrix pm must be double and sparse");
      }

      setPointers();
      setSizes();
    }

    // Construct our own (in MATLAB's memory, I believe, but not positive).
    SparseMat(size_t rows,
              size_t cols,
              size_t nzMax,
              bool owned_=false) : BaseMat(rows, cols, owned_) {
      pm = mxCreateSparse(rows, cols, nzMax, mxREAL);
      setPointers();
    }

    // TODO: Determine whether mexCallMATLAB can alter the arguments. I don't think so.
    // Construct from iVec/jVec/wVec format.
    SparseMat(Mat<double> &iVec,
              Mat<double> &jVec,
              Mat<double> &wVec,
              size_t N=-1,
              size_t M=-1,
              bool owned_=false) : BaseMat(0, 0, owned_) {
      mxArray *lhs[1];
      int ret;
      if (N == -1 || M == -1) {
        // 3 argument call
        mxArray *rhs[] = { iVec.pm, jVec.pm, wVec.pm };
        ret = mexCallMATLAB(1, lhs, 3, rhs, "sparse");
      } else {
        mxArray *rhs[] = { iVec.pm, jVec.pm, wVec.pm, scalar(double(N)), scalar(double(M)) };
        ret = mexCallMATLAB(1, lhs, 5, rhs, "sparse");
        // Destroy the new new scalars we made
        mxDestroyArray(rhs[3]);
        mxDestroyArray(rhs[4]);
      }

      mxAssert(ret == 0, "mexCallMATLAB failed.");
      pm = lhs[0];

      setPointers();
      setSizes();
    }

  private:
    void setPointers() {
      pr = mxGetPr(pm);
      ir = mxGetIr(pm);
      jc = mxGetJc(pm);
    }

    void setSizes() {
      N = mxGetN(pm);
      M = mxGetM(pm);
      nzMax = mxGetNzmax(pm);
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
