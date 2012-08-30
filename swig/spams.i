%define DOCSTRING
"This module gives access to some functions of the spams C++ library.
The functions defined here should not be called directly.
Use of spams functions should only be done through module spams."
%enddef

#ifdef SWIGPYTHON
%module(docstring=DOCSTRING) spams_wrap
#endif
#ifdef SWIGR
%module spams
#endif

%{
#define SWIG_FILE_WITH_INIT

#include "spams.h"
#ifdef DEBUG
#include "spams-tst.h"
#endif
%}

%define argout_vectors(ctype)
	Vector<ctype> **omiter
%enddef

/*
   Following macros define which typemap must be applied to arguments of C++ functions are converted.
   args of type inplace_* may containt input data or are ready (allocated)
    to receive output data. Typemaps of type 'in' are applied to them.
   args of type argout_* are used when multiple return values are needed.
    they are of type **; their storage is allocated on return of the C++ function.
     Typemaps of type 'argout' are applied to them.
  The third type ('out') of typemaps used is applied to return values of C++ functions.
*****

*/
// list of arguments of type INPLACE_MATRIX
%define inplace_bool_matrices
	Matrix<bool> *B
%enddef
%define inplace_bool_spmatrices
    SpMatrix<bool> *groups,
    SpMatrix<bool> *groups_var
%enddef

%define inplace_matrices(ctype)
     Matrix<ctype> *A,
     Matrix<ctype> *B,
     Matrix<ctype> *D,
     Matrix<ctype> *D1,
     Matrix<ctype> *Q,
     Matrix<ctype> *q,
     Matrix<ctype> *U,
     Matrix<ctype> *V,
     Matrix<ctype> *X,
     Matrix<ctype> *m_A,
     Matrix<ctype> *m_B,
     Matrix<ctype> *XAt,
     Matrix<ctype> *Y,
     Matrix<ctype> *XY,
     Matrix<ctype> *AAt,
     Matrix<ctype> *alpha0,
     Matrix<ctype> *alpha,
     Matrix<ctype> *W
%enddef
%define argout_matrices(ctype)
     Matrix<ctype> **path,
     Matrix<ctype> **path2,
     Matrix<ctype> **path3,
     Matrix<ctype> **omA,
     Matrix<ctype> **omB
%enddef
%define inplace_vectors(ctype)
    Vector<ctype> *v,
    Vector<ctype> *b,
    Vector<ctype> *x,
    Vector<ctype> *inner_weights,
    Vector<ctype> *L,
    Vector<ctype> *eps,
    Vector<ctype> *Lambda,
    Vector<ctype> *eta_g,
    Vector<ctype> *own_variables,
    Vector<ctype> *N_own_variables,
    Vector<ctype> *groups
%enddef
%define inplace_spmatrices(ctype)
    SpMatrix<ctype> *A,
    SpMatrix<ctype> *alpha,
    SpMatrix<ctype> *groups
%enddef
%define inplace_dspmatrices(ctype)
    AbstractMatrixB<ctype> *D
%enddef

%define inplace_datamatrices(ctype)
    Data<ctype> *X
%enddef

#ifdef SWIGPYTHON
%include "py_typemaps.i"
%init %{
    import_array();
%}
#endif

#ifdef SWIGR
%include "R_typemaps.i"
#endif
%include "spamstools.i"
%include "exception.i"


%include <spams.h>
void im2col_sliding(Matrix<double>  *,Matrix<double>  *,int,int,bool);
#ifdef DEBUG
%include <spams-tst.h>
Matrix<double> *tst(Matrix<double> **,bool,AbstractMatrixB<double> *);
//int xtst(Matrix<double> **,bool );
SpMatrix<double> *xtst(Matrix<double> **,bool );
Matrix<double> *ztst(Data<double> *,Matrix<double> **,Matrix<double> **,Vector<int> **,bool);
#endif

// linalg
INSTANTIATE_DATA(sort)
INSTANTIATE_DATA(mult)
INSTANTIATE_DATA(AAt)
INSTANTIATE_DATA(XAt)
INSTANTIATE_DATA(applyBayerPattern)
INSTANTIATE_DATA(conjugateGradient)
INSTANTIATE_DATA(invSym)
INSTANTIATE_DATA(normalize)

/**** decomp ****/
enum constraint_type { L1COEFFS, L2ERROR, PENALTY, SPARSITY, L2ERROR2, PENALTY2};

INSTANTIATE_DATA(sparseProject)
INSTANTIATE_DATA(lassoD)
INSTANTIATE_DATA(lassoQq)
INSTANTIATE_DATA(lassoMask)
INSTANTIATE_DATA(lassoWeighted)
INSTANTIATE_DATA(omp)
INSTANTIATE_DATA(ompMask)
INSTANTIATE_DATA(somp)
INSTANTIATE_DATA(cd)
INSTANTIATE_DATA(l1L2BCD)
/**** dictLearn ****/
enum constraint_type_D { L2,  L1L2, L1L2FL, L1L2MU};
INSTANTIATE_DATA(alltrainDL)

/**** prox ****/
INSTANTIATE_DATA(fistaFlat)
INSTANTIATE_DATA(fistaTree)
INSTANTIATE_DATA(fistaGraph)
INSTANTIATE_DATA(proximalFlat)
INSTANTIATE_DATA(proximalTree)
INSTANTIATE_DATA(proximalGraph)
