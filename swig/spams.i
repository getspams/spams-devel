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
/* !!!!!!!!!!!
*  HACK en attendant une modif de Julien : les sorties sur cerr font segfault
En fait c'est un pb plus general de non init de donnees statiques de la libc++
(on le retrouve avec les exceptions levées par le c++)
il faut mettre -lstdc++ en t^ete des libs au link ou faire LD_PRELOAD=libstdc++.so.6
*/
//#ifdef SWIGR
//#define cerr cout
//#endif
/* !!!! */

#include "spams.h"
#ifdef DEBUG
#include "spams-tst.h"
#endif
%}

%define argout_vectors(ctype)
	Vector<ctype> **omiter
%enddef

// list of arguments of type INPLACE_MATRIX
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
     Matrix<ctype> *alpha
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
    Vector<ctype> *inner_weights
%enddef
%define inplace_spmatrices(ctype)
    SpMatrix<ctype> *A
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
enum constraint_type { L1COEFFS, L2ERROR, PENALTY, SPARSITY, PENALTY2};

INSTANTIATE_DATA(sparseProject)
INSTANTIATE_DATA(lassoD)
INSTANTIATE_DATA(lassoQq)

/**** dictLearn ****/
enum constraint_type_D { L2,  L1L2, L1L2FL, L1L2MU};
INSTANTIATE_DATA(alltrainDL)

/**** prox ****/
INSTANTIATE_DATA(fistaFlat)
INSTANTIATE_DATA(proximalFlat)
