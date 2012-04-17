

/* ***************************
   macros to apply typemaps
************************** */


%define APPLY_INPLACE_VECTOR( ctype )
%apply Vector<ctype> *INPLACE_VECTOR {
       inplace_vectors(ctype)
};
%enddef
%define APPLY_ARGOUT_VECTOR( ctype )
%apply Vector<ctype> **ARGOUT_VECTOR {
       argout_vectors(ctype)
};
%enddef


%define APPLY_INPLACE_MATRIX( ctype )
%apply Matrix<ctype> *INPLACE_MATRIX {
       inplace_matrices(ctype)
};
%enddef
%define APPLY_ARGOUT_MATRIX( ctype )
%apply Matrix<ctype> **ARGOUT_MATRIX {
       argout_matrices(ctype)
};
%enddef

%define APPLY_INPLACE_SPMATRIX( ctype )
%apply SpMatrix<ctype> *INPLACE_SPMATRIX {
       inplace_spmatrices(ctype)
};
%enddef
%define APPLY_INPLACE_DSPMATRIX( ctype )
%apply AbstractMatrixB<ctype> *INPLACE_DSPMATRIX {
       inplace_dspmatrices(ctype)
};
%enddef
%define APPLY_INPLACE_DATAMATRIX( ctype )
%apply Data<ctype> *INPLACE_DATAMATRIX {
       inplace_datamatrices(ctype)
};
%enddef

APPLY_INPLACE_VECTOR(float)
APPLY_INPLACE_VECTOR(double)
APPLY_INPLACE_VECTOR(int)
APPLY_ARGOUT_VECTOR(int)

APPLY_INPLACE_MATRIX(float)
APPLY_INPLACE_MATRIX(double)
APPLY_INPLACE_MATRIX(bool)
APPLY_ARGOUT_MATRIX(double)

APPLY_INPLACE_SPMATRIX(float)
APPLY_INPLACE_SPMATRIX(double)
APPLY_INPLACE_SPMATRIX(bool)

APPLY_INPLACE_DSPMATRIX(float)
APPLY_INPLACE_DSPMATRIX(double)

APPLY_INPLACE_DATAMATRIX(double)
