%{
extern "C" {
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <Rembedded.h>
#include <Rinterface.h>
#include <R_ext/RS.h>
#include <R_ext/Error.h>
}
#include "spams.h"

/* ********************
* this is a hack to handle multiple values outrput
* return value is a vector initailized by swig in typemap out
* R_result_pos is the index in the vector : it must be set to 0 in typemap(out)
* typemap(argout) cannot be used without typemap(ou) in a function
*/

static int R_result_pos = 0;

SEXP appendOutput(SEXP value,SEXP result) {
     R_result_pos++;
     if(LENGTH(result) > R_result_pos)
     	SET_VECTOR_ELT(result,R_result_pos,value);
     return result;
}
/*    end of hack */
%}

#define myerr(msg,n) {Rf_error(msg,n); return R_NilValue;}
#define MYADD_OUTPUT_ARG(result, value)  r_ans = appendOutput(value, R_OutputValues);

%typemap(throws) const char * %{
	Rf_error("Runtime Error %s",$1); 
	return R_NilValue;
%}


/* One dimensional input arrays */
%define %vector_typemaps(R_TYPE,R_CAST,DATA_TYPE)

%typemap(typecheck, precedence=SWIG_TYPECHECK_POINTER)
    (Vector<DATA_TYPE> *INPLACE_VECTOR)
{
    $1 = (TYPEOF($input) == R_TYPE && Rf_isVector($input) ) ? 1 : 0;
}

%typemap(in) (Vector<DATA_TYPE> *INPLACE_VECTOR)
{
    SEXP rvec=$input;
    if (TYPEOF(rvec) != R_TYPE || ! Rf_isVector(rvec))
    {
        /*SG_ERROR("Expected Double Vector as argument %d\n", m_rhs_counter);*/
        myerr("Expected DATA_TYPE Vector as argument %d",$argnum);
    }

    $1 = new Vector<DATA_TYPE>((DATA_TYPE*) R_CAST(rvec), LENGTH(rvec));
}
%typemap(out) (Vector<DATA_TYPE> *) 
{
    R_result_pos = 0;
    R_len_t n = result->n();
    DATA_TYPE *data = result->rawX();
    SEXP rmat;
    PROTECT(rmat = Rf_allocMatrix(REALSXP,1,n));
   if (!rmat) error_return("Cannot alloc R matrix for return value");
   DATA_TYPE *rdata = (DATA_TYPE *)R_CAST(rmat);
   memcpy(rdata,data,n * sizeof(DATA_TYPE));
   delete result;
   $result = rmat;
   UNPROTECT(1);
 
}
%typemap(freearg) (Vector<DATA_TYPE> *INPLACE_VECTOR)
{
	delete arg$argnum;
}
//  ARGOUT
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY) (Vector<DATA_TYPE> **ARGOUT_VECTOR)
{
	$1 = 1;

}
%typemap(in) (Vector<DATA_TYPE> **ARGOUT_VECTOR)
	     (Vector<DATA_TYPE>  *data_temp)
{
    $1 = &data_temp;
}
%typemap(argout) (Vector<DATA_TYPE> **ARGOUT_VECTOR ) 
{
	  # test argout
	  if(data_temp$argnum != NULL) {
	    R_len_t n = data_temp$argnum->n();
	    DATA_TYPE *data = data_temp$argnum->rawX();
	    SEXP rmat;
	    PROTECT(rmat = Rf_allocMatrix(R_TYPE,1,n));
	   if (!rmat) myerr("Cannot alloc R matrix for arg %d",$argnum);
	   DATA_TYPE *rdata = (DATA_TYPE *)R_CAST(rmat);
	   memcpy(rdata,data,n * sizeof(DATA_TYPE));
	   delete data_temp$argnum;
	   MYADD_OUTPUT_ARG($result,rmat);
	   UNPROTECT(1);
        }
}

%enddef /* %vector_typemaps */

%define map_matrix(DATA_TYPE,R_TYPE,R_CAST)
     SEXP rmat=$input;
     SEXP dims = Rf_getAttrib(rmat,Rf_install("dim"));
    if (TYPEOF(rmat) != R_TYPE || LENGTH(dims) != 2)	
    {	
        /*SG_ERROR("Expected Double Vector as argument %d\n", m_rhs_counter);*/
        myerr("Expected Double dense matrix as argument %d",$argnum);
    }
    $1 = new Matrix<DATA_TYPE> ((DATA_TYPE *)R_CAST(rmat),Rf_nrows(rmat),Rf_ncols(rmat));

%enddef /* map_matrix */

/* full matrix input/output */
%define %matrix_typemaps(R_TYPE,R_CAST,DATA_TYPE)
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
           (Matrix<DATA_TYPE> *INPLACE_MATRIX)
{
	$1 = (TYPEOF($input) == R_TYPE) ? 1 : 0;

}
%typemap(in) (Matrix<DATA_TYPE> *INPLACE_MATRIX)
{
     map_matrix(DATA_TYPE,R_TYPE,R_CAST)
}
%typemap(out) (Matrix<DATA_TYPE> *) 
{
    R_result_pos = 0;
    R_len_t m = result->m();
    R_len_t n = result->n();
    DATA_TYPE *data = result->rawX();
    SEXP rmat;
    PROTECT(rmat = Rf_allocMatrix(REALSXP,m,n));
   if (!rmat) error_return("Cannot alloc R matrix for return value");
   DATA_TYPE *rdata = (DATA_TYPE *)R_CAST(rmat);
   memcpy(rdata,data,m * n * sizeof(DATA_TYPE));
   delete result;
   $result = rmat;
   UNPROTECT(1);
 
}

%typemap(freearg)
  (Matrix<DATA_TYPE> *INPLACE_MATRIX)
{
	delete arg$argnum;
}
//  ARGOUT
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY) (Matrix<DATA_TYPE> **ARGOUT_MATRIX)
{
	$1 = 1;

}
%typemap(in) (Matrix<DATA_TYPE> **ARGOUT_MATRIX)
	     (Matrix<DATA_TYPE>  *data_temp)
{
    $1 = &data_temp;
}
%typemap(argout) (Matrix<DATA_TYPE> **ARGOUT_MATRIX ) 
{
	  # test argout
	  if(data_temp$argnum != NULL) {
	    R_len_t m = data_temp$argnum->m();
	    R_len_t n = data_temp$argnum->n();
	    DATA_TYPE *data = data_temp$argnum->rawX();
	    SEXP rmat;
	    PROTECT(rmat = Rf_allocMatrix(REALSXP,m,n));
	   if (!rmat) myerr("Cannot alloc R matrix for arg %d",$argnum);
	   DATA_TYPE *rdata = (DATA_TYPE *)R_CAST(rmat);
	   memcpy(rdata,data,m * n * sizeof(DATA_TYPE));
	   delete data_temp$argnum;
	   MYADD_OUTPUT_ARG($result,rmat);
	   UNPROTECT(1);
        }
}

%enddef /* %matrix_typemaps */

%fragment("DSp_Check","header")
{
    const int check_sparse(SEXP arg) {
	SEXP ivec = Rf_getAttrib(arg,Rf_install("i"));
	SEXP xvec = Rf_getAttrib(arg,Rf_install("x"));
	SEXP dims = Rf_getAttrib(arg,Rf_install("Dim"));
//	return (Rf_isInteger(ivec) && (TYPEOF(xvec) == R_TYPE) && (LENGTH(dims) == 2) ? 1 : 0);
	return (Rf_isInteger(ivec) && (TYPEOF(xvec) == REALSXP) && (LENGTH(dims) == 2) ? 1 : 0);

    }
    const int check_matrix(SEXP arg) {
    	  SEXP dims = Rf_getAttrib(arg,Rf_install("dim"));
	  return ((TYPEOF(arg) == REALSXP) && (LENGTH(dims) == 2) ? 1 : 0);
    }
}
%define map_sparse(DATA_TYPE,R_TYPE,R_CAST)
     SEXP spmat=$input;
     SEXP ivec = Rf_getAttrib(spmat,Rf_install("i"));
     SEXP pvec = Rf_getAttrib(spmat,Rf_install("p"));
     SEXP xvec = Rf_getAttrib(spmat,Rf_install("x"));
     SEXP dims = Rf_getAttrib(spmat,Rf_install("Dim"));

     if (TYPEOF(xvec) != R_TYPE || LENGTH(dims) != 2)	
    {	
        /*SG_ERROR("Expected Double Vector as argument %d\n", m_rhs_counter);*/
        myerr("Expected DATA_TYPE sparse matrix as argument %d",$argnum);
    }
     int *pB = INTEGER(pvec);
     int *pE = pB + 1;
     int *idims = INTEGER(dims);
     $1 = new SpMatrix<DATA_TYPE> ((DATA_TYPE *)R_CAST(xvec),INTEGER(ivec),pB,pE,idims[0],idims[1],LENGTH(xvec));
	
%enddef /* map_sparse */

%define %spmatrix_typemaps(R_TYPE,R_CAST,DATA_TYPE)
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
	fragment="DSp_Check")
        (SpMatrix<DATA_TYPE> *INPLACE_SPMATRIX)
{
	$1 = check_sparse((SEXP) $input);

}
%typemap(in) (SpMatrix<DATA_TYPE> *INPLACE_SPMATRIX)
{
     map_sparse(DATA_TYPE,R_TYPE,R_CAST)
}
%typemap(freearg)
  (SpMatrix<DATA_TYPE> *INPLACE_SPMATRIX)
{
	delete arg$argnum;
}

%typemap(out) (SpMatrix<DATA_TYPE> *) 
{
    R_result_pos = 0;
    R_len_t m = result->m();
    R_len_t n = result->n();
    R_len_t nzmax = result->nzmax();
    SEXP indptr, indices, vdata, dims, output;
    PROTECT(indptr = Rf_allocVector(INTSXP,n + 1));
    PROTECT(indices = Rf_allocVector(INTSXP,nzmax));
    PROTECT(vdata = Rf_allocVector(REALSXP,nzmax));
    PROTECT(dims = Rf_allocVector(VECSXP,2));
    SET_VECTOR_ELT(dims,0,Rf_ScalarInteger(m));
    SET_VECTOR_ELT(dims,1,Rf_ScalarInteger(n));

    DATA_TYPE *xdata = result->v();
    memcpy(R_CAST(vdata),xdata,nzmax * sizeof(DATA_TYPE));
    int *pB = result->pB();
    int *r = result->r();
    memcpy(INTEGER(indices),r,nzmax * sizeof(int));
    memcpy(INTEGER(indptr),pB,(n + 1) * sizeof(int));
    
    PROTECT(output = Rf_allocVector(VECSXP,4));
    SET_VECTOR_ELT(output,0,indptr);
    SET_VECTOR_ELT(output,1,indices);
    SET_VECTOR_ELT(output,2,vdata);
    SET_VECTOR_ELT(output,3,dims);
    delete result;
    $result = output;
    UNPROTECT(5);
}
%enddef /* %spmatrix_typemaps */

%define %dspmatrix_typemaps(R_TYPE,R_CAST,DATA_TYPE)
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
	fragment="DSp_Check")
        (AbstractMatrixB<DATA_TYPE> *INPLACE_DSPMATRIX)
{
	SEXP pvec = Rf_getAttrib($input,Rf_install("p"));
	if(pvec == R_NilValue) {
	    $1 = check_matrix((SEXP) $input);
	 } else {
	    $1 = check_sparse((SEXP) $input);
	 }

}
%typemap(in) (AbstractMatrixB<DATA_TYPE> *INPLACE_DSPMATRIX)
{
    SEXP pvec = Rf_getAttrib($input,Rf_install("p"));
    if(pvec ==  R_NilValue) {
       map_matrix(DATA_TYPE,R_TYPE,R_CAST)
    } else {
       map_sparse(DATA_TYPE,R_TYPE,R_CAST)
    }
}
%typemap(freearg)
  (AbstractMatrixB<DATA_TYPE> *INPLACE_DSPMATRIX)
{
	delete arg$argnum;
}

%enddef /* %dspmatrix_typemaps */

%define %datamatrix_typemaps(R_TYPE,R_CAST,DATA_TYPE)
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
	fragment="DSp_Check")
        (Data<DATA_TYPE> *INPLACE_DATAMATRIX)
{
	SEXP pvec = Rf_getAttrib($input,Rf_install("p"));
	if(pvec == R_NilValue) {
	    $1 = check_matrix((SEXP) $input);
	 } else {
	    $1 = check_sparse((SEXP) $input);
	 }

}
%typemap(in) (Data<DATA_TYPE> *INPLACE_DATAMATRIX)
{
    SEXP pvec = Rf_getAttrib($input,Rf_install("p"));
    if(pvec ==  R_NilValue) {
       map_matrix(DATA_TYPE,R_TYPE,R_CAST)
    } else {
       map_sparse(DATA_TYPE,R_TYPE,R_CAST)
    }
}
%typemap(freearg)
  (Data<DATA_TYPE> *INPLACE_DATAMATRIX)
{
	delete arg$argnum;
}

%enddef /* %datamatrix_typemaps */

/* ********************** */

%vector_typemaps(REALSXP,REAL,float)
%vector_typemaps(REALSXP,REAL,double)
%vector_typemaps(INTSXP,INTEGER,int)

%matrix_typemaps(REALSXP,REAL,float)
%matrix_typemaps(REALSXP,REAL,double)
%matrix_typemaps(LGLSXP,LOGICAL,bool)

%spmatrix_typemaps(REALSXP,REAL,float)
%spmatrix_typemaps(REALSXP,REAL,double)
%spmatrix_typemaps(LGLSXP,LOGICAL,bool)

%dspmatrix_typemaps(REALSXP,REAL,float)
%dspmatrix_typemaps(REALSXP,REAL,double)

%datamatrix_typemaps(REALSXP,REAL,float)
%datamatrix_typemaps(REALSXP,REAL,double)

// must instatiate templates with different names, else generated R file will do bad dispatching
%define INSTANTIATE_DATA( f_name )
%template(## f_name) _ ## f_name<double>;
//%template(f_ ## f_name) _ ## f_name<float>;
%enddef
