#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>

///////////////////////// Definition of the structures, mainly used as parameters for functions ///////////////////////


///////////Meaning of Words:////////////////////
//
//  coef ---> coefficient
//  transl ---> tranlation
//  stddev ---> standard deviation
//  magnif ---> magnification of the variable, 
//  e.g.: sin(2x) where 2 is the magnification
//
//  in all functions the combination of magnification 
//  and translation is interpreted as
//  f( magnif( x - transl) ).
//
////////////////////////////////////////////////

//struct of chained list describing a compound function of several terms multiplying together, where ID refers to the type of funcion (in other words the 'true' struct the pointer 'target' is referring to) . 'flag' signifies this term is 'multiplying' (value 1) or 'dividing' (other value) the whole term.

typedef struct MulMap
{
    void * NEXT;
    void * target;
    int ID;
    int flag;
}MulMap;

//struct of chained list describing a series of multiplied compound function summed up together, where flag means the term is summing (value 1) or subtracting (other value) to the whole term (though coef in each struct describing the functions has the ability to describe this)

typedef struct SumMap
{
    void *NEXT;
    void *MulMaphead;
    int flag;
}SumMap;

typedef struct intresult
{
    double result;
    double error;
}intresult;


//ID: 1
typedef struct PolyFunc
{
    int length;
    double* coef;
    double transl;
    double magnif;
}PolyFunc;



//ID: 2
typedef struct GaussianFunc
{
    double mean;
    double stddev;
    double coef;
}GaussianFunc;

//ID: 3
typedef struct SinCosFunc
{
    double* sincoef;
    double* coscoef;
    double transl;
    double magnif;
    int length;
}SinCosFunc;

//ID: 4
typedef struct Exp
{
    double coef;
    double transl;
    double magnif;
}Exp;

//ID: 5
typedef struct Log
{
    double coef;
    double magnif;
    double base;
    double transl;
}Log;

//ID: 6
typedef struct LegendreP
{
    int l;
    double coef;
    double magnif;
    double transl;
}LegendreP;

double IDMapping(double x, int ID, void * param);
double FunctionMultiply(double x, void * HeadMulMap);
double FunctionSum(double x,void * summap);

double Combination(int n,int r);

double PolyFunction(double x, void * param);
double GaussianFunction(double x,void * param);
double SinCosFunction(double x,void * param);
double Exponential(double x,void * param);
double Logarithm(double x,void * param);

//Series of functions to help create the structs
void FuncStructConstruct(FILE * inputfile,void * target, int ID);
void MulMapConstruct(FILE * inputfile, MulMap * headpointer);
void SumMapConstruct(FILE * inputfile, SumMap * headpointer);


//Series of functions to help free the memory taken by the structs
void FuncStructDelete(void* target, int ID);
void MulMapDelete(void * target);
void SumMapDelete(void * target);


void OrthoPolyCoefGen(gsl_matrix * target,int rows,int columns,double a,double b);
void PolynomialTranslationCoef(gsl_vector * target,PolyFunc * param,int vectorlength);
void LegendrePCoefficient(PolyFunc * target, LegendreP * param);


//GSL trifles
void gsl_vector_complex_convert(gsl_vector * source, gsl_vector_complex * target, int length);
void gsl_matrix_complex_convert(gsl_matrix * source, gsl_matrix_complex * target, int rows, int columns);
void gsl_vector_poly_convert(PolyFunc * target,gsl_vector * source, int vectorlength);
void gsl_vector_complex_extract(gsl_vector_complex * source, gsl_vector * real, gsl_vector * imag, int length);
void gsl_matrix_complex_extract(gsl_vector_complex * source, gsl_matrix * real,gsl_matrix * imag, int rows, int columns);
void gsl_vector_complex_combine(gsl_vector * real, gsl_vector * imag, gsl_vector_complex * target);
void gsl_matrix_complex_combine(gsl_matrix * real,gsl_matrix * imag, gsl_matrix_complex * target);
void gsl_matrix_diag(gsl_matrix * target, gsl_vector * diag, int length);
void gsl_matrix_complex_diag(gsl_matrix_complex * target, gsl_vector_complex * diag, int length);
void gsl_matrix_mul(gsl_matrix * A, gsl_matrix *B, gsl_matrix * Result,int Acolumn,int Arow,int Bcolumn);
void gsl_matrix_complex_mul(gsl_matrix_complex * A, gsl_matrix_complex *B, gsl_matrix_complex * Result,int Acolumn,int Arow,int Bcolumn);
double gsl_vector_inner_product(gsl_vector * A, gsl_vector * B,int length);
gsl_complex gsl_vector_complex_product(gsl_vector_complex * A, gsl_vector_complex * B, int length);
gsl_complex gsl_vector_complex_inner_product(gsl_vector_complex * A, gsl_vector_complex * B, int length);
void gsl_vector_transform(gsl_vector * vec,gsl_matrix * trf,int length);
void gsl_vector_complex_transform(gsl_vector_complex * vec, gsl_matrix_complex * trf,int length);
void gsl_matrix_unitmatrix(gsl_matrix * m,int length);
void gsl_matrix_complex_unitmatrix(gsl_matrix_complex * m,int length);
void gsl_vector_complex_conjugate(gsl_vector_complex * v, int length);
void gsl_matrix_complex_conjugate(gsl_matrix_complex * m, int rows, int columns);
