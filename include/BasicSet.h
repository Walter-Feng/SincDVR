#include <iostream>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>

double InnerProduct(gsl_vector* ,gsl_vector*,int,int,double,double);
double PolynomialIntegral(gsl_vector *, int, double, double);

void SchmidtOrthogonalize(gsl_matrix*,gsl_matrix*,int,int,double,double);
void OverlapMatrix(gsl_matrix*, gsl_matrix *,int, int, double, double);
void TransformationMatrix(gsl_matrix * target, gsl_matrix * CoefMatrix, gsl_vector * quadpoints, gsl_vector * weights, int length);

void PolynomialMul(gsl_vector*, gsl_vector *, gsl_vector *, int, int);
void PolynomialDifferential(gsl_vector* result, gsl_vector * input, int length); 
void ModuleVector(gsl_vector*, gsl_matrix *,int,int,double,double);

void Normalize(gsl_matrix*,int,int,double,double);
void QuadraturePoint(gsl_matrix *,gsl_vector *,gsl_vector *,int, int, double, double);
void QuadraturePoint_withBasis(gsl_matrix *,gsl_vector *,gsl_vector *,gsl_matrix *,int, int, double, double);
