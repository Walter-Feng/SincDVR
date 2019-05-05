#include "../include/BasicSet.h"

using namespace std;

double PolynomialIntegral(gsl_vector * coefficients, int n, double min, double max)
{
	int i;

	double result = 0;

    for(i=0;i<n;i++){
    	result += gsl_vector_get(coefficients,i) / ((double) i+1.0) * gsl_pow_int(max,i+1);
    	result -= gsl_vector_get(coefficients,i) / ((double) i+1.0) * gsl_pow_int(min,i+1);
    }

    return result;

}

void PolynomialMul(gsl_vector * result, gsl_vector * vector1, gsl_vector * vector2, int n1, int n2)
{
	int i,j;

	gsl_vector_set_zero(result);
	
	for(i=0;i<n1;i++)
		for(j=0;j<n2;j++)
			gsl_vector_set(result,i+j,gsl_vector_get(result,i+j) + gsl_vector_get(vector1,i) * gsl_vector_get(vector2,j));
}

void PolynomialDifferential(gsl_vector * result,gsl_vector * input, int length)
{
	int i;
	gsl_vector_set_zero(result);

	for(i=0;i<length-1;i++)
		gsl_vector_set(result,i,(double) (i+1) * gsl_vector_get(input,i+1));
	
}

//The definition of the inner product utilized in the program, which follows integrations of the polynomials 
double InnerProduct(gsl_vector* vector1, gsl_vector* vector2,int length1, int length2, double a,double b)
{
    gsl_vector * tempvector;

    double result;

	tempvector = gsl_vector_calloc(length1+length2-1);

    PolynomialMul(tempvector,vector1,vector2,length1,length2);

    result = PolynomialIntegral(tempvector,length1+length2-1,a,b);

	gsl_vector_free(tempvector);

    return result;
}

double gsl_vector_complex_module(gsl_vector_complex * vector, int length,double min, double max)
{
    int i;
    double result;
    gsl_complex complextemp;
    double temp;

    gsl_vector * vectortemp;

    vectortemp = gsl_vector_calloc(length);

    result = 0;

    for(i=0;i<length;++i){
        complextemp = gsl_vector_complex_get(vector,i);
        temp = gsl_complex_abs(complextemp);
        gsl_vector_set(vectortemp,length,temp);
    }

    result = InnerProduct(vectortemp,vectortemp,length,length,min,max);

    gsl_vector_free(vectortemp);

    return result;
}

double gsl_vector_complex_Euclidmodule(gsl_vector_complex * vector, int length)
{
    int i;
    double result;
	
    result = 0;

    for(i=0;i<length;i++)
        result += gsl_complex_abs2(gsl_vector_complex_get(vector,i));

        return result;
}