#include "../include/BasicSet.h"
#include "../include/FunctionLib.h"

using namespace std;

void OverlapMatrix(gsl_matrix * result,gsl_matrix * vectors,int rows, int columns, double a, double b)
{
    gsl_vector * vectortemp1;
    gsl_vector * vectortemp2;

    vectortemp1 = gsl_vector_calloc(columns);
    vectortemp2 = gsl_vector_calloc(columns);

    double tempresult;
    int i,j;
    for(i=0;i<columns;i++)
    {
        gsl_matrix_get_col(vectortemp1,vectors,i);

        for(j=0;j<columns;j++){
            gsl_matrix_get_col(vectortemp2,vectors,j);
            tempresult = InnerProduct(vectortemp1,vectortemp2,rows,rows, a,b);
            gsl_matrix_set(result,i,j,tempresult);
        }
    }

    gsl_vector_free(vectortemp1);
    gsl_vector_free(vectortemp2);

}

void ModuleVector(gsl_vector * result, gsl_matrix * vectors,int rows,int columns,double a,double b)
{
    gsl_vector * temp;

    gsl_vector_set_zero(result);

    temp = gsl_vector_calloc(rows);

    int i;

    for(i=0;i<columns;i++)
    {
        gsl_matrix_get_col(temp,vectors,i);
        
        gsl_vector_set(result,i,InnerProduct(temp,temp,rows,rows,a,b));
    }

    gsl_vector_free(temp);
}

void SchmidtOrthogonalize(gsl_matrix* result,gsl_matrix* Vectors, int rows, int columns,double a,double b)
{
    //define, and allocate memories for the variables
    gsl_vector* vectortemp1;
    gsl_vector* vectortemp2;


    vectortemp1 = gsl_vector_alloc(rows);
    vectortemp2 = gsl_vector_alloc(rows);

    int i,j;

    double innerproductresult,projresult;

    double innerproduct;


    //schmidt orthogonalization

    for(i=0;i<columns;i++){
        //store the original vector of the i-th column
        gsl_matrix_get_col(vectortemp1,Vectors,i);

        for(j=0;j<i;j++){
            //read the processed vector and do the projection
            gsl_matrix_get_col(vectortemp2,result,j);

            innerproductresult = InnerProduct(vectortemp1,vectortemp2,rows,rows,a,b);

            projresult = innerproductresult;

            innerproductresult = InnerProduct(vectortemp2,vectortemp2,rows,rows,a,b);

            projresult = projresult / innerproductresult;

            gsl_vector_scale(vectortemp2,projresult);

            gsl_vector_sub(vectortemp1,vectortemp2);
        };

        //Normalize the vector

        innerproduct = InnerProduct(vectortemp1,vectortemp1,rows,rows,a,b);

        gsl_vector_scale(vectortemp1,1.0/sqrt(innerproduct));

        gsl_matrix_set_col(result,i,vectortemp1);
    };

    //free the memory
    gsl_vector_free(vectortemp1);
    gsl_vector_free(vectortemp2);
}


void Normalize(gsl_matrix* vectors, int rows, int columns,double a,double b)
{
    int i;
    gsl_vector* vectortemp1;
    double innerproduct;

    vectortemp1 = gsl_vector_calloc(rows);


    for(i=0;i<columns;i++){
        gsl_matrix_get_col(vectortemp1,vectors,i);
        innerproduct = InnerProduct(vectortemp1,vectortemp1,rows,rows,a,b);
        gsl_vector_scale(vectortemp1,1.0/sqrt(innerproduct));

        gsl_matrix_set_col(vectors,i,vectortemp1);
    };

    gsl_vector_free(vectortemp1);
}

void QuadraturePoint(gsl_matrix * m,gsl_vector * eval,gsl_vector * weight,int rows, int columns, double a, double b)
{

    //using the Golub-Welsch algorithm
    gsl_vector * tempvector1, * tempvector2;
    gsl_vector * xmultiplied1;
    gsl_vector * weighttemp;

    gsl_matrix * Jacobimatrix;
    gsl_matrix * Eigenvectors;
    gsl_eigen_symmv_workspace * eigenworkspace;

    double temp1,temp2;
    double cell;

    int i,j,k;

    Jacobimatrix = gsl_matrix_calloc(rows,columns);
    eigenworkspace = gsl_eigen_symmv_alloc(rows);
    tempvector1 = gsl_vector_calloc(rows);
    tempvector2 = gsl_vector_calloc(rows);
    Eigenvectors = gsl_matrix_calloc(rows,columns);
    weighttemp = gsl_vector_calloc(columns);


    //multiplied with x
    xmultiplied1 = gsl_vector_calloc(rows+1);


    for(i=0;i<rows;i++){
        for(j=0;j<i+1;j++){
                gsl_matrix_get_col(tempvector1,m,i);
                gsl_matrix_get_col(tempvector2,m,j);
                gsl_vector_set_zero(xmultiplied1);
                for(k=0;k<rows;k++){
                    gsl_vector_set(xmultiplied1,k+1,gsl_vector_get(tempvector1,k));
                }

                cell = InnerProduct(xmultiplied1,tempvector2,rows+1,rows,a,b);
                
                gsl_matrix_set(Jacobimatrix,i,j,cell);

        }
    }

    gsl_eigen_symmv(Jacobimatrix,eval,Eigenvectors,eigenworkspace);

    gsl_eigen_symmv_sort(eval,Eigenvectors,GSL_EIGEN_SORT_VAL_ASC);
    //Weight function over the domain
    temp1 = 1.0/(b-a);

    for(i=0;i<columns;i++){
        temp2 = gsl_pow_2(gsl_matrix_get(Eigenvectors,0,i));
        gsl_vector_set(weighttemp,i,temp1*temp2);
    }
    
    gsl_vector_memcpy(weight,weighttemp);

    gsl_matrix_free(Jacobimatrix);
    gsl_matrix_free(Eigenvectors);
    gsl_vector_free(tempvector1);
    gsl_vector_free(tempvector2);
    gsl_vector_free(xmultiplied1);
    gsl_vector_free(weighttemp);

    gsl_eigen_symmv_free(eigenworkspace);

}

void QuadraturePoint_withBasis(gsl_matrix * m,gsl_vector * eval,gsl_vector * weight, gsl_matrix * trsfm, int rows, int columns, double a, double b)
{

    //using the Golub-Welsch algorithm
    gsl_vector * tempvector1, * tempvector2;
    gsl_vector * xmultiplied1;
    gsl_vector * weighttemp;

    gsl_matrix * Jacobimatrix;
    gsl_matrix * Eigenvectors;
    gsl_eigen_symmv_workspace * eigenworkspace;

    double temp1,temp2;
    double cell;

    int i,j,k;

    Jacobimatrix = gsl_matrix_calloc(rows,columns);
    eigenworkspace = gsl_eigen_symmv_alloc(rows);
    tempvector1 = gsl_vector_calloc(rows);
    tempvector2 = gsl_vector_calloc(rows);
    Eigenvectors = gsl_matrix_calloc(rows,columns);
    weighttemp = gsl_vector_calloc(columns);


    //multiplied with x
    xmultiplied1 = gsl_vector_calloc(rows+1);


    for(i=0;i<rows;i++){
        for(j=0;j<i+1;j++){
                gsl_matrix_get_col(tempvector1,m,i);
                gsl_matrix_get_col(tempvector2,m,j);
                gsl_vector_set_zero(xmultiplied1);
                for(k=0;k<rows;k++){
                    gsl_vector_set(xmultiplied1,k+1,gsl_vector_get(tempvector1,k));
                }

                cell = InnerProduct(xmultiplied1,tempvector2,rows+1,rows,a,b);
                
                gsl_matrix_set(Jacobimatrix,i,j,cell);

        }
    }

    gsl_eigen_symmv(Jacobimatrix,eval,Eigenvectors,eigenworkspace);

    gsl_eigen_symmv_sort(eval,Eigenvectors,GSL_EIGEN_SORT_VAL_ASC);
    //Weight function over the domain
    temp1 = 1.0/(b-a);

    for(i=0;i<columns;i++){
        temp2 = gsl_pow_2(gsl_matrix_get(Eigenvectors,0,i));
        gsl_vector_set(weighttemp,i,temp1*temp2);
    }
    
    gsl_vector_memcpy(weight,weighttemp);
    gsl_matrix_memcpy(trsfm,Eigenvectors);

    gsl_matrix_free(Jacobimatrix);
    gsl_matrix_free(Eigenvectors);
    gsl_vector_free(tempvector1);
    gsl_vector_free(tempvector2);
    gsl_vector_free(xmultiplied1);
    gsl_vector_free(weighttemp);

    gsl_eigen_symmv_free(eigenworkspace);

}

void TransformationMatrix(gsl_matrix * target, gsl_matrix * CoefMatrix, gsl_vector * quadpoints, gsl_vector * weights,int length)
{
    int i,j,k;
    gsl_vector * CoefTemp;
    double sum;

    CoefTemp = gsl_vector_alloc(length);
    sum = 0.0;

    for(i=0;i<length;i++){
        for(j=0;j<length;j++){
            gsl_matrix_get_col(CoefTemp,CoefMatrix,j);
            for(k=0;k<length;k++)
                sum = sum + gsl_vector_get(CoefTemp,k)*gsl_pow_int(gsl_vector_get(quadpoints,i),k);

                /////////////////////////////////////////////////////
                ////
                //// The transformation matrix that transforms matrices between real space and DVR grids
                //// 
                //// Rewrite the following code for deploying weight function
                ////
                /////////////////////////////////////////////////////

            gsl_matrix_set(target,i,j,sqrt(gsl_vector_get(weights,i))*sum);

            sum = 0.0;

        }
    }

    gsl_vector_free(CoefTemp);
}


