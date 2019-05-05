#include "../include/BasicSet.h"
#include "../include/FunctionLib.h"
#include "../include/fileprint.h"

#define planck 1.0
#define mass 1.0

using namespace std;

int main(int argc, char const *argv[])
{

    //Handling Files
    FILE * input;
    FILE * output;

    //strings, chars
    char str[100] = "";


    //integers
    int grades;
    int recursiondepth;

    //loop int
    int i,j;

    int step;

    //flags
    int inputflag,outputflag;
    int wavedefflag,potentialdefflag;

    //doubles
    double min,max;
    double dt;
    double dx;

    //structs
    SumMap * wavefunction;
    SumMap * potential;

    PolyFunc * polyfunctemp;



    //gsl stuffs
    gsl_matrix * basis_kinetics;
    gsl_matrix * gsl_matrix_temp;
    gsl_matrix * hamiltonianmatrix;

    gsl_matrix_complex * gsl_matrix_complex_temp1;
    gsl_matrix_complex * gsl_matrix_complex_temp2;
    gsl_matrix_complex * hamiltonianmatrixconverted;


    gsl_vector * gridpoints;
    gsl_vector * weights;
    gsl_vector * gsl_vector_temp;
    gsl_vector * gsl_vector_temp2;
    gsl_vector * potentialgrid;
    gsl_vector * coefficient_real;
    gsl_vector * coefficient_imag;

    gsl_vector_complex * gridvalues1;
    gsl_vector_complex * gridvalues2;
    
    gsl_complex I;
    gsl_complex unit;

    GSL_SET_IMAG(&I,1);
    GSL_SET_REAL(&I,0);
    GSL_SET_REAL(&unit,1);
    GSL_SET_IMAG(&unit,0);

    gsl_permutation * gsl_permutation_temp1;


    //Initialization
    inputflag = 0;
    outputflag = 0;
    wavedefflag = 0;
    potentialdefflag = 0;

    wavefunction = new SumMap;
    potential = new SumMap;

    //when assigning memory to the struct pointers of 'function series', you need to set the pointer in these structs to NULL before you actually use them.
    //YES YES YES I AM NOOB I KNOW IT BUT I DON'T WANT TO (JAILED)
    polyfunctemp = new PolyFunc;
    polyfunctemp->coef = NULL;


    //Scanning the arguments, mainly getting the input/output file name/directory
    for(i=0;i<argc;i++)
    {
        if(strcmp(argv[i],"-i")==0){
            input = fopen(argv[i+1],"r");
            inputflag = 1;
        }
        else if(strcmp(argv[i],"-o")==0){
            output = fopen(argv[i+1],"r");
            strcpy(str,argv[i+1]);
            while(output!=NULL){
                strcat(str,".swp");
                fclose(output);
                output = fopen(str,"r");
            }
            fclose(output);
            output = fopen(str,"w");
            outputflag = 1;
        }
    }


    if(inputflag == 0){
        printf("NO INPUT FILE!\n");
        return 1;
    }
    
    if(outputflag == 0){
        output = fopen("output/NONAME.txt","r");
        strcpy(str,"output/NONAME.txt");
        while(output!=NULL){
            strcat(str,".swp");
            fclose(output);
            output = fopen(str,"r");
        }
        fclose(output);
        output = fopen(str,"w");
    }

    // argument scanning end

    // scanning input file

    while(fscanf(input,"%s",str) != EOF){
    // scanning basic constants

        if(strcmp(str,"min:")==0)
            fscanf(input,"%lf",&min);
        if(strcmp(str,"max:")==0)
            fscanf(input,"%lf",&max);
        if(strcmp(str,"grades:")==0)
            fscanf(input,"%d",&grades);
        if(strcmp(str,"dx:")==0)
            fscanf(input,"%lf",&dx);
        if(strcmp(str,"dt:")==0)
            fscanf(input,"%lf",&dt);
        if(strcmp(str,"recursiondepth:")==0)
            fscanf(input,"%d",&recursiondepth);

    //storing functions
        if(strcmp(str,"INIWAVEDEF")==0){
            while(fscanf(input,"%s",str)){
                if(strcmp(str,"SUMMAPCALL")==0){
                    SumMapConstruct(input,wavefunction);
                    inputflag = 1;
                }
                if(strcmp(str,"INIWAVEDEFEND")==0){
                    inputflag = 2;
                    break;
                }
            }
        }

        if(strcmp(str,"POTENTIALDEF")==0){
            while(fscanf(input,"%s",str)){
                if(strcmp(str,"SUMMAPCALL")==0){
                    SumMapConstruct(input,potential);
                    inputflag = 1;
                }
                if(strcmp(str,"POTENTIALDEFEND")==0){
                    inputflag = 2;
                    break;
                }
            }
        }
    }

    //assigning the memories & Initializations

    basis_kinetics = gsl_matrix_calloc(2*grades+1,2*grades+1);
    gsl_matrix_temp = gsl_matrix_calloc(2*grades+1,2*grades+1);
    hamiltonianmatrix = gsl_matrix_calloc(2*grades+1,2*grades+1);

    gsl_matrix_complex_temp1 = gsl_matrix_complex_calloc(2*grades+1,2*grades+1);
    gsl_matrix_complex_temp2 = gsl_matrix_complex_calloc(2*grades+1,2*grades+1);
    hamiltonianmatrixconverted = gsl_matrix_complex_calloc(2*grades+1,2*grades+1);

    gridvalues1 = gsl_vector_complex_calloc(2*grades+1);
    gridvalues2 = gsl_vector_complex_calloc(2*grades+1);

    gridpoints = gsl_vector_calloc(2*grades+1);
    weights = gsl_vector_calloc(2*grades+1);    
    gsl_vector_temp = gsl_vector_calloc(2*grades+1);
    gsl_vector_temp2 = gsl_vector_calloc(2*grades+1);
    coefficient_real = gsl_vector_calloc(2*grades+1);
    coefficient_imag = gsl_vector_calloc(2*grades+1);
    potentialgrid = gsl_vector_calloc(2*grades+1);

    gsl_permutation_temp1 = gsl_permutation_calloc(2*grades+1);

    //int series;
    int permutation_signum[2*grades+1];

    // Generate the orthonormal basis for the system

    // Store the grid points
    for(i=0;i<2*grades+1;i++)
        gsl_vector_set(gridpoints,i,((double) (i-grades))*dx);


    //convert the wavefunction to discrete point sets
    for(i=0;i<2*grades+1;i++){
        gsl_vector_set(gsl_vector_temp,i,FunctionSum(gsl_vector_get(gridpoints,i),(void *) wavefunction));
    }

    gsl_vector_complex_convert(gsl_vector_temp,gridvalues1,2*grades+1);


    //store the gridpoints for the potential
    for(i=0;i<2*grades+1;i++){
        gsl_vector_set(potentialgrid,i,FunctionSum(gsl_vector_get(gridpoints,i),potential));
    }
        
    //calculate the translation energy matrix for the basis
    for(i=0;i<2*grades+1;i++){
        for(j=0;j<2*grades+1;j++){
            if(i==j){
                gsl_matrix_set(basis_kinetics,i,j,gsl_pow_2(M_PI)/(6.0 * gsl_pow_2(dx)));
            }
            else
            {
                gsl_matrix_set(basis_kinetics,i,j,gsl_pow_int(-1.0,i-j)/(gsl_pow_2(dx) * gsl_pow_2(i-j)));
            }
            
        }
    }

    //obtain the hamiltonian matrix in DVR
    gsl_matrix_diag(gsl_matrix_temp,potentialgrid,2*grades+1);
    gsl_matrix_add(gsl_matrix_temp,basis_kinetics);
    gsl_matrix_memcpy(hamiltonianmatrix,gsl_matrix_temp);

    //calculate I \hbar + \frac{i dt}{2} H stored in gsl_matrix_complex_temp1, and 
    // I \hbar - \frac{i dt}{2} H stored in gsl_matrix_complex_temp2 
    gsl_matrix_complex_convert(hamiltonianmatrix,hamiltonianmatrixconverted,2*grades+1,2*grades+1);
    gsl_matrix_complex_scale(hamiltonianmatrixconverted,gsl_complex_mul_real(I,dt/2.0));
    gsl_matrix_complex_unitmatrix(gsl_matrix_complex_temp1,2*grades+1);
    gsl_matrix_complex_scale(gsl_matrix_complex_temp1,gsl_complex_mul_real(unit,planck));
    gsl_matrix_complex_add(gsl_matrix_complex_temp1,hamiltonianmatrixconverted);
    gsl_matrix_complex_unitmatrix(gsl_matrix_complex_temp2,2*grades+1);
    gsl_matrix_complex_scale(gsl_matrix_complex_temp2,gsl_complex_mul_real(unit,planck));
    gsl_matrix_complex_sub(gsl_matrix_complex_temp2,hamiltonianmatrixconverted);
    //decompose the gsl_matrix_complex_temp1 to prepare solving the linear system
    gsl_linalg_complex_LU_decomp(gsl_matrix_complex_temp1,gsl_permutation_temp1,permutation_signum);

    //Evolution

    //Initialize the discrete wavefunction 
    for(i=0;i<2*grades+1;i++){
        gsl_vector_set(gsl_vector_temp,i,FunctionSum(gsl_vector_get(gridpoints,i),wavefunction));
    }
    gsl_vector_complex_convert(gsl_vector_temp,gridvalues1,2*grades+1);

    //print the information
    fprintf(output,"\nHamiltonian Matrix:\n");
    gsl_matrix_fprint(output,hamiltonianmatrix,2*grades+1,2*grades+1,"%20.6f");
    fprintf(output,"\nGrid points:\n");
    gsl_vector_fprint(output,gridpoints,2*grades+1,"%20.6f");

    fprintf(output,"\nSIMULATION START\n");

    for(step=0;step<recursiondepth;step++){
        fprintf(output,"step = %d\ntime = %10.2f",step,((double) step)*dt);
        fprintf(output,"\ngrid values:\n");
        gsl_vector_complex_fprint(output,gridvalues1,2*grades+1,"%20.6f");
        gsl_vector_complex_memcpy(gridvalues2,gridvalues1);
        gsl_vector_complex_transform(gridvalues2,gsl_matrix_complex_temp2,2*grades+1);
        gsl_linalg_complex_LU_solve(gsl_matrix_complex_temp1,gsl_permutation_temp1,gridvalues2,gridvalues1);
    }

    fprintf(output,"\nFINALE\n");
    fprintf(output,"step = %d\ntime = %10.2f",step,((double) step)*dt);
    fprintf(output,"\ngrid values:\n");
    gsl_vector_complex_fprint(output,gridvalues1,2*grades+1,"%20.6f");

    fclose(input);
    fclose(output);

    return 0;
}
