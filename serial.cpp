#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>   // only used for setting random seed, for creating matrix A, vectors x and B.
#include <math.h>
#include <iostream>
#include "common.h"
using namespace std;

#define DEFAULT_SIZE 64  // square matrix A side length
#define SHIFTSIZE 0.00000000001



int main(int argc, char *argv[]){
    int SIZE = DEFAULT_SIZE;
    if(argc==2){ 
	SIZE = atoi(argv[1]);
    }

    input_test(SIZE);
	
    srand(time(NULL));
    double *A_mat = (double *)calloc(SIZE*3,sizeof(double));
    double *x_vec = (double *)malloc(SIZE*sizeof(double));
    double *B_vec = (double *)malloc(SIZE*sizeof(double));

    create_tridiagonal(A_mat,SIZE);

    if(SIZE<=16){ 
	print_Tridiag(A_mat,SIZE);
    }else{
	print_tridiagonal(A_mat,SIZE);
    }

    make_answer_to_solve(A_mat,x_vec,B_vec,SIZE);

    mainSPIKE(A_mat,B_vec,SIZE);

    vector_subtract(x_vec, B_vec, SIZE);
}//end of main
