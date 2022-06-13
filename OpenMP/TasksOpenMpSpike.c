/*
 * MUST COMPILE WITH:
 *
 * module load gcc/5.4.0
 * gcc -fopenmp -O3 -Wall -pedantic TasksOpenmpSpike.c -o spike
 *
 * This program preallocates threads at the beginning of the program
 * and uses tasks to complete the calculations in each partition.
 *
 * benchmarking occurs in the program for the spike and thomas algorithms.
 * results printed at the end. 
 *
 * Use the respective sbatch file to run for different matrix sizes.
 * ./spike ##1 ##2, ##1 = size of matrix. must be a size 2^n. ##2 = number of processors
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>  
#include <math.h>

#ifdef _OPENMP
#include "omp.h"
#endif


#define DEFAULT_SIZE 64  
#define SHIFTSIZE 0.00000000001
#define DEPTH 10

void input_test_size(int size){
	if(size<4){ fprintf(stderr, "Matrix size too small.\n"); exit(0);}
	int temp = size;
	while(temp!=2){
	if(temp%2!=0){
        	fprintf(stderr, "Program only designed to run sizes 2^n\n");
		exit(0);
	}
	temp/=2;
	}
}


void input_test_numthreads(int nthreads){
int maxthreads = 1020;
if(nthreads>maxthreads) printf("Too many threads.\n\n");
}

void print_tridiagonal(double *A, int size){
    for(int i =0; i<3*size; i++){
	printf("%f ",A[i]);
	if((i+1)%3==0) printf("\n");
    }
}


void printlocalmatrix(double *V, double *W, int size,int blocksize){
	int vlower = 0;
	int vupper = blocksize-1;
	int vxpos = 2;
	int wxpos = -1;
	int wlower = 0;
	int wupper = 1;
	for(int down = 0; down<size; down++){
		for(int across = 0; across<size;across++){
			if(across== down){ printf("1.000000 ");
			}else if(across ==wxpos && down>=wlower&&down<=wupper){ printf("%f ",W[down]);
			}else if(across == vxpos && down>=vlower&&down<=vupper){ printf("%f ",V[down]);
				if(down==vupper){ vlower+=blocksize; vupper+=blocksize;vxpos+=blocksize;
				wlower+=blocksize; wupper+=blocksize;wxpos+=blocksize;
				}
			}else{
				printf("-------- ");
			}
		}printf("\n");
	}
}

void print_Tridiag(double *A,int size){
	printf("%f %f ",A[1],A[2]);
	for(int i =0; i<size-2;i++) printf("-------- "); printf("\n");
	for(int i =1;i<size-1; i++){
		for(int j = 0;j<size; j++){
			if(j==i-1){ 
				printf("%f %f %f ",A[3*i],A[3*i+1],A[3*i+2]);
				j+=2;
				continue;
			}
		printf("-------- ");
		}
		printf("\n");
	}
	for(int i =0; i<size-2;i++) printf("-------- "); 
	printf("%f %f ",A[3*size-3],A[3*size-2]);
	printf("\n");
}



void create_tridiagonal(double *A, int size){
    A[0] = 0;
    int i;
    for(i = 1; i<3*size-1; i++){
//	A[i] = (double)(rand()%9+1);
	A[i] = (double)rand()/(double)RAND_MAX;
// 	A[i] = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
   }
    A[i]=0;
}



void make_answer_to_solve(double *A, double *x, double *B,int size){
    for(int index = 0; index<size; index++) x[index] = (double)(rand()%9+1); 
    B[0] = A[1]*x[0] + A[2]*x[1];
    for(int i = 1; i<size;i++){
	B[i] = A[3*i]*x[i-1]+A[3*i+1]*x[i]+A[3*i+2]*x[i+1];	
    }
    B[size-1] = A[3*size-3]*x[size-2]+A[3*size-2]*x[size-1];
}


void vector_subtract(double *x, double *guess, int height){
    double accuracy_level=0.0001;
    double diff[height];
    double max=0; 
	for(int i =0; i<height; i++){
            diff[i] = x[i]-guess[i];
	    if(fabs(diff[i])>max) max = fabs(diff[i]);
                if(fabs(diff[i])>accuracy_level) printf("position:%d, \tresult: %f, \tdifference: %f, \tactual:%f\n",i,guess[i],diff[i],x[i]);
        }
	if(max<0.000001) printf("MAXDIFF: %f\n",max);
}


void Thomas(double *A, double *G, int size){
	A[2]/=A[1]; 
	G[0]/=A[1];  
	A[1]=1;	     
	double prevB=A[2];
	int i=0;
	for(i = 1; i<size-1; i++){
		A[3*i+1]-=A[3*i]*prevB;
		G[i]-=A[3*i]*G[i-1];
		A[3*i]=0;
		A[3*i+2]/=A[3*i+1];
		prevB=A[3*i+2];
		G[i]/=A[3*i+1];
		A[3*i+1]=1;
	}
	A[3*i+1]-=prevB*A[3*i]; 
	G[i]-=G[i-1]*A[3*i]; 
	A[3*i]=0;
	G[i]/=A[3*i+1];
	A[3*i+1]=1;
	for(i; i>0;i--){
		G[i-1]-=A[3*i-1]*G[i];
	}
}

void twoxtwocalc(double *A, double *G, double *V, double *W,int type){
double det;
double a = A[1];
double b = A[2];
double c = A[3];
double d = A[4];
det = a*d-b*c;

if(det==0){
	double shift = 0.00000000001;
	//double shift = SHIFTSIZE;
	if( (fabs(a)+fabs(c))>=(fabs(b)+fabs(d)) )
	a+=shift*(fabs(a)+fabs(c));
	else 
	a+=shift*(fabs(b)+fabs(d));
    det = a*d-b*c;
//    printf("diagonal boosting occurs.\n");
}
if(type==0||type==1){
    double w = A[0];
    W[0] = d*w/det;
    W[1] = -c*w/det;
}
if(type==0||type==-1){
    double v = A[5];
    V[0] = -v*b/det;
    V[1] = v*a/det;
}
double G0 = G[0];
double G1 = G[1];
G[0] = (d*G0-b*G1)/det;
G[1] = (a*G1-c*G0)/det;
}



void rSpike(double* A,double *B, double *V, double *W, int type, int size, int depth){
	int partition_size = size/2;
	if(partition_size ==2){
	    if(type==0){
	    twoxtwocalc(A,B,V,W,0);
	    twoxtwocalc(A+6,B+2,V+2,W+2,0);
	    }else{
		    if(type==-1){
	  		    twoxtwocalc(A,B,V,W,-1);
	   		    twoxtwocalc(A+6,B+2,V+2,W+2,0);
           	    }else{// if(type==1){
	   		    twoxtwocalc(A,B,V,W,0);
	    		    twoxtwocalc(A+6,B+2,V+2,W+2,1);
		    } 
	   }
	}else{		
		if(type==0){
	    		#pragma omp task if(depth>0)
			rSpike(A,B,V,W,0,partition_size,depth-1);
	    		#pragma omp task if(depth>0)
			rSpike(A+3*partition_size,B+partition_size,V+partition_size,W+partition_size,0,partition_size,depth-1);	
		}else{
			if(type==-1){
				#pragma omp task if(depth>0)
				rSpike(A,B,V,W,-1,partition_size,depth-1);	    		
				#pragma omp task if(depth>0)
				rSpike(A+3*partition_size,B+partition_size,V+partition_size,W+partition_size,0,partition_size,depth-1);
			}else{// if(type==1){
				#pragma omp task if(depth>0)
				rSpike(A,B,V,W,0,partition_size,depth-1);
				#pragma omp task if(depth>0)
				rSpike(A+3*partition_size,B+partition_size,V+partition_size,W+partition_size,1,partition_size,depth-1);
			}
		}
	}

#pragma omp taskwait


	double *localVspike = (double *)malloc(partition_size*sizeof(double));
	double *localWspike = (double *)malloc(partition_size*sizeof(double));
//	double localVspike[partition_size];
//	double localWspike[partition_size];		
	double centerA = B[partition_size-1];
	double centerB = B[partition_size];
	for(int i = 0; i<partition_size;i++){ localVspike[i]=V[i]; localWspike[i]=W[i+partition_size];}
	double closestV = localVspike[partition_size-1];
	double closestW = localWspike[0];
	double factor;
	if(1-closestV*closestW == 0){
	//printf("Larger block boost occurs.\n");
	//double shift = SHIFTSIZE;
		double shift = 0.000000001;
		double tmp;
		if(fabs(closestV)>=fabs(closestW)) tmp = fabs(closestV)+1;
		else tmp = fabs(closestW)+1;
		closestV/= tmp*shift+1;
		localVspike[partition_size-1] = closestV;
		V[partition_size-1]=closestV;
		factor = 1.0+(1.0/(tmp*shift));
	}else{
		factor = 1.0/(1-closestV*closestW);
	}

#pragma omp task untied
	{
	int index; 
	double temp = factor*(closestW*centerA-centerB);
	for(index = 0; index<partition_size-1;index++) B[index]+= localVspike[index]*temp;
	B[index] = (centerA-closestV*centerB)*factor;
	B[index+1] = -temp;
	temp = factor*(closestV*centerB-centerA);
	for(index+=2;index<partition_size*2;index++) B[index]+= localWspike[index-partition_size]*temp;
	}

#pragma omp task untied
	{
	V[partition_size] *=factor;
	double temp = V[partition_size];
	for(int i =0; i<partition_size;i++) V[i] *= -temp;
	for(int i =1; i<partition_size;i++) V[partition_size+i]+= localWspike[i]*closestV*temp;
	W[partition_size-1]*=factor;
	temp = W[partition_size-1];
	for(int i =0; i<partition_size-1;i++) W[i] += localVspike[i]*closestW*temp;
	for(int i = 0; i< partition_size; i++) W[i+partition_size]*=-temp;	
	}


//#pragma omp taskwait
}


//called when only one set of spikes remains.
//updates G for the last time, giving us our answer.
void finish(double *V, double *W, double *G,int size){
	size/=2;
	double centerA = G[size-1];
	double centerB = G[size];
	double cV = V[size-1];
	double cW = W[size];
	double factor = 1.0/(1-cV*cW);
	int index;
	double temp = factor*(cW*centerA-centerB);
	for(index = 0; index<size-1;index++){
	    G[index]+=V[index]*temp;
	}
	G[index] = (centerA-cV*centerB)*factor;
	G[index+1] = -temp;
	temp = factor*(cV*centerB-centerA);
	for(index+=2; index<size*2;index++){
	   G[index]+=W[index]*temp;
	}
}



void mainSPIKE(double *A, double *B, int size,int depth){
	//double V[size];
	//double W[size];
	double *V = (double *)malloc(size*sizeof(double));
	double *W = (double *)malloc(size*sizeof(double));
	if(size == 4){
	    twoxtwocalc(A,B,V,W,-1);
	    twoxtwocalc(A+6,B+2,V+2,W+2,1);
	    printlocalmatrix(V, W,size,2);
	}else{
	   #pragma omp task if(depth>0)
	    rSpike(A,B,V,W,-1,size/2,depth-1);
           #pragma omp task if(depth>0)
	    rSpike(A+3*(size/2),B+(size/2),V+(size/2),W+(size/2),1,size/2,depth-1);
	   
	#pragma omp taskwait
	}

	finish(V,W,B,size);
}



int main(int argc, char *argv[]){
    int SIZE = DEFAULT_SIZE;

    int nthreads=4;
    for(int i = SIZE; i>2; i/=2){
	nthreads+=2;
    }
    if(argc==2){ 
	SIZE = atoi(argv[1]);
    	 nthreads=2; for(int i = SIZE; i>2; i/=2) nthreads+=2;
    }else if(argc==3){
	SIZE = atoi(argv[1]);
	nthreads = atoi(argv[2]);
    }
    input_test_size(SIZE);
    input_test_numthreads(nthreads);
    srand(time(NULL)); //seed for random number creation.
#ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_nested(1);
//    omp_set_max_threads(nthreads);
    omp_set_num_threads(nthreads);
//	printf("%d\n",omp_get_num_threads());
#endif
  //  printf("nthreads: %d",nthreads);

    double A_mat[3*SIZE];
    double x_vec[SIZE];
    double B_vec[SIZE];
   double G_vec[SIZE];

    create_tridiagonal(A_mat,SIZE);
//    if(SIZE<=16){ print_Tridiag(A_mat,SIZE); }else if(SIZE<=256){ print_tridiagonal(A_mat,SIZE); }
    make_answer_to_solve(A_mat,x_vec,B_vec,SIZE);
    for(int index = 0; index<SIZE; index++) G_vec[index] = B_vec[index];

    int depth=0;
    for(int i=SIZE; i>2; i/=2){
	depth++;
    }
    double simulation_time = omp_get_wtime();
    
#pragma omp parallel num_threads(nthreads)
{
#pragma omp single	
mainSPIKE(A_mat,G_vec,SIZE,depth-4);
}
    simulation_time = omp_get_wtime()-simulation_time;
    printf("Simulation time: %g\n",simulation_time);


    simulation_time = omp_get_wtime();
    Thomas(A_mat,B_vec,SIZE);
    simulation_time = omp_get_wtime()-simulation_time;
    printf("Simulation time, Thomas: %g\n",simulation_time);

//   print_tridiagonal(A_mat,SIZE);
   vector_subtract(x_vec, G_vec, SIZE);

}//end of main
