#include <stdlib.h>
#include<iostream>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "common.h"
#include<vector>
#include<algorithm>
#include <bits/stdc++.h>
using namespace std;
double size;

#define DEFAULT_SIZE 64  // square matrix A side length
#define SHIFTSIZE 0.00000000001

//
//  timer
//
double read_timer( )
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

//
//Solver Routines
//

void input_test(int size){
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
				//if(down==wupper){ wlower+=blocksize; wupper+=blocksize;wxpos+=blocksize;}
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
	if (i%3==1)
		A[i] = 2;
	else
		A[i]= 1;
    }
    A[i]=0;
}

void make_answer_to_solve(double *A, double *x, double *B,int size){
    for(int index = 0; index<size; index++) x[index] = 1;//(double)(rand()%9+1); 
    B[0] = A[1]*x[0] + A[2]*x[1];
    for(int i = 1; i<size;i++){
	B[i] = A[3*i]*x[i-1]+A[3*i+1]*x[i]+A[3*i+2]*x[i+1];	
    }
    B[size-1] = A[3*size-3]*x[size-2]+A[3*size-2]*x[size-1];
}

double max(double *a, double *b,int size){
	double suma=0;
	double sumb=0;
	for(int i =0; i<size; i++){
   	    suma+=fabs(a[i]);
   	    sumb+=fabs(b[i]);
	}
	if(suma>=sumb) return suma+1; else return sumb+1;
}

void vector_subtract(double *x, double *guess, int height){
    double accuracy_level=0.005;
    double *diff = (double *)malloc(height*sizeof(double));
        for(int i =0; i<height; i++){
            diff[i] = x[i]-guess[i];
                if(fabs(diff[i])>accuracy_level) printf("position:%d, \tresult: %f, \tdifference: %f, \tactual:%f\n",i,guess[i],diff[i],x[i]);
        }
}





/*
void copy_vector(double *one, double *two,int size){
for(int i = 0; i<size; i++)two[i]  = one[i];
}
*/

void twoxtwocalc(double *A, double *G, double *V, double *W,int type){
double det;
double a = A[1];
double b = A[2];
double c = A[3];
double d = A[4];
det = a*d-b*c;

if(det==0){
	//double shift = 0.00000000001;
	double shift = SHIFTSIZE;
	if( (fabs(a)+fabs(c))>=(fabs(b)+fabs(d)) )
	a+=shift*(fabs(a)+fabs(c));
	else 
	a+=shift*(fabs(b)+fabs(d));
    det = a*d-b*c;
    printf("diagonal boosting occurs.\n");
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

void rSpike(double *A,double *B, double *V, double *W, int type, int size){
	int partition_size = size/2;
	if(partition_size ==2){
	    
	    if(type==0){
            //one processor/thread
	    twoxtwocalc(A,B,V,W,0);
	    //second processor
	    twoxtwocalc(A+6,B+2,V+2,W+2,0);
	    }else{
		    if(type==-1){
        		    //one processor
	  		    twoxtwocalc(A,B,V,W,-1);
	  		    //second processor
	   		    twoxtwocalc(A+6,B+2,V+2,W+2,0);
           	    }else{// if(type==1){
         		    //one processor
	   		    twoxtwocalc(A,B,V,W,0);
	    		    //second processor
	    		    twoxtwocalc(A+6,B+2,V+2,W+2,1);
		    } 
	   }
	}else{	
		if(type==0){
			//one processor
			rSpike(A,B,V,W,0,partition_size);
			//another processor
			rSpike(A+3*partition_size,B+partition_size,V+partition_size,W+partition_size,0,partition_size);
		}else{

		if(type==-1){
			rSpike(A,B,V,W,-1,partition_size);
			rSpike(A+3*partition_size,B+partition_size,V+partition_size,W+partition_size,0,partition_size);
		}else{// if(type==1){
			rSpike(A,B,V,W,0,partition_size);
			rSpike(A+3*partition_size,B+partition_size,V+partition_size,W+partition_size,1,partition_size);
		}
		}
	}
	double *localVspike= (double *)malloc(partition_size*sizeof(double));	
	double *localWspike= (double *)malloc(partition_size*sizeof(double));	
	double centerA = B[partition_size-1];
	double centerB = B[partition_size];
	for(int i = 0; i<partition_size;i++){ localVspike[i]=V[i]; localWspike[i]=W[i+partition_size];}
	double closestV = localVspike[partition_size-1];
	double closestW = localWspike[0];
	double factor;

	if(1-closestV*closestW == 0){
	printf("Larger block boost occurs.\n");
		double shift = SHIFTSIZE;
		//double shift = 0.00000000001;
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

	int index; 
	double temp = factor*(closestW*centerA-centerB);
	for(index = 0; index<partition_size-1;index++) B[index]+= localVspike[index]*temp;
	B[index] = (centerA-closestV*centerB)*factor;
	B[index+1] = -temp;
	temp = factor*(closestV*centerB-centerA);
	for(index+=2;index<partition_size*2;index++){
	B[index]+= localWspike[index-partition_size]*temp;
	}

	V[partition_size] *=factor;
	temp = V[partition_size];
	for(int i =0; i<partition_size;i++){
	    V[i] *= -temp;
	}
	for(int i =1; i<partition_size;i++){
	    V[partition_size+i]+= localWspike[i]*closestV*temp;
	}
	W[partition_size-1]*=factor;
	temp = W[partition_size-1];
	for(int i =0; i<partition_size-1;i++){
	    W[i] += localVspike[i]*closestW*temp;
	}
	for(int i = 0; i< partition_size; i++){
	    W[i+partition_size]*=-temp;
	}
}


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

void mainSPIKE(double *A, double *B, int size){
	double *V = (double *)calloc(size,sizeof(double));
	double *W = (double *)calloc(size,sizeof(double));
	if(size == 4){
            //one processor
	    twoxtwocalc(A,B,V,W,-1);
	    //second processor
	    twoxtwocalc(A+6,B+2,V+2,W+2,1);
	    printlocalmatrix(V, W,size,2);
	}else{
	    //one processor
	    rSpike(A,B,V,W,-1,size/2);
	    //second processor
	    rSpike(A+3*(size/2),B+(size/2),V+(size/2),W+(size/2),1,size/2);
	}
	//this processor
	finish(V,W,B,size);
}

void delete_alternative (vector<vector<double> > &vec, int start){
    int sz = vec.size()/2;
    for(int rp=start; rp<sz+start; ++rp) {
	    vec.erase(vec.begin()+rp);
     }
}


//
//  command line option processing
//
int find_option( int argc, char **argv, const char *option )
{
    for( int i = 1; i < argc; i++ )
        if( strcmp( argv[i], option ) == 0 )
            return i;
    return -1;
}

int read_int( int argc, char **argv, const char *option, int default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return atoi( argv[iplace+1] );
    return default_value;
}

char *read_string( int argc, char **argv, const char *option, char *default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return argv[iplace+1];
    return default_value;
}
//
//MPI Specific functions
//

