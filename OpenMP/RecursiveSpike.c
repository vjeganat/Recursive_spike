/*

This is the serial recursive version.
compile with


module load gcc/5.4.0
g++ -O3 RecursiveSpike.c -o rspike

can run on single compute node, 
with 
./rspike ##1: ##1 = size of array, size 2^n

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



/*
 * Checks to make sure user input is 2^n
*/
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

// for testing time
double read_timer(){
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if(!initialized){
        gettimeofday(&start,NULL);
        initialized = true;
    }
    gettimeofday(&end,NULL);
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}


//prints in 3 columns
void print_tridiagonal(double *A, int size){
    for(int i =0; i<3*size; i++){
	printf("%f ",A[i]);
	if((i+1)%3==0) printf("\n");
    }
}

//prints as a full square matrix. only useful for smaller matrix sizes.
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

// nicely prints spike matrix. will print spike matrix at any level, but must be sent block size.
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


//makes tridiagonal matrix.  
//Format is 1D matrix, in row ordering, skipping zeroes.  
//First and last element are 0, since there are only two elements in those rows.
// *Updated to be doubles between 0 and 1.
void create_tridiagonal(double *A, int size){
    A[0] = 0;
    int i;
    for(i = 1; i<3*size-1; i++){
//	A[i] = (double)(rand()%9+1);
//	A[i] = (double)rand()/(double)RAND_MAX;
 	A[i] = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
   }
    A[i]=0;
}

//creates vector, does matrix vector multiply to create B.
//currently the x elements are whole numbers.
void make_answer_to_solve(double *A, double *x, double *B,int size){
    for(int index = 0; index<size; index++) x[index] = (double)(rand()%9+1); 
    B[0] = A[1]*x[0] + A[2]*x[1];
    for(int i = 1; i<size;i++){
	B[i] = A[3*i]*x[i-1]+A[3*i+1]*x[i]+A[3*i+2]*x[i+1];	
    }
    B[size-1] = A[3*size-3]*x[size-2]+A[3*size-2]*x[size-1];
}

// for testing results against original x vector.
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

// unoptimized thomas algorithm. serial.
void Thomas(double *A, double *B, int size){
	A[2]/=A[1]; 
	B[0]/=A[1];  
	A[1]=1;	     
	double prevB=A[2];
	int i=0;
	for(i = 1; i<size-1; i++){
		A[3*i+1]-=A[3*i]*prevB;
		B[i]-=A[3*i]*B[i-1];
		A[3*i]=0;
		A[3*i+2]/=A[3*i+1];
		prevB=A[3*i+2];
		B[i]/=A[3*i+1];
		A[3*i+1]=1;
	}
	A[3*i+1]-=prevB*A[3*i]; 
	B[i]-=B[i-1]*A[3*i]; 
	A[3*i]=0;
	B[i]/=A[3*i+1];
	A[3*i+1]=1;
	for(i; i>0;i--){
		B[i-1]-=A[3*i-1]*B[i];
	}
}


// lowest level. calculates the 2x2 matrix. updates W,V,G
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


//Recursive Loop
void rSpike(double *A,double *B, double *V, double *W, int type, int size, int depth){

	int partition_size = size/2;

	if(partition_size ==2){
//center blocks
	    if(type==0){
	    twoxtwocalc(A,B,V,W,0);
	    twoxtwocalc(A+6,B+2,V+2,W+2,0);
	    }else{
//left blocks
		    if(type==-1){
	  		    twoxtwocalc(A,B,V,W,-1);
	   		    twoxtwocalc(A+6,B+2,V+2,W+2,0);
//right blocks
           	    }else{// if(type==1){
	   		    twoxtwocalc(A,B,V,W,0);
	    		    twoxtwocalc(A+6,B+2,V+2,W+2,1);
		    } 
	   }
	}else{	
//center blocks
		if(type==0){
			rSpike(A,B,V,W,0,partition_size,depth-1);
			rSpike(A+3*partition_size,B+partition_size,V+partition_size,W+partition_size,0,partition_size,depth-1);	
		}else{
//left blocks
			if(type==-1){
				rSpike(A,B,V,W,-1,partition_size,depth-1);
				rSpike(A+3*partition_size,B+partition_size,V+partition_size,W+partition_size,0,partition_size,depth-1);
//right blocks
			}else{// if(type==1){
				rSpike(A,B,V,W,0,partition_size,depth-1);
				rSpike(A+3*partition_size,B+partition_size,V+partition_size,W+partition_size,1,partition_size,depth-1);
			}
		}
	}

//
// This section starts being processed after 2x2 is reached and calculated.
// Calculates new spikes, updates G, returns to parent/ goes back up one level.
//
	double localVspike[partition_size];
	double localWspike[partition_size];		
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
	int index; 
	double temp = factor*(closestW*centerA-centerB);
	for(index = 0; index<partition_size-1;index++) B[index]+= localVspike[index]*temp;
	B[index] = (centerA-closestV*centerB)*factor;
	B[index+1] = -temp;
	temp = factor*(closestV*centerB-centerA);
	for(index+=2;index<partition_size*2;index++) B[index]+= localWspike[index-partition_size]*temp;
	V[partition_size] *=factor;
	temp = V[partition_size];
	for(int i =0; i<partition_size;i++) V[i] *= -temp;
	for(int i =1; i<partition_size;i++) V[partition_size+i]+= localWspike[i]*closestV*temp;
	W[partition_size-1]*=factor;
	temp = W[partition_size-1];
	for(int i =0; i<partition_size-1;i++) W[i] += localVspike[i]*closestW*temp;
	for(int i = 0; i< partition_size; i++) W[i+partition_size]*=-temp;	
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


//Start of recursive spike function.
void mainSPIKE(double *A, double *B, int size,int depth){
	double V[size];
	double W[size];
// skips straight to 2x2 if the input matrix is a 4x4.
	if(size == 4){
	    twoxtwocalc(A,B,V,W,-1);
	    twoxtwocalc(A+6,B+2,V+2,W+2,1);
//	    printlocalmatrix(V, W,size,2);
	}else{
	    rSpike(A,B,V,W,-1,size/2,depth-1);
	    rSpike(A+3*(size/2),B+(size/2),V+(size/2),W+(size/2),1,size/2,depth-1);
	}
	finish(V,W,B,size);
}



int main(int argc, char *argv[]){
    int SIZE = DEFAULT_SIZE;

    // easy to test different matrix sizes. default is at top, currently 64. 
    // if you want to run a different size, ./spike "matrixsize"
    // returns error if not 2^n
    if(argc==2){ 
	SIZE = atoi(argv[1]);
    }    
    input_test_size(SIZE);
    srand(time(NULL)); //seed for random number creation.

//At least with OpenMp, dynamically allocated arrays aren't as good for efficiency.
//Since the sizes don't change at all, its better to define them this way.
    double A_mat[3*SIZE];
    double x_vec[SIZE];
    double B_vec[SIZE];
    double G_vec[SIZE];

    create_tridiagonal(A_mat,SIZE);
//    if(SIZE<=16){ print_Tridiag(A_mat,SIZE); }else if(SIZE<=256){ print_tridiagonal(A_mat,SIZE); }
    make_answer_to_solve(A_mat,x_vec,B_vec,SIZE);
    for(int index = 0; index<SIZE; index++) G_vec[index] = B_vec[index];

    double simulation_time = read_timer();
    //Thomas(A_mat,G_vec,SIZE);
    int depth=0;
    for(int i=SIZE; i>2; i/=2){
	depth++;
    }
//    printf("%d\n",depth);
    mainSPIKE(A_mat,G_vec,SIZE,depth-2);
    simulation_time = read_timer()-simulation_time;
    printf("Simulation time: %g\n",simulation_time);

//   print_tridiagonal(A_mat,SIZE);
   vector_subtract(x_vec, G_vec, SIZE);


    simulation_time = read_timer();
    Thomas(A_mat,B_vec,SIZE);
    simulation_time = read_timer()-simulation_time;
    printf("Simulation time, Thomas: %g\n",simulation_time);

//   print_tridiagonal(A_mat,SIZE);
   vector_subtract(x_vec, G_vec, SIZE);




}//end of main
