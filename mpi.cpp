#include <mpi.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>   // only used for setting random seed, for creating matrix A, vectors x and B.
#include <math.h>
#include <iostream>
#include <vector>
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
    // Setup MPI
    int n_proc, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //Common Variable to all processors
    int p; //Partition size (varies with every step)
    int np; //number of partitions (varies with every step)
    int n= log(SIZE)/log(2);//number of steps required to solve the matrix
    int row_per_proc= SIZE/n_proc;//number of elements that each processor gets
    int ser_steps= log(n_proc)/log(2)-1;
    int par_steps= n-2-ser_steps;

    if (rank==0)
    cout<<"serial step "<<ser_steps<<endl<< "par step"<<par_steps<<endl;
    
    //Setup Data Partitioning across processors
    int * arow_per_proc= (int *) malloc(n_proc* sizeof(int));//Array of no. of rows of B to send
    int * row_offsets= (int*) calloc((n_proc+1), sizeof (int));//Array of offset of rows of B to pick to send to proc
    int * aele_per_proc= (int *) malloc(n_proc* sizeof(int));//Array of no. of elements of A to send
    int * ele_offsets= (int*) calloc((n_proc+1), sizeof (int));//Array of offset of elements of A to send
    double *B_vec = (double *)malloc(SIZE*sizeof(double));
    double *A_mat = (double *)calloc(SIZE*3,sizeof(double));
    double *x_vec = (double *)malloc(SIZE*sizeof(double));
    vector<vector <double>> V(row_per_proc/2, vector<double>(0));
    vector<vector <double>> W(row_per_proc/2, vector<double>(0));
    //Setup the number of rows array to do scatterv
     for (int i=0; i<n_proc; i++){
	arow_per_proc[i]=row_per_proc;
 	aele_per_proc[i]=3*row_per_proc;
    }
    //Setup the offset number of rows array to do scatterv
    for (int i=1; i<n_proc+1; i++){
    	row_offsets[i]= row_offsets[i-1]+row_per_proc;
	ele_offsets[i]= 3*(row_offsets[i-1]+row_per_proc);
    }
    
    //if (rank==0)
    //	for (int i=0; i<n_proc+1; i++)
    //		cout<<row_offsets[i]<<endl;
    //

    //Allocate Storage for local partition
    double *A = (double*) malloc( 3*row_per_proc* sizeof(double));
    double *G = (double*) malloc(row_per_proc* sizeof (double));

    if (rank == 0){
             
    create_tridiagonal(A_mat,SIZE);
   
   // if(SIZE<=16){ 
   //     print_Tridiag(A_mat,SIZE);
   // }else{
   //     print_tridiagonal(A_mat,SIZE);
   // }

    make_answer_to_solve(A_mat,x_vec,B_vec,SIZE);
     
   }
   
   MPI_Scatterv(B_vec,arow_per_proc,row_offsets, MPI_DOUBLE, G, row_per_proc,MPI_DOUBLE, 0, MPI_COMM_WORLD );
   MPI_Scatterv(A_mat,aele_per_proc,ele_offsets, MPI_DOUBLE, A, 3*row_per_proc,MPI_DOUBLE, 0, MPI_COMM_WORLD );
 

  // if (rank==0){
  //      cout<<"Global "<<endl;
  //      print_tridiagonal(A_mat,SIZE);}
  //      MPI_Barrier(MPI_COMM_WORLD);
   //if (rank==0){
   //     cout<<"rank0  "<<endl;
   //     print_tridiagonal(A,row_per_proc);}
   //     MPI_Barrier(MPI_COMM_WORLD);
  // if (rank==1){
  //      cout<<"rank1  "<<endl;
  //      print_tridiagonal(A,row_per_proc);}
  //      MPI_Barrier(MPI_COMM_WORLD);
  // if (rank==2){
  //      cout<<"rank2  "<<endl;
  //      print_tridiagonal(A,row_per_proc);}
  //      MPI_Barrier(MPI_COMM_WORLD);
  // if (rank==3){
  //      cout<<"rank3  "<<endl;
  //      print_tridiagonal(A,row_per_proc);}
  //      MPI_Barrier(MPI_COMM_WORLD);
  //if (rank==0){
  //	cout<<"Global RHS"<<endl;
  //      for (int i=0; i<SIZE; i++)
  //      	cout<< B_vec[i]<<endl;
  //      		}
  // MPI_Barrier(MPI_COMM_WORLD);


// if (rank==1){
// 	cout<<"G for  "<<rank <<endl;
// 	for (int i=0; i<row_per_proc; i++)
//		cout<< G[i]<<endl;
//	  		}
//  MPI_Barrier(MPI_COMM_WORLD);

  //Forming initial Spikes in each processor
  //Completely Parallel
   int partition_proc= row_per_proc/2;
   int partition_size= row_per_proc/partition_proc;
   for (int i=0; i<partition_proc; i++){
   	double det;
	double a = A[1+6*i];
	double b = A[2+6*i];
	double c = A[3+6*i];
	double d = A[4+6*i];
	det = a*d-b*c;


	//double shift = 0.00000000001;
	double shift = SHIFTSIZE;
	if( (fabs(a)+fabs(c))>=(fabs(b)+fabs(d)) )
		a+=shift*(fabs(a)+fabs(c));
	else 
		a+=shift*(fabs(b)+fabs(d));
    	det = a*d-b*c;

   	double w = A[0+6*i];
	W[i].insert(W[i].end(),{d*w/det,-c*w/det});
    	double v = A[5+6*i];
	V[i].insert(V[i].end(),{-v*b/det,v*a/det});

	double G0 = G[0+2*i];
	double G1 = G[1+2*i];
	G[0+2*i] = (d*G0-b*G1)/det;
	G[1+2*i] = (a*G1-c*G0)/det;
    }
//if (rank==1){
// 	cout<<"G after first step for  "<<rank <<endl;
// 	for (int i=0; i<row_per_proc; i++)
//		cout<< G[i]<<endl;
//	  		}
//  MPI_Barrier(MPI_COMM_WORLD);
//
//    if (rank==1){
//        cout<<"Rank "<<rank<<" after pre-processing"<<endl;
//        for (int j=0; j<row_per_proc/2; j++){
//        	cout<<"W "<<j<<endl;
//        	for (int i=0; i<W[j].size(); i++)
//        		cout<< W[j][i]<<endl;
//          		}
//        	}
    MPI_Barrier(MPI_COMM_WORLD);

    // Generate spikes for next step

    // Adding zeros in alternative V and W to serve as future spikes
   for (int step=0; step<par_steps; step++){ 
    
    
    vector <double> zrs(V[0].size(),0);
    for (int i=1; i<partition_proc; i=i+2)
    	V[i].insert(V[i].begin(),zrs.begin(),zrs.end());

    for (int i=0; i<partition_proc; i=i+2)
    	W[i].insert(W[i].end(),zrs.begin(),zrs.end());

//    if (rank==0){
// 	cout<<"Rank "<<rank<<" after zeros"<<endl;
//	for (int j=1; j<V.size()/V[0].size(); j++){
//		cout<<"V "<<j<<endl;
// 		for (int i=0; i<W[j].size(); i++)
//			cout<< W[j][i]<<endl;
//	  		}
//		}
//	MPI_Barrier(MPI_COMM_WORLD);
    partition_proc= partition_proc/2;
    partition_size= row_per_proc/partition_proc;

     for (int i=0; i<partition_proc; i++){
	//Finding center elements of spikes
	double det;
	double a = 1;
	double b = V[i*2].back();
	double d = 1;
	double c = W[2*i+1].front();
	det = a*d-b*c;
	
	double w0= W[2*i][partition_size/2-1];
	double w1 = W[2*i][partition_size/2];
	W[2*i][partition_size/2-1]=(d*w0-b*w1)/det;
	W[2*i][partition_size/2]= (a*w1-c*w0)/det;
	
	for (int j=1; j<partition_size/2; j++){
		W[2*i][partition_size/2-1-j]= W[2*i][partition_size/2-1-j]-V[i*2].end()[-j-1]*W[2*i][partition_size/2];
		W[2*i][partition_size/2+j]=W[2*i][partition_size/2+j]-W[2*i+1][j]*W[2*i][partition_size/2-1];
	}
	double v0= V[2*i+1][partition_size/2-1];
    	double v1 = V[2*i+1][partition_size/2];
	V[2*i+1][partition_size/2-1]= (d*v0-b*v1)/det;
	V[2*i+1][partition_size/2]= (a*v1-c*v0)/det;
	for (int j=1; j<partition_size/2; j++){
		V[2*i+1][partition_size/2-1-j]=V[2*i+1][partition_size/2-1-j]-V[i*2].end()[-j-1]*V[2*i+1][partition_size/2];
		V[2*i+1][partition_size/2+j]=V[2*i+1][partition_size/2+j]-W[2*i+1][j]*V[2*i+1][partition_size/2-1];
		}
	double G0 = G[partition_size/2+partition_size*i-1];
	double G1 = G[partition_size/2+partition_size*i];
	G[partition_size/2+partition_size*i-1] = (d*G0-b*G1)/det;
	G[partition_size/2+partition_size*i] = (a*G1-c*G0)/det;

	for (int j=1; j<partition_size/2; j++){
		G[partition_size/2+partition_size*i-1-j]=G[partition_size/2+partition_size*i-1-j]-V[i*2].end()[-j-1]*G[partition_size/2+partition_size*i];
		G[partition_size/2+partition_size*i+j]=G[partition_size/2+partition_size*i+j]-W[2*i+1][j]*G[partition_size/2+partition_size*i-1];
		}


     }
//if (rank==1){
// 	cout<<"G after second step for  "<<rank <<endl;
// 	for (int i=0; i<row_per_proc; i++)
//		cout<< G[i]<<endl;
//	  		}
//  MPI_Barrier(MPI_COMM_WORLD);


//    if (rank==0){
// 	cout<<"Rank "<<rank<<" after second step"<<endl;
//	for (int j=1; j<2*partition_proc; j++){
//		cout<<"W "<<j<<endl;
// 		for (int i=0; i<W[j].size(); i++)
//			cout<< W[j][i]<<endl;
//	  		}
//		}
//	MPI_Barrier(MPI_COMM_WORLD);

    delete_alternative(V,0);
    delete_alternative(W,1);
   

        }

  for (int proc=0; proc<n_proc; proc++){
  if (rank==proc){
 	cout<<"Rank "<<rank<<"W after end of parallel"<<endl;
	for (int j=0; j<partition_proc; j++){
		cout<<"W "<<j<<endl;
 		for (int i=0; i<W[j].size(); i++)
			cout<< W[j][i]<<endl;
	  		}
		}
  if (rank==proc){
 	cout<<"Rank "<<rank<<"V after end of parallel"<<endl;
	for (int j=0; j<partition_proc; j++){
		cout<<"V "<<j<<endl;
 		for (int i=0; i<V[j].size(); i++)
			cout<< V[j][i]<<endl;
	  		}
		}

  if (rank==proc){
 	cout<<"Rank"<<rank <<"G after end of parallel" <<endl;
 	for (int i=0; i<row_per_proc; i++)
		cout<< G[i]<<endl;
	  		}
  MPI_Barrier(MPI_COMM_WORLD);
}



//Serial - parallel portion
    for (int step=0; step<1; step++){
	vector <proc_pair> pair(0);

	for (int i= pow(2,step)-1; i<n_proc; i=i+pow(2, step+1)){
		proc_pair temp;
		temp.up=i;
		temp.down=i+1;
		pair.push_back(temp);
		}
		MPI_Barrier(MPI_COMM_WORLD);
//Communication against processors to solve the remaining steps
//Getting the center pieces from nearby processors

    double send_intr[3];
    double recv_intr[3];

    for (int proc_intr=0; proc_intr< 1; proc_intr++){
        
    if (rank==pair[proc_intr].up){
    send_intr[0]=V[0].back();
    send_intr[1]=W[0].back();
    send_intr[2]= G[row_per_proc-1];
  
    MPI_Sendrecv(send_intr, 3, MPI_DOUBLE, pair[proc_intr].down, 0, recv_intr, 3, MPI_DOUBLE, pair[proc_intr].down, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    	}

    if (rank==pair[proc_intr].down){
    send_intr[0]=V[0].front();
    send_intr[1]=W[0].front();
    send_intr[2]= G[0];
    MPI_Sendrecv(send_intr, 3, MPI_DOUBLE, pair[proc_intr].up, 0, recv_intr, 3, MPI_DOUBLE, pair[proc_intr].up, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    }
    }

//Solving the center pieces together

   for (int proc_intr=0; proc_intr< 1; proc_intr++){
   
   if (rank==pair[proc_intr].up){
	double det;
	double a = 1;
	double b = V[0].back();
	double d = 1;
	double c = recv_intr[1];
	det = a*d-b*c;
	
	double w0= W[0].end()[-1];
	double w1 = 0.0;
	double wx0=(d*w0-b*w1)/det;
	double wx1= (a*w1-c*w0)/det;
	
	double v0= 0.0;
    	double v1 = recv_intr[0];
	double vx0= (d*v0-b*v1)/det;
	double vx1= (a*v1-c*v0)/det;
	
	double G0 = G[row_per_proc-1];
	double G1 = recv_intr[2];
	double gx0= (d*G0-b*G1)/det;
	double gx1= (a*G1-c*G0)/det;

	cout<<"Rank "<<rank <<" Vx0 "<<vx0<<" vx1 "<<vx1<<endl;
	cout<<"Rank "<<rank <<" Wx0 "<<wx0<<" Wx1 "<<wx1<<endl;
	cout<<"Rank "<<rank <<" Gx0 "<<gx0<<" Gx1 "<<gx1<<endl;



   

   }
   MPI_Barrier(MPI_COMM_WORLD);

   if (rank==pair[proc_intr].down){
	double det;
	double a = 1;
	double b = recv_intr[0];
	double d = 1;
	double c = W[0].front();
	det = a*d-b*c;
	
	double w0= recv_intr[1];
	double w1 = 0.0;
	double wx0=(d*w0-b*w1)/det;
	double wx1= (a*w1-c*w0)/det;
	
	double v0= 0.0;
    	double v1 = V[0].front();
	double vx0= (d*v0-b*v1)/det;
	double vx1= (a*v1-c*v0)/det;
	
	double G0 = recv_intr[2];
	double G1 = G[0];
	double gx0= (d*G0-b*G1)/det;
	double gx1= (a*G1-c*G0)/det;

	cout<<"Rank "<<rank <<" Vx0 "<<vx0<<" vx1 "<<vx1<<endl;
	cout<<"Rank "<<rank <<" Wx0 "<<wx0<<" Wx1 "<<wx1<<endl;
	cout<<"Rank "<<rank <<" Gx0 "<<gx0<<" Gx1 "<<gx1<<endl;



   }

   MPI_Barrier(MPI_COMM_WORLD);
}

}

    if (rank == 0){
	mainSPIKE(A_mat,B_vec,SIZE);
    	vector_subtract(x_vec, B_vec, SIZE);
    }

 
    MPI_Finalize();
}//end of main
