#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__
#include<vector>
using namespace std;
inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }
//
//  saving parameters
//
//
// Processor exchange pairs
//
typedef struct{
int up;
int down;
}proc_pair;
//
//  timing routines
//
double read_timer( );

//
//  solver routines
//
void input_test(int size);
void print_tridiagonal(double *A, int size);
void printlocalmatrix(double *V, double *W, int size,int blocksize);
void print_Tridiag(double *A,int size);
void create_tridiagonal(double *A, int size);
void make_answer_to_solve(double *A, double *x, double *B,int size);
double max(double *a, double *b,int size);
void vector_subtract(double *x, double *guess, int height);
void twoxtwocalc(double *A, double *G, double *V, double *W,int type);
void finish(double *V, double *W, double *G,int size);
void mainSPIKE(double *A, double *B, int size);
void delete_alternative (vector<vector<double> > &vec, int start);
//
//  I/O routines
//

//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );

#endif
