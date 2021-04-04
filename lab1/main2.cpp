#include <iostream>
#include <cmath>
#include <mpi/mpi.h>
using namespace std;

bool iflast(int size, int rank){
    return size == rank + 1;
}

void multiplication_vector(double* matrix, double* vector, double* result, int N, int NumLines){
    for (int i = 0; i < N; i++){
        result[i] = 0.;
        for (int j = 0; j < NumLines; j++){
            result[i] += matrix[j * N + i] * vector[j];
        }
    }
    for (int i = 0; i < N; i++) {
        MPI_Allreduce(&result[i], &result[i], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
}

void multiplication_tay(double* vector,  int NumLines){
    static double tay = 5e-5;
    for(int i = 0; i < NumLines; i++){
        vector[i] *= tay;
    }
}

void subtraction(double* deductible, double* subtractor, int NumLines){
    for(int i = 0; i < NumLines; i++){
        deductible[i] -= subtractor[i];
    }
}

double vector_length( double* vector, int N){
    double result = 0.0f;
    for (int i = 0; i < N; i++){
        result += vector[i] * vector[i];
    }
    return result;
}


bool check( double *result, int NumLines, double *B, double norma_B, double norma) {
    static double  e = 1e-5;
    norma = vector_length(result, NumLines);
    norma_B = vector_length(B, NumLines);
    MPI_Allreduce(&norma, &norma, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&norma_B, &norma_B, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double f = norma / norma_B;
    return f < e * e;
}


int Number_lines(int N, int size){
    int lines = N / size;
    if (N % size != 0)
        lines++;
    return lines;
}

void subtractionNormal(double* X,double* result,int  NumLines){
    for (int i = 0; i < NumLines; i++){
        X[i] -= result[i];
    }
}

int main(int argc, char *argv[]) {
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);



    int N = 10;
    int NormalLines = Number_lines(N, size);
    int NumLines = NormalLines;

    if (iflast(size, rank))
        NumLines = (N - rank * NumLines);


   double* A = new double[N * NumLines];
    for(int i = 0; i < N * NumLines; i++){
        A[i] = 1;
    }
    for( int i = 0; i < NumLines ; i++){
        A[ i * N + rank * NormalLines + i ] = 2;
    }

    double* B = new double[NumLines];
    double* X = new double[NumLines];
    double* result = new double[N];
    double norma_B = 0;
    double norma = 0;

    fill (X, X + NumLines, 0);
    fill (B, B + NumLines, N + 1);
    fill (result, result + N, 0.0f);

    do {
        multiplication_vector(A, X, result, N,  NumLines);
        subtraction(result, B, NumLines);
        multiplication_tay(result, NumLines);
        subtractionNormal(X, result, NumLines);

    }while(!check(result, NumLines, B, norma_B, norma));

   for (int i = 0; i < size; i ++){
       if (rank == i){
           for (int j = 0; j < NumLines; j++)
               cout << X[j] << " ";
       }
   }
    delete [] A;
    delete [] B;
    delete [] X;
    delete [] result;
    MPI_Finalize();
    return 0;
}