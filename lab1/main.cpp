#include <iostream>
#include <cmath>
#include <mpi/mpi.h>
using namespace std;

int iflast(int size, int rank){
    if (size == rank + 1)
        return 1;
    return 0;
}

void multiplication_vector(double* matrix, double* vector, double* result, int N, int rank, int NumLines, int NormalLines, int * CountELEM, int * shift){

    for(int i = 0; i < NumLines; i++){
        result[i] = 0.;
        for (int j = 0; j < N; j++){
            result[i]  += matrix[i * N + j] * vector[j];
        }
    }
    MPI_Allgatherv(result, CountELEM[rank], MPI_DOUBLE, result, CountELEM, shift, MPI_DOUBLE, MPI_COMM_WORLD);


}

void multiplication_tay(double* vector, int N, int rank, int NormalLines){
    double tay = 5e-5;
    for(int i = 0; i < N; i++){
        vector[i] *= tay;
    }
}

void subtraction(double* deductible, double* subtractor, int N, int rank, int NumLines, int NormalLines){
    for(int i = 0; i < N; i++){
        deductible[i] -= subtractor[i];
    }
}

double vector_length(double* vector, int N){
    double result = 0.0f;
    for (int i = 0; i < N; i++){
        result += vector[i] * vector[i];
    }
    return result;
}

bool check(double* result, int N, double* B, int rank){
    double  e = 1e-5;
    double f = vector_length(result, N) / vector_length(B, N);
    if ( f < e * e){
        return true;}
    else {
        return false;
    }
}



int Number_lines(int N, int size){
    int lines = N / size;
    if (N % size != 0)
        lines++;
    return lines;
}

void subtractionNormal(double* X,double* result,int  N){
    for (int i = 0; i < N; i++){
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

    double* A= (double*)malloc(N * NumLines * sizeof(double*));
    for(int i = 0; i < N * NumLines; i++){
        A[i] = 1;
    }
    for( int i = 0; i < NumLines ; i++){
        A[ i * N + rank * NormalLines + i ] = 2;
    }

    double* B = (double*)malloc(N * sizeof(double*));
    double* X = (double*)malloc(N * sizeof(double*));
    double* result = (double*)malloc( N * sizeof(double*));
    int* CountELEM = (int*)malloc(size * sizeof(int*));
    int* shift = (int*)malloc(size * sizeof(int*));
    for(int i = 0; i < size; i++){
        shift[i] = i * NormalLines;
    }
    fill (X, X + N, 0);
    fill (B, B + N, N + 1);
    fill (result, result + N, 0.0f);

    fill (CountELEM, CountELEM + size, NormalLines);
    CountELEM[size - 1] = N - NormalLines * (size - 1);
    int k = 0;
    do {
        multiplication_vector(A, X, result, N, rank, NumLines, NormalLines, CountELEM, shift);
        subtraction(result, B, N, rank, NumLines, NormalLines);
        multiplication_tay(result, N, rank, NormalLines);
        subtractionNormal(X, result, N);

    }while(!check(result, N, B, rank));
    MPI_Barrier(MPI_COMM_WORLD);
    if (iflast(size, rank)) {
        for (int i = 0; i < N; ++i)
            cout << X[i] << " ";
    }

    MPI_Finalize();

    return 0;
}