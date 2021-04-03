#include <iostream>
#include <cmath>
#include <mpi/mpi.h>
using namespace std;

bool iflast(int size, int rank){
    return size == rank + 1;
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
    static double tay = 5e-5;
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
    static double  e = 1e-5;
    double f = vector_length(result, N) / vector_length(B, N);
    return f < e * e;
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

    double* A = new double[N * NumLines];
    for(int i = 0; i < N * NumLines; i++){
        A[i] = 1;
    }
    for( int i = 0; i < NumLines ; i++){
        A[ i * N + rank * NormalLines + i ] = 2;
    }

    double* B = new double[N];
    double* X = new double[N];
    double* result = new double[N];
    int* CountELEM = new int[size];
    int* shift = new int[size];
    for(int i = 0; i < size; i++){
        shift[i] = i * NormalLines;
    }
    fill (X, X + N, 0);
    fill (B, B + N, N + 1);
    fill (result, result + N, 0.0f);

    fill (CountELEM, CountELEM + size, NormalLines);
    CountELEM[size - 1] = N - NormalLines * (size - 1);

    do {
        multiplication_vector(A, X, result, N, rank, NumLines, NormalLines, CountELEM, shift);
        subtraction(result, B, N, rank, NumLines, NormalLines);
        multiplication_tay(result, N, rank, NormalLines);
        subtractionNormal(X, result, N);

    }while(!check(result, N, B, rank));
    if (iflast(size, rank)) {
        for (int i = 0; i < N; ++i)
            cout << X[i] << " ";
    }
    delete [] A;
    delete [] B;
    delete [] X;
    delete [] result;
    delete [] CountELEM;
    delete [] shift;
    MPI_Finalize();

    return 0;
}