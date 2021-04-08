#include <iostream>
#include <cmath>
#include <omp.h>
using namespace std;

void multiplication_vector(double* matrix, double* vector, double* result, int N){
#pragma omp parallel for
    for(int i = 0; i < N; i++){
        result[i] = 0.;
        for (int j = 0; j < N; j++){
            result[i]  += matrix[i * N + j] * vector[j];
        }
    }
}

void multiplication_tay(double* vector, int N){
    static double tay = 5e-5;
#pragma omp parallel for
    for(int i = 0; i < N; i++){
        vector[i] *= tay;
    }
}

void subtraction(double* deductible, double* subtractor, int N){
#pragma omp parallel for
    for(int i = 0; i < N; i++){
        deductible[i] -= subtractor[i];
    }
}

double vector_length( double* vector, int N){
    double result = 0.0f;
#pragma omp parallel for
    for (int i = 0; i < N; i++){
        result += vector[i] * vector[i];
    }
    return result;
}


bool check( double *result, int N, double *B) {
    static double  e = 1e-5;
    static double tay = 5e-5;
    double f = vector_length(result, N) / (vector_length(B, N) * tay);
    return f < e * e;
}



int main() {
    int N = 10;
    double* A = new double[N * N];
    double* B = new double[N];
    double* X = new double[N];
    double* result = new double[N];
    fill(A, A + N * N, 1);
    for (int i = 0; i < N; i++)
        A[i * N + i] = 2;
    fill (X, X + N, 0);
    fill (B, B + N, N + 1);
    fill (result, result + N, 0.0f);
    do {
        multiplication_vector(A, X, result, N);
        subtraction(result, B, N);
        multiplication_tay(result, N);
        subtraction(X, result, N);
    } while (!check(result, N, B));
    for (int i = 0; i < N; ++i)
        cout << X[i] << " ";

    delete [] A;
    delete [] B;
    delete [] X;
    delete [] result;
    return 0;
}