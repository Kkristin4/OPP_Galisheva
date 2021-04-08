#include <iostream>
#include <cmath>
#include <omp.h>
using namespace std;

int main() {
    int N = 10;
    static double  e = 1e-5;
    static double tay = 5e-5;
    double vector_norm = 0;
    double B_norm = 0;
    double* A = new double[N * N];
    double* B = new double[N];
    double* X = new double[N];
    double* result = new double[N];
    fill(A, A + N * N, 1);
    for (int i = 0; i < N; i++)
        A[i * N + i] = 2;
    fill (X, X + N, 0);
    fill (B, B + N, N + 1);
    B_norm = (N + 1) * (N + 1) * N;
    fill (result, result + N, 0.0f);
#pragma omp parallel
    do {
#pragma omp for
        for(int i = 0; i < N; i++){
            result[i] = 0.;
            for (int j = 0; j < N; j++){
                result[i]  += A[i * N + j] * X[j];
            }
        }
#pragma omp for
        for(int i = 0; i < N; i++){
            result[i] -= B[i];
        }
#pragma omp for
        for(int i = 0; i < N; i++){
            X[i] -= result[i] * tay;
        }
        vector_norm = 0;
#pragma omp for
        for (int i = 0; i < N; i++){
            vector_norm += result[i] * result[i];
        }
    } while (vector_norm / B_norm > e * e);
    for (int i = 0; i < N; ++i)
        cout << X[i] << " ";
    delete [] A;
    delete [] B;
    delete [] X;
    delete [] result;
    return 0;
}