#include <iostream>
#include <cmath>
#include <mpi/mpi.h>
using namespace std;
/*
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
 */
void  multiplication_vector(double * A, double  * B, double * C, int n1, int n2, int n3){
    double number;
    double* BT = new double[n2 * n3];
    for (int i = 0; i < n2; i++){
        for (int j = 0; j < n3; j++){
            BT[i + j * n2] = B[i * n3 + j];
        }
    }
    for (int i = 0; i < n3; i++){
        for (int j = 0; j < n1; j++){
            number = 0.;
            for (int k = 0; k < n2; k++){
                number += A[j * n2 + k] * BT[i * n3 + k];
            }
            C[j * n3 + i] = number;
        }
    }
    for (int i = 0; i < n1 * n3; i++){
        std::cout <<  C[i];
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);                         //презентация
    int dims[2] = {2, 2};
    int periods[2] = {0, 0};
    int coords[2];
    int reorder = 0;
    int size,rank,sizey,sizex,ranky,rankx;
    int prevy,prevx,nexty,nextx;
    MPI_Comm comm2d;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Dims_create(size,2,dims);
    sizey = dims[0]; sizex = dims[1];
    MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,reorder,&comm2d);
    MPI_Comm_rank(comm2d,&rank);
    MPI_Cart_get(comm2d,2,dims,periods,coords);
    ranky=coords[0]; rankx=coords[1];
    MPI_Cart_shift(comm2d,0,1,&prevy,&nexty);
    MPI_Cart_shift(comm2d,1,1,&prevx,&nextx);
    int colrank, linerank;
    MPI_Comm Col;
    MPI_Comm_split(MPI_COMM_WORLD, rankx, ranky, &Col);
    MPI_Comm_rank(Col,&colrank);
    MPI_Comm Line;
    MPI_Comm_split(MPI_COMM_WORLD, ranky, rankx, &Line);
    MPI_Comm_rank(Line,&linerank);

    int n1 = 2;
    int n2 = 4;
    int n3 = 6;
    if ((n1 % sizey != 0)||(n3 % sizex !=0)){
        MPI_Finalize();
        return 0;
    }
    int local_n1 = n1 / sizey;
    int local_n3 = n3 / sizex;
    double* local_A = new double [local_n1 * n2];
    fill(local_A, local_A + local_n1 * n2, 0);
    double* local_B = new double [local_n3 * n2];
    fill(local_B, local_B + local_n3 * n2, 0);
    double* local_C = new double [local_n1 * local_n3];
    fill(local_C, local_C + local_n1 * local_n3, 0);

    double* A;
    double* B;
    double* C;
    if (rank == 0) {
        A = new double[n1 * n2];
        for (int i = 0; i < n1 * n2; i++){
            if (i < 4) {A[i] = 1;}
            else {A[i] = 0;}
        }
        B = new double [n3 * n2];
        for (int i = 0; i < n3 * n2; i++){
            B[i] = 2;
            if (i < 3) B[i] = 1;
            if (i > 20) B[i] = 1;
        }
        C = new double [n1 * n3];
        fill(C, C + n1 * n3, 0.);

    }

    if (rankx == 0){
        MPI_Scatter(A, local_n1 * n2, MPI_DOUBLE, local_A, local_n1 * n2, MPI_DOUBLE, 0, Col);

    }
    MPI_Bcast(local_A, local_n1 * n2, MPI_DOUBLE, 0, Line);
    if (ranky == 0){
        MPI_Scatter(B, n2 * local_n3, MPI_DOUBLE, local_B, n2 * local_n3, MPI_DOUBLE, 0, Line);
    }
    MPI_Bcast(local_B, n2 * local_n3, MPI_DOUBLE, 0, Col);
if (rank == 0) {
    std::cout << local_A[0]<<'\n';
    multiplication_vector(local_A, local_B, local_C, local_n1, n2, local_n3);
}



  /*


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

  */
    MPI_Finalize();
    return 0;
}