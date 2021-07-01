#include <iostream>
#include <cmath>
#include <mpi/mpi.h>
using namespace std;

void  multiplication_vector(double * A, double  * B, double * C, int n1, int n2, int n3){
    double number;
    for (int i = 0; i < n3; i++) {
        for (int j = 0; j < n1; j++) {
            number = 0.;
            for (int k = 0; k < n2; k++) {
                number += A[j * n2 + k] * B[i * n2 + k];
            }
            C[j * n3 + i] = number;
        }
    }
}

void CT_v_C(double *CT, double * C, int n1, int n3, int local_n1, int local_n3, int sizex){
   int rank = 0;
   for ( int y  = 0; y < n1; y++){
       for (int x = 0; x < n3; x++){
           rank  = (y / local_n1) * sizex + x / local_n3;
           C[y*n1 + x] = CT[rank * n3 + (x % local_n3) + (y % local_n1) * local_n3];
       }
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

    int n1 = 4;
    int n2 = 4;
    int n3 = 4;
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
    double * CT;
    double* BT;
    int*  count;
    int* displs;
    count = new int[size];
    fill(count, count + size, local_n1 * local_n3);
    displs = new int[size];
    for(int i = 0; i < size; i++){
        displs[i] = i * local_n3 * local_n1;
    }
    if (rank == 0) {
        A = new double[n1 * n2];
        for (int i = 0; i < n1 * n2; i++){
            A[i] = i;
        }
        B = new double [n3 * n2];
        for (int i = 0; i < n3 * n2; i++){
            B[i] = i;
        }
        BT = new double[n2 * n3];
        for (int i = 0; i < n2; i++){
            for (int j = 0; j < n3; j++){
                BT[i + j * n2] = B[i * n3 + j];
            }
        }
        C = new double [n1 * n3];
        fill(C, C + n1 * n3, 0.);
        CT = new double [n1 * n3];
        fill(C, C + n1 * n3, 0.);
    }

    if (rankx == 0){
        MPI_Scatter(A, local_n1 * n2, MPI_DOUBLE, local_A, local_n1 * n2, MPI_DOUBLE, 0, Col);
    }
    MPI_Bcast(local_A, local_n1 * n2, MPI_DOUBLE, 0, Line);

    if (ranky == 0){
        MPI_Scatter(BT, n2 * local_n3, MPI_DOUBLE, local_B, n2 * local_n3, MPI_DOUBLE, 0, Line);
    }
    MPI_Bcast(local_B, n2 * local_n3, MPI_DOUBLE, 0, Col);
    multiplication_vector(local_A, local_B, local_C, local_n1, n2, local_n3);
    MPI_Gatherv(local_C, local_n1 * local_n3, MPI_DOUBLE, CT,count,displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (rank == 0){
      CT_v_C(CT, C, n1, n3, local_n1, local_n3, sizex);
      for (int i = 0; i < n1 * n3; i++){
          std::cout << C[i] << " ";
      }
    }

    MPI_Finalize();
    return 0;
}