#include <iostream>
#include <cmath>
#include <cstdio>
#include <mpi/mpi.h>

struct Area {
    Area(int Nx, int Ny, int Nz) : Nx(Nx), Ny(Ny), Nz(Nz) {
        area = new double[(Nx +2) * (Ny + 2) * (Nz +2)];
        std::fill(area, area + (Nx +2) * (Ny + 2) * (Nz +2), 0);
    }
    double get(int x, int y, int z) {
        return area[x * Nz * Ny + y * Nz + z];
    }
    void set(int x,int  y,int  z, double value){
        area[x * Nz * Ny + y * Nz + z] = value;
    }
    void swap(Area &newarea){
        std::swap(area, newarea.area);
    }
private:
    int Nx, Ny, Nz;
    double* area;
};

double next_phi(Area &area, int i,int  j,int k, int rank, int size){
    const double Dx = 2, Dy = 2, Dz = 2;
    const int Nx = 8 / size, Ny = 8, Nz = 8;
    const double hx = Dx / (Nx - 1), hy = Dy / (Ny - 1), hz = Dz / (Nz - 1);
    const int a = 100000;

    const double x0 = -1 + Dx / size * rank, y0 = -1, z0 = -1;
    double hx_2 = hx * hx, hy_2 = hy *hy, hz_2 = hz * hz;
    double cof = 1/(2 / hx_2 + 2/ hy_2 + 2/ hz_2 + a);
    double first = (area.get(i+1, j, k) + area.get(i-1, j, k)) / hx_2;
    double second = (area.get(i, j+1, k) + area.get(i, j-1, k)) / hy_2;
    double third = (area.get(i, j, k+1) + area.get(i, j, k-1))  / hz_2;
    double x = x0 + i *hx, y = y0 + j * hy, z = z0 + k * hz;
    double ro = 6 - a* (x *x + y* y + z *z);
    return cof * (first + second +third - ro);
}

void check_max(double & max, double first, double second){
    double diff = std::abs(first - second);
    if (max < diff)
        max = diff;

}

int main(int argc, char *argv[]) {
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    const double e = 1e-10;
    const int phi0 = 0;

    const int Nx = 8 / size, Ny = 8, Nz = 8;

  /*  double * area = new double [Nx * Ny * Nz];
    std::fill(area, area + Nz *Nx* Ny, 0);
    double *newarea = new double[Nx * Ny * Nz];*/
    Area area(Nx, Ny, Nz);
    Area newarea(Nx, Ny, Nx);
    double max = 0;
    double * buf  = new double[Nz * Ny];
    MPI_Request request;
   do {
       
       max = 0;
       if (rank != size - 1) {
           MPI_Isend(buf, cout, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &request);
       }
       for (int i = 1; i < Nx +1; i++) {
           for (int j = 1; j < Ny+1; j++) {
               for (int k = 1; k < Nz+1; k++) {
                   newarea.set(i, j, k, next_phi(area, i, j, k, rank, size));
                   check_max(max, newarea.get(i, j, k), area.get(i, j, k));
               }
           }
       }
       area.swap(newarea);
       MPI_Status status;
       MPI_Wait(&request, &status);
   } while (max > e);


    MPI_Finalize();

    return 0;
}