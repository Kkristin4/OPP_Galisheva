#include <iostream>
#include <cmath>
#include <cstdio>
#include <mpi/mpi.h>

struct Area {
    Area(int Nx, int Ny, int Nz) : Nx(Nx), Ny(Ny), Nz(Nz) {
        area = new double[(Nx + 2) * (Ny + 2) * (Nz + 2)];
        std::fill(area, area + (Nx + 2) * (Ny + 2) * (Nz + 2), 0);
    }

    double get(int x, int y, int z) {
        return area[x * (Nz + 2) * (Ny +2) + y * (Nz + 2) + z];
    }

    void set(int x, int y, int z, double value) {
        area[x * (Nz +2) * (Ny + 2) + y * (Nz + 2)+ z] = value;
    }

    void swap(Area &newarea) {
        std::swap(area, newarea.area);
    }

    double *link(int x) {
        return area + (x * (Nz + 2) * (Ny + 2));
    }

    ~Area() {
        delete[] area;
    }

private:
    int Nx, Ny, Nz;
    double *area;
};

double next_phi(Area &area, int i, int j, int k, int rank, int size) {
    const double Dx = 2, Dy = 2, Dz = 2;
    const int Nx = 8 / size, Ny = 8, Nz = 8;
    const double hx = Dx / (Nx - 1), hy = Dy / (Ny - 1), hz = Dz / (Nz - 1);
    const int a = 100000;

    const double x0 = -1 + Dx / size * rank, y0 = -1, z0 = -1;
    double hx_2 = hx * hx, hy_2 = hy * hy, hz_2 = hz * hz;
    double cof = 1 / (2 / hx_2 + 2 / hy_2 + 2 / hz_2 + a);
    double first = (area.get(i + 1, j, k) + area.get(i - 1, j, k)) / hx_2;
    double second = (area.get(i, j + 1, k) + area.get(i, j - 1, k)) / hy_2;
    double third = (area.get(i, j, k + 1) + area.get(i, j, k - 1)) / hz_2;
    double x = x0 + i * hx, y = y0 + j * hy, z = z0 + k * hz;
    double ro = 6 - a * (x * x + y * y + z * z);
    return cof * (first + second + third - ro);
}

void check_max(double &max, double first, double second) {
    double diff = std::abs(first - second);
    if (max < diff)
        max = diff;

}

int main(int argc, char *argv[]) {
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::cout << std::endl;

    const double e = 1e-14;
    const int phi0 = 0;

    const int Nx = 8 / size, Ny = 8, Nz = 8;
    Area area(Nx, Ny, Nz);
    Area newarea(Nx, Ny, Nz);
    double max = 0;
    MPI_Request hight;
    MPI_Request low;
    MPI_Request h;
    MPI_Request l;
    int cout = (Ny + 2) * (Nz + 2);
    do {
        max = 0;
        if (rank != size - 1) {
            for (int j = 2; j < Ny; j++) {
                for (int k = 2; k < Nz; k++) {
                    newarea.set(Nx, j, k, next_phi(area, Nx , j, k, rank, size));
                    check_max(max, newarea.get(Nx, j, k), area.get(Nx, j, k));
                }
            }
            MPI_Isend(area.link(Nx), cout, MPI_DOUBLE, rank + 1, 123, MPI_COMM_WORLD, &hight);
            MPI_Irecv(newarea.link(Nx + 1), cout, MPI_DOUBLE, rank + 1, 1234, MPI_COMM_WORLD, &l);
        }
        if (rank != 0) {
            for (int j = 2; j < Ny; j++) {
                for (int k = 2; k < Nz; k++) {
                    newarea.set(1, j, k, next_phi(area, 1, j, k, rank, size));
                    check_max(max, newarea.get(1, j, k), area.get(1, j, k));

                }
            }
            MPI_Isend(area.link(1), cout, MPI_DOUBLE, rank - 1, 1234, MPI_COMM_WORLD, &low);
            MPI_Irecv(newarea.link(0), cout, MPI_DOUBLE, rank - 1, 123, MPI_COMM_WORLD, &h);

        }

        for (int i = 2; i < Nx ; i++) {
            for (int j = 1; j < Ny + 1; j++) {
                for (int k = 1; k < Nz + 1; k++) {
                    newarea.set(i, j, k, next_phi(area, i, j, k, rank, size));
                    check_max(max, newarea.get(i, j, k), area.get(i, j, k));
                }
            }
        }
        if (rank != size - 1) {
            MPI_Status status;
            MPI_Wait(&hight, &status);

        }
        if (rank != 0) {
            MPI_Status status;
            MPI_Wait(&low, &status);

        }
        if (rank != 0) {
            MPI_Status status;
            MPI_Wait(&h, &status);

        }
        if (rank != size - 1) {
            MPI_Status status;
            MPI_Wait(&l, &status);

        }

        MPI_Allreduce(&max, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    // if (rank == 0){   std::cout  << max << std::endl;}
        area.swap(newarea);
    } while (max > e );
   // std::cout << "out " << rank << "  "  << max<< '\n' << std::endl;

    MPI_Finalize();

    return 0;
}