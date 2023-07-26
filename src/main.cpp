#include <iostream>
#include <mpi.h>
#include <exception>
#include <cstdio>
#include <cmath>
#include "sample.h"

double F(int x, int y, int L) {

    if (x < 0 || x >= L || y < 0 || y >= L) {
        throw std::logic_error("grid coordinates are out of bounds");
    }

    return 0;
}

double G(int x, int y, int L) {

    if (x < 0 || x >= L || y < 0 || y >= L) {
        throw std::logic_error("grid coordinates are out of bounds");
    }

    auto delta = (2.0) * M_PI / L;
    if (y == 0 || y == L - 1) {
        // horizontal boundary
        return sin(delta * x);
    } else if (x == 0 || x == L - 1) {
        // vertical boundary
        return sin(delta * y);
    }

    throw std::logic_error("grid coordinates is not on a boundary");
}

double rms(int *U, int *V, int L) {

    double sum = 0;
    for (int x = 0; x < L; x++) {
        for (int y = 0; y < L; y++) {
            sum += pow(U[y * L + x] - V[y * L + x], 2);
        }
    }

    auto n = (double) (L * L);
    return sqrt((1.0 / n) * sum);
}

int main(int argc, char *argv[]) {

    // Let the domain of u,f, and g be 0 <= x <= 2 * pi and 0 <= y <= 2 * pi.
    // del2 u(x,y) = f(x, y) inside boundary
    // u(x, y) = g(x, y) on boundary
    // https://py-pde.readthedocs.io/en/latest/examples_gallery/laplace_eq_2d.html
    int L = 60;
    int W = 20;
    int H = 15;

    int num_process;
    int world_rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_process);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Create Cartesian Communicator
    int coords[2];
    int dim_size[2] = {L / W, L / H};
    int periods[2] = {0, 0};

    MPI_Comm world = MPI_COMM_WORLD;
    MPI_Comm cartesian;
    int cart_rank;

    MPI_Cart_create(world, 2, dim_size, periods, 1, &cartesian);
    MPI_Comm_rank(cartesian, &cart_rank);
    MPI_Cart_coords(cartesian, cart_rank, 2, coords);

    printf("world_rank=%d cart_rank=%d coords=(%d,%d)\n", world_rank, cart_rank, coords[0], coords[1]);

    // auto f = [L](int x, int y) { return F(x, y, L); };
    auto g = [L](int x, int y) { return G(x, y, L); };

    auto border = [L](int x, int y) {
        return x == 0 || x == L - 1 || y == 0 || y == L - 1;
    };

    // Create U and V
    auto U = (double *) malloc(sizeof(double) * (L * L));
    // auto V = (double *) malloc(sizeof(double) * GRID_SIZE);

    for (int x = 0; x < L; x++) {
        for (int y = 0; y < L; y++) {
            U[y * L + x] = border(x, y) ? g(x, y) : 0;
        }
    }

    auto x = [&](int i) {
        return coords[0] * W + (i - 1);
    };

    auto y = [&](int j) {
        return coords[1] * H + (j - 1);
    };

    // Create sub-grid for current process with halo
    auto u = (double *) malloc(sizeof(double) * (W + 2) * (H + 2));
    auto v = (double *) malloc(sizeof(double) * (W + 2) * (H + 2));

    // fill horizontal halo
    for (int i = 0; i <= (W + 1); i++) {
        u[0 * (W + 2) + i] = 0;
        u[(H + 1) * (W + 2) + i] = 0;
    }

    // fill vertical halo
    for (int j = 0; j <= (H + 1); j++) {
        u[j * (W + 2) + 0] = 0;
        u[j * (W + 2) + (W + 1)] = 0;
    }

    // fill
    for (int i = 1; i <= W; i++) {
        for (int j = 1; j <= H; j++) {
            u[j * (W + 2) + i] = border(x(i), y(j)) ? g(x(i), y(j)) : 0;
        }
    }

    for (int i = 0; i <= W + 1; i++) {
        for (int j = 0; j <= H + 1; j++) {
            v[j * (W + 2) + i] = 0;
        }
    }

    /* ----- */

    if (coords[0] == 0 && coords[1] == 0) {
        for (int j = H + 1; j >= 0; j--) {
            for (int i = 0; i < (W + 2); i++) {
                std::cout << u[j * (W + 2) + i] << " ";
            }
            std::cout << std::endl;
        }
    }
    /* ----- */

    // compute
    int step = 1;
    int maxstep = 5 * L;
    while (step <= maxstep) {

        // if (globalRMS < someValue) {break;}

        // swap halo

        // update v base on u

        // all reduce rms

        // copy u = v

        step++;
    }

    // create U with all zeroes except for the corresponding sub-grid of this process

    // All reduce U of world_rank = 0

    // write U to bitmap

    return 0;
}