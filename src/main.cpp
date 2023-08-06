#include <iostream>
#include <mpi.h>
#include <cmath>
#include <cstring>
#include <cfloat>
#include "bitmap.h"

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
    if (y == L - 1) {
        // horizontal boundary
        return sin(delta * x);
    } else if (y == 0) {
        return sin(delta * x);
    } else if (x == 0) {
        // vertical boundary
        return sin(delta * y);
    } else if (x == L - 1) {
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

    int root = 0;
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

    auto f = [L](int x, int y) { return F(x, y, L); };
    auto g = [L](int x, int y) { return G(x, y, L); };

    auto border = [L](int x, int y) {
        return x == 0 || x == L - 1 || y == 0 || y == L - 1;
    };

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


    auto rank = [&](int cartX, int cartY) {

        if (cartX < 0 || cartX >= dim_size[0] ||
            cartY < 0 || cartY >= dim_size[1]) {
            return MPI_PROC_NULL;
        }

        int rank;
        int coords[2] = {cartX, cartY};
        MPI_Cart_rank(cartesian, coords, &rank);
        return rank;
    };

    MPI_Datatype HORIZONTAL_HALO;
    MPI_Datatype VERTICAL_HALO;

    // count, blocklength, stride
    MPI_Type_vector(1, W, W, MPI_DOUBLE, &HORIZONTAL_HALO);
    MPI_Type_vector(H, 1, W + 2, MPI_DOUBLE, &VERTICAL_HALO);
    MPI_Type_commit(&HORIZONTAL_HALO);
    MPI_Type_commit(&VERTICAL_HALO);

    auto left = rank(coords[0] - 1, coords[1]);
    auto right = rank(coords[0] + 1, coords[1]);
    auto up = rank(coords[0], coords[1] + 1);
    auto down = rank(coords[0], coords[1] - 1);

    int tag = 1;

    // Receiving buffers for vertical halos
    auto bufferLeft = (double *) malloc(H * sizeof(double));
    auto bufferRight = (double *) malloc(H * sizeof(double));

    /* -----


    if (coords[0] == 0 && coords[1] == 0) {
        for (int j = H + 1; j >= 0; j--) {
            for (int i = 0; i < (W + 2); i++) {
                std::cout << u[j * (W + 2) + i] << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "------------" << std::endl;
        std::cout << "------------" << std::endl;
    }
    ----- */

    // compute
    int step = 1;
    int maxstep = 40 * L;
    while (step <= maxstep) {

        // if (globalRMS < someValue) {break;}

        // swap halo
        MPI_Status status;
        MPI_Request rightSend, leftSend, upSend, downSend;
        MPI_Issend(&u[1 * (W + 2) + 1], 1, VERTICAL_HALO, left, tag, cartesian, &leftSend);
        MPI_Issend(&u[1 * (W + 2) + W], 1, VERTICAL_HALO, right, tag, cartesian, &rightSend);
        MPI_Issend(&u[1 * (W + 2) + 1], 1, HORIZONTAL_HALO, down, tag, cartesian, &downSend);
        MPI_Issend(&u[H * (W + 2) + 1], 1, HORIZONTAL_HALO, up, tag, cartesian, &upSend);

        // Receive horizontal halos directly onto the current sub-grid
        MPI_Recv(&u[(H + 1) * (W + 2) + 1], W, MPI_DOUBLE, up, tag, cartesian, &status);
        MPI_Recv(&u[0 * (W + 2) + 1], W, MPI_DOUBLE, down, tag, cartesian, &status);

        // Receive the vertical halos to buffer arrays
        MPI_Recv(&bufferLeft[0], H, MPI_DOUBLE, left, tag, cartesian, &status);
        MPI_Recv(&bufferRight[0], H, MPI_DOUBLE, right, tag, cartesian, &status);

        if (left != MPI_PROC_NULL) {
            for (int j = 1; j <= H; j++) {
                u[j * (W + 2) + 0] = bufferLeft[j - 1];
            }
        }

        if (right != MPI_PROC_NULL) {
            for (int j = 1; j <= H; j++) {
                u[j * (W + 2) + W + 1] = bufferRight[j - 1];
            }
        }

        auto h2 = pow(1.0 / (L + 1.0), 2.0);

        // update v base on u
        for (int i = 1; i <= W; i++) {
            for (int j = 1; j <= H; j++) {

                if (border(x(i), y(j))) {
                    v[j * (W + 2) + i] = u[j * (W + 2) + i];
                } else {

                    v[j * (W + 2) + i] = (0.25) * (u[j * (W + 2) + (i - 1)] +
                                                   u[j * (W + 2) + (i + 1)] +
                                                   u[(j - 1) * (W + 2) + i] +
                                                   u[(j + 1) * (W + 2) + i] -
                                                   h2 * f(x(i), y(j)));
                }

            }
        }

        // all reduce rms

        /* Copies contents of v to u */
        std::memcpy(u, v, ((W + 2) * (H + 2)) * sizeof(double));

        step++;
    }

    /* -----

    if (coords[0] == 0 && coords[1] == 0) {
        for (int j = H + 1; j >= 0; j--) {
            for (int i = 0; i < (W + 2); i++) {
                std::cout << u[j * (W + 2) + i] << " ";
            }
            std::cout << std::endl;
        }
    }
    ----- */

    // Create U
    auto U = (double *) malloc(sizeof(double) * (L * L));
    auto V = (double *) malloc(sizeof(double) * (L * L));
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            U[j * L + i] = 0;
            V[j * L + i] = 0;
        }
    }


    // Paste u on the correct block of U
    for (int i = 1; i <= W; i++) {
        for (int j = 1; j <= H; j++) {
            U[y(j) * L + x(i)] = u[j * (W + 2) + i];
        }
    }


    // All reduce U of world_rank = 0
    // send buffer, receive buffer, count, datatype, op, root, comm
    MPI_Reduce(&U[0], &V[0], L * L, MPI_DOUBLE, MPI_SUM, root, cartesian);

    // write U to bitmap
    if (world_rank == root) {

        Bitmap bitmap(L, L);
        auto maximum = -DBL_MAX;
        auto minimum = DBL_MAX;
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                maximum = std::max(maximum, V[j * L + i]);
                minimum = std::min(minimum, V[j * L + i]);
            }
        }

        std::cout << minimum << std::endl;
        std::cout << maximum << std::endl;

        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                auto value = std::round(std::abs((V[j * L + i] - minimum) / (maximum - minimum)) * 255.0);
                bitmap.setPixel(i, j, value, value, value);
            }
        }

        std::cout << "---" << std::endl;

        std::cout << bitmap.write("test.bmp") << std::endl;
    }

    MPI_Type_free(&VERTICAL_HALO);
    MPI_Type_free(&HORIZONTAL_HALO);
    MPI_Finalize();
    free(V);
    free(U);
    free(u);
    free(v);

    return 0;
}