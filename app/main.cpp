#include <iostream>
#include <mpi.h>
#include <cmath>
#include <cstring>
#include <cfloat>
#include <functional>
#include "parameters.h"
#include "bitmap.h"
#include "util.h"

double squares(double *u, double *v, int W, int H) {

    double sum2 = 0;
    for (int i = 1; i <= W; i++) {
        for (int j = 1; j <= H; j++) {
            sum2 += pow(u[j * (W + 2) + i] - v[j * (W + 2) + i], 2);
        }
    }

    //auto n = (double) (L * L);
    //return sqrt((1.0 / n) * sum2);
    return sum2;
}

int main(int argc, char *argv[]) {

    // Let the domain of A,f, and g be 0 <= x <= 2 * pi and 0 <= y <= 2 * pi.
    // del2 A(x,y) = f(x, y) inside boundary
    // A(x, y) = g(x, y) on boundary
    // https://py-pde.readthedocs.io/en/latest/examples_gallery/laplace_eq_2d.html

    int root = 0;
    int num_process;
    int world_rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_process);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Create Cartesian Communicator
    int cartesian_coords[2];
    int cartesian_dimension[2] = {L / W, L / H};
    int periods[2] = {0, 0};

    MPI_Comm world = MPI_COMM_WORLD;
    MPI_Comm cartesian;
    int cartesian_rank;

    MPI_Cart_create(world, 2, cartesian_dimension, periods, 1, &cartesian);
    MPI_Comm_rank(cartesian, &cartesian_rank);
    MPI_Cart_coords(cartesian, cartesian_rank, 2, cartesian_coords);

    MPI_Datatype HORIZONTAL_HALO;
    MPI_Datatype VERTICAL_HALO;

    // count, blocklength, stride
    MPI_Type_vector(1, W, W, MPI_DOUBLE, &HORIZONTAL_HALO);
    MPI_Type_vector(H, 1, W + 2, MPI_DOUBLE, &VERTICAL_HALO);
    MPI_Type_commit(&HORIZONTAL_HALO);
    MPI_Type_commit(&VERTICAL_HALO);

    int cartesian_x = cartesian_coords[0];
    int cartesian_y = cartesian_coords[1];
    printf("world_rank=%d cartesian_rank=%d coords=(%d,%d)\n", world_rank, cartesian_rank, cartesian_x, cartesian_y);

    auto x = [&](int i) { return cartesian_x * W + (i - 1); };
    auto y = [&](int j) { return cartesian_y * H + (j - 1); };
    auto isBorder = [&](int x, int y) { return x == 0 || x == L - 1 || y == 0 || y == L - 1; };

    auto rank = std::bind(coordToRank, std::placeholders::_1, std::placeholders::_2, cartesian_dimension, cartesian);
    auto left = rank(cartesian_x - 1, cartesian_y);
    auto right = rank(cartesian_x + 1, cartesian_y);
    auto up = rank(cartesian_x, cartesian_y + 1);
    auto down = rank(cartesian_x, cartesian_y - 1);

    // Create sub-globalGrid for current process with halo
    // The solution is computed alternately in the array A and B.
    // The iteration is terminated when the difference between the successive approximations to the solution
    // is less than ...
    auto smallGridA = (double *) malloc(sizeof(double) * (W + 2) * (H + 2));
    auto smallGridB = (double *) malloc(sizeof(double) * (W + 2) * (H + 2));

    auto A = [&](int i, int j) -> double & { return smallGridA[j * (W + 2) + i]; };
    auto B = [&](int i, int j) -> double & { return smallGridB[j * (W + 2) + i]; };

    // fill
    for (int i = 1; i <= W; i++) {
        for (int j = 1; j <= H; j++) {
            A(i, j) = isBorder(x(i), y(j)) ? g(x(i), y(j)) : 0;
            B(i, j) = A(i, j) + 100 * TARGET_RMS;
        }
    }

    int tag = 1;

    // Receiving buffers for vertical halos
    auto bufferLeft = (double *) malloc(H * sizeof(double));
    auto bufferRight = (double *) malloc(H * sizeof(double));

    // compute
    int step = 1;
    while (step <= MAX_ITERATION) {

        // if (globalRMS < someValue) {break;}

        // swap halo
        MPI_Request rightSend, leftSend, upSend, downSend;
        MPI_Issend(&A(1, 1), 1, VERTICAL_HALO, left, tag, cartesian, &leftSend);
        MPI_Issend(&A(W, 1), 1, VERTICAL_HALO, right, tag, cartesian, &rightSend);
        MPI_Issend(&A(1, 1), 1, HORIZONTAL_HALO, down, tag, cartesian, &downSend);
        MPI_Issend(&A(1, H), 1, HORIZONTAL_HALO, up, tag, cartesian, &upSend);

        // Receive horizontal halos directly onto the current sub-globalGrid
        MPI_Recv(&A(1, H + 1), W, MPI_DOUBLE, up, tag, cartesian, MPI_STATUS_IGNORE);
        MPI_Recv(&A(1, 0), W, MPI_DOUBLE, down, tag, cartesian, MPI_STATUS_IGNORE);

        // Receive the vertical halos to buffer arrays
        MPI_Recv(&bufferLeft[0], H, MPI_DOUBLE, left, tag, cartesian, MPI_STATUS_IGNORE);
        MPI_Recv(&bufferRight[0], H, MPI_DOUBLE, right, tag, cartesian, MPI_STATUS_IGNORE);

        if (left != MPI_PROC_NULL) {
            for (int j = 1; j <= H; j++) {
                A(0, j) = bufferLeft[j - 1];
            }
        }

        if (right != MPI_PROC_NULL) {
            for (int j = 1; j <= H; j++) {
                A(W + 1, j) = bufferRight[j - 1];
            }
        }

        auto h2 = pow(1.0 / (L + 1.0), 2.0);

        // update B base on A
        for (int i = 1; i <= W; i++) {
            for (int j = 1; j <= H; j++) {

                if (isBorder(x(i), y(j))) {
                    B(i, j) = A(i, j);
                } else {

                    B(i, j) = (0.25) * (A(i - 1, j) +
                                        A(i + 1, j) +
                                        A(i, j - 1) +
                                        A(i, j + 1) -
                                        h2 * f(x(i), y(j)));
                }

            }
        }

        // all reduce rms
        // &sendbuf, &recvbuf, count, datatype, operation, comm
        double globalSquares = 0;
        double localSquares = squares(smallGridA, smallGridB, W, H);
        MPI_Allreduce(&localSquares, &globalSquares, 1, MPI_DOUBLE, MPI_SUM, cartesian);

        auto n = (double) (L * L);
        auto rms = sqrt((1.0 / n) * globalSquares);

        /* Copies contents of B to A */
        std::memcpy(smallGridA, smallGridB, ((W + 2) * (H + 2)) * sizeof(double));

        if (rms < TARGET_RMS) {
            break;
        }

        step++;
    }

    auto localGrid = (double *) malloc(sizeof(double) * (L * L));
    auto globalGrid = (double *) malloc(sizeof(double) * (L * L));
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            localGrid[j * L + i] = 0;
            globalGrid[j * L + i] = 0;
        }
    }

    // Paste A on the correct block of localGrid
    for (int i = 1; i <= W; i++) {
        for (int j = 1; j <= H; j++) {
            localGrid[y(j) * L + x(i)] = A(i, j);
        }
    }

    // All reduce localGrid of world_rank = 0
    // send buffer, receive buffer, count, datatype, op, root, comm
    MPI_Reduce(&localGrid[0], &globalGrid[0], L * L, MPI_DOUBLE, MPI_SUM, root, cartesian);

    // write localGrid to bitmap
    if (world_rank == root) {

        Bitmap bitmap(L, L);
        auto maximum = -DBL_MAX;
        auto minimum = DBL_MAX;
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                maximum = std::max(maximum, globalGrid[j * L + i]);
                minimum = std::min(minimum, globalGrid[j * L + i]);
            }
        }

        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                auto pixel = std::round(std::abs((globalGrid[j * L + i] - minimum) / (maximum - minimum)) * 255.0);
                bitmap.setPixel(i, j, pixel, pixel, pixel);
            }
        }

        std::cout << bitmap.write(OUTPUT_BMP) << std::endl;
    }

    MPI_Type_free(&VERTICAL_HALO);
    MPI_Type_free(&HORIZONTAL_HALO);
    MPI_Finalize();
    free(globalGrid);
    free(localGrid);
    free(smallGridA);
    free(smallGridB);

    return 0;
}