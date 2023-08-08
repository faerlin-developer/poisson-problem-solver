#include <mpi.h>
#include <cmath>
#include <cstring>
#include <cfloat>
#include "main.h"
#include "parameters.h"
#include "bitmap.h"
#include "csv.h"


int main(int argc, char *argv[]) {

    // ---- BEGIN INITIALIZE MPI ENVIRONMENT ----
    int num_process;
    int world_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_process);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm world = MPI_COMM_WORLD;

    // Create cartesian communicator with non-periodic boundaries and a dimension derived from global L, W, and H.
    int cartesian_rank;
    int cartesian_coords[2];
    int cartesian_periodicity[2] = {0, 0};
    int cartesian_dimension[2] = {L / W, L / H};

    MPI_Comm cartesian;
    MPI_Cart_create(world, 2, cartesian_dimension, cartesian_periodicity, 1, &cartesian);
    MPI_Comm_rank(cartesian, &cartesian_rank);
    MPI_Cart_coords(cartesian, cartesian_rank, 2, cartesian_coords);

    //
    int cartesian_x = cartesian_coords[0];
    int cartesian_y = cartesian_coords[1];
    auto left = coordToRank(cartesian_x - 1, cartesian_y, cartesian_dimension, cartesian);
    auto right = coordToRank(cartesian_x + 1, cartesian_y, cartesian_dimension, cartesian);
    auto up = coordToRank(cartesian_x, cartesian_y + 1, cartesian_dimension, cartesian);
    auto down = coordToRank(cartesian_x, cartesian_y - 1, cartesian_dimension, cartesian);

    // Create vector types for ghost points
    MPI_Datatype HORIZONTAL_HALO;
    MPI_Datatype VERTICAL_HALO;
    MPI_Type_vector(1, W, W, MPI_DOUBLE, &HORIZONTAL_HALO);    // count=1, blocklength=W, stride=W
    MPI_Type_vector(H, 1, W + 2, MPI_DOUBLE, &VERTICAL_HALO);  // count=H, blocklength=1, stride=W+2
    MPI_Type_commit(&HORIZONTAL_HALO);
    MPI_Type_commit(&VERTICAL_HALO);

    printf("Running process with rank %d and cartesian coords (%d,%d)\n", cartesian_rank, cartesian_x, cartesian_y);
    // ---- END INITIALIZE MPI ENVIRONMENT ----

    // Create sub-globalGrid for current process with halo
    // The solution is computed alternately in the array A and B.
    // The iteration is terminated when the residuals between the successive approximations to the solution
    // is less than ...
    auto smallGridA = (double *) malloc(sizeof(double) * (W + 2) * (H + 2));
    auto smallGridB = (double *) malloc(sizeof(double) * (W + 2) * (H + 2));
    auto bufferLeft = (double *) malloc(H * sizeof(double));
    auto bufferRight = (double *) malloc(H * sizeof(double));

    auto A = [&](int i, int j) -> double & { return smallGridA[j * (W + 2) + i]; };
    auto B = [&](int i, int j) -> double & { return smallGridB[j * (W + 2) + i]; };

    auto x = [&](int i) { return cartesian_x * W + (i - 1); };
    auto y = [&](int j) { return cartesian_y * H + (j - 1); };
    auto isBorder = [&](int x, int y) { return x == 0 || x == L - 1 || y == 0 || y == L - 1; };

    // fill
    for (int i = 1; i <= W; i++) {
        for (int j = 1; j <= H; j++) {
            A(i, j) = isBorder(x(i), y(j)) ? g(x(i), y(j)) : 0;
            B(i, j) = A(i, j) + 10 * TARGET_RMSE;
        }
    }

    // ---- BEGIN JACOBI METHOD ----
    // Apply Jacobi Method until target RMS is reached or maximum of iterations is exceeded
    auto step = 1;
    auto rmse = DBL_MAX;
    while (step <= MAX_ITERATION && rmse > TARGET_RMSE) {

        // ---- BEGIN SWAP HALO ----
        MPI_Request rightSend, leftSend, upSend, downSend;
        MPI_Issend(&A(1, 1), 1, VERTICAL_HALO, left, 1, cartesian, &leftSend);
        MPI_Issend(&A(W, 1), 1, VERTICAL_HALO, right, 1, cartesian, &rightSend);
        MPI_Issend(&A(1, 1), 1, HORIZONTAL_HALO, down, 1, cartesian, &downSend);
        MPI_Issend(&A(1, H), 1, HORIZONTAL_HALO, up, 1, cartesian, &upSend);

        // Receive horizontal halos and store them directly onto the small grid
        MPI_Recv(&A(1, H + 1), W, MPI_DOUBLE, up, MPI_ANY_TAG, cartesian, MPI_STATUS_IGNORE);
        MPI_Recv(&A(1, 0), W, MPI_DOUBLE, down, MPI_ANY_TAG, cartesian, MPI_STATUS_IGNORE);

        // Receive the vertical halos to buffer arrays
        MPI_Recv(&bufferLeft[0], H, MPI_DOUBLE, left, MPI_ANY_TAG, cartesian, MPI_STATUS_IGNORE);
        MPI_Recv(&bufferRight[0], H, MPI_DOUBLE, right, MPI_ANY_TAG, cartesian, MPI_STATUS_IGNORE);

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
        // ---- END SWAP HALO ----

        // ---- BEGIN JACOBI UPDATE ----
        auto h2 = pow(1.0 / (L + 1.0), 2.0);
        for (int i = 1; i <= W; i++) {
            for (int j = 1; j <= H; j++) {
                //
                if (isBorder(x(i), y(j))) {
                    B(i, j) = A(i, j);
                } else {
                    B(i, j) = 0.25 * (A(i - 1, j) + A(i + 1, j) + A(i, j - 1) + A(i, j + 1) - h2 * f(x(i), y(j)));
                }
            }
        }
        // ---- END JACOBI UPDATE ----

        // ---- BEGIN COMPUTE GLOBAL RMSE ----
        double globalResiduals = 0;
        double localResiduals = residuals(smallGridA, smallGridB, W, H);
        MPI_Allreduce(&localResiduals, &globalResiduals, 1, MPI_DOUBLE, MPI_SUM, cartesian);
        rmse = sqrt((1.0 / (L * L)) * globalResiduals);
        // ---- END COMPUTE GLOBAL RMSE ----

        // ...
        std::memcpy(smallGridA, smallGridB, ((W + 2) * (H + 2)) * sizeof(double));
        step++;
    }
    // ---- END JACOBI METHOD ----

    // ---- BEGIN GATHER APPROXIMATE SOLUTIONS FROM ALL PROCESSES ----
    //
    auto localGrid = (double *) malloc(sizeof(double) * (L * L));
    auto globalGrid = (double *) malloc(sizeof(double) * (L * L));
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            localGrid[j * L + i] = 0;
            globalGrid[j * L + i] = 0;
        }
    }

    // ...
    for (int i = 1; i <= W; i++) {
        for (int j = 1; j <= H; j++) {
            localGrid[y(j) * L + x(i)] = A(i, j);
        }
    }

    // ...
    MPI_Reduce(&localGrid[0], &globalGrid[0], L * L, MPI_DOUBLE, MPI_SUM, ROOT_RANK, cartesian);
    // ---- END GATHER APPROXIMATE SOLUTIONS FROM ALL PROCESSES ----

    // Write approximate solution on a csv file and draw 2D visualization on a bitmap file
    if (world_rank == ROOT_RANK) {
        Csv::write(CSV_FILENAME, globalGrid, L, L);
        Bitmap::write(BMP_FILENAME, globalGrid, L, L);
    }

    // Finalize MPI Environment
    MPI_Type_free(&VERTICAL_HALO);
    MPI_Type_free(&HORIZONTAL_HALO);
    MPI_Finalize();
    free(globalGrid);
    free(localGrid);
    free(smallGridB);
    free(smallGridA);

    return 0;
}

int coordToRank(int cartesian_x, int cartesian_y, const int cartesian_dimension[2], MPI_Comm cartesian) {

    if (cartesian_x < 0 || cartesian_x >= cartesian_dimension[0] ||
        cartesian_y < 0 || cartesian_y >= cartesian_dimension[1]) {
        return MPI_PROC_NULL;
    }

    int rank;
    int coords[2] = {cartesian_x, cartesian_y};
    MPI_Cart_rank(cartesian, coords, &rank);
    return rank;

}

double residuals(double *u, double *v, int width, int height) {

    double sum2 = 0;
    for (int i = 1; i <= width; i++) {
        for (int j = 1; j <= height; j++) {
            sum2 += pow(u[j * (width + 2) + i] - v[j * (width + 2) + i], 2);
        }
    }

    return sum2;
}