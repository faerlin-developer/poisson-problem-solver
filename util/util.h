#pragma once

/**
 *
 * @param cartesian_x
 * @param cartesian_y
 * @param cartesian_dimension
 * @param cartesian
 * @return
 */
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