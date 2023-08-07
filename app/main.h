#pragma once

/**
 *
 * @param u
 * @param v
 * @param W
 * @param H
 * @return
 */
double squares(double *u, double *v, int W, int H);

/**
 *
 * @param cartesian_x
 * @param cartesian_y
 * @param cartesian_dimension
 * @param cartesian
 * @return
 */
int coordToRank(int cartesian_x, int cartesian_y, const int cartesian_dimension[2], MPI_Comm cartesian);