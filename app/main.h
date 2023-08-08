#pragma once

/**
 *
 * @param u
 * @param v
 * @param width
 * @param height
 * @return
 */
double residuals(double *u, double *v, int width, int height);

/**
 *
 * @param cartesian_x
 * @param cartesian_y
 * @param cartesian_dimension
 * @param cartesian
 * @return
 */
int coordToRank(int cartesian_x, int cartesian_y, const int cartesian_dimension[2], MPI_Comm cartesian);