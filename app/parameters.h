#pragma once

//
int L = 60;

//
int W = 20;
int H = 15;

//
std::string OUTPUT_BMP = "solution.bmp";

//
double TARGET_RMS = 1e-7;

// Maximum number of iteration for the Jacobi method
int MAX_ITERATION = 50 * L;

/**
 *
 * @param x
 * @param y
 * @param L
 * @return
 */
double f(int x, int y) {

    if (x < 0 || x >= L || y < 0 || y >= L) {
        throw std::logic_error("grid coordinates are out of bounds");
    }

    return 0;
}

/**
 *
 * @param x
 * @param y
 * @param L
 * @return
 */
double g(int x, int y) {

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
