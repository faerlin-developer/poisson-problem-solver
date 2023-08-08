#pragma once

#include <stdexcept>
#include <string>
#include <cmath>

//
int L = 60;

//
int W = 20;
int H = 15;

//
std::string BMP_FILENAME = "solution.bmp";
std::string CSV_FILENAME = "solution.csv";

//
double TARGET_RMSE = 1e-6;

// Maximum number of iteration for the Jacobi method
int MAX_ITERATION = 50 * L;

//
int ROOT_RANK = 0;

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
    if (y == 0 || y == L - 1) {
        return sin(delta * x);
    } else if (x == 0 || x == L - 1) {
        return sin(delta * y);
    } else {
        throw std::logic_error("grid coordinates is not on a boundary");
    }

}
