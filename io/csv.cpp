#include <fstream>
#include <string>
#include "csv.h"

/**
 *
 * @param grid
 */
void Csv::write(const std::string &filename, double *grid, int W, int H) {

    std::ofstream file;
    file.open(filename, std::ios::trunc);
    for (int j = H - 1; j >= 0; j--) {
        file << std::to_string(grid[j * W + 0]);
        for (int i = 1; i < W; i++) {
            file << "," << std::to_string(grid[j * W + i]);
        }
        file << std::endl;
    }
    file.close();
}
