#include <fstream>
#include <cfloat>
#include <cmath>
#include "bitmap.h"
#include "bitmapinfoheader.h"
#include "bitmapfileheader.h"

// Bitmap coordinates is cartesian coordinates
void Bitmap::setPixel(int x, int y, uint8_t red, uint8_t green, uint8_t blue, uint8_t *pixels, int width) {

    pixels += (y * 3) * width + x * 3;

    // Bitmap expects little endian format
    pixels[2] = red;
    pixels[1] = green;
    pixels[0] = blue;
}

bool Bitmap::write(const std::string &filename, double *grid, int width, int height) {

    auto maximum = -DBL_MAX;
    auto minimum = DBL_MAX;
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            maximum = std::max(maximum, grid[j * width + i]);
            minimum = std::min(minimum, grid[j * width + i]);
        }
    }

    auto pixels = std::unique_ptr<uint8_t>(new uint8_t[width * height * 3]{});
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            auto value = std::round(std::abs((grid[j * width + i] - minimum) / (maximum - minimum)) * 255.0);
            Bitmap::setPixel(i, j, value, 0, 0, pixels.get(), width);
        }
    }

    BitmapFileHeader fileHeader;
    BitmapInfoHeader infoHeader;

    fileHeader.fileSize = sizeof(BitmapFileHeader) + sizeof(BitmapInfoHeader) + width * height * 3;
    fileHeader.dataOffset = sizeof(BitmapFileHeader) + sizeof(BitmapInfoHeader);

    infoHeader.width = width;
    infoHeader.height = height;

    std::ofstream file;
    file.open(filename, std::ios::out | std::ios::binary);
    if (!file) {
        return false;
    }

    file.write((char *) &fileHeader, sizeof(fileHeader));
    file.write((char *) &infoHeader, sizeof(infoHeader));
    file.write((char *) pixels.get(), width * height * 3);

    file.close();
    if (!file) {
        return false;
    }

    return true;
}
