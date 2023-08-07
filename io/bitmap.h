#pragma once

#include <string>
#include <cstdint>
#include <memory>

class Bitmap {
private:
    std::string filename;

public:
    explicit Bitmap(std::string filename);

    bool write(double *grid, int width, int height);

private:
    static void setPixel(int x, int y, uint8_t red, uint8_t green, uint8_t blue, uint8_t *pixels, int width);
};
