#pragma once

#include <string>
#include <cstdint>
#include <memory>

class Bitmap {
public:
    static bool write(const std::string &filename, double *grid, int width, int height);

private:
    static void setPixel(int x, int y, uint8_t red, uint8_t green, uint8_t blue, uint8_t *pixels, int width);
};

// Align members of BitmapFileHeader struct to 2-byte boundaries or natural alignment boundary, whichever is less
// For BitmapFileHeader struct in particular, there will be no padding.
#pragma pack(push, 2)
struct BitmapFileHeader {
    char header[2]{'B', 'M'};
    int32_t fileSize;
    int32_t reserved{0};
    int32_t dataOffset;
};
#pragma pack(pop)

/*
 * (https://www.geeksforgeeks.org/structure-member-alignment-padding-and-data-packing/)
 * - #pragma pack(n) means align members of struct to n-byte boundaries or natural alignment boundary, whichever is
 *  less.
 * - Aligning to an n-byte boundary means having an address that is a multiple of n.
 * - For example, a member aligned on a 1-byte boundary means it can have any address. A member aligned on a 2-byte
 * - boundary means it can only have address that is a multiple of 2.
 * - Additional padding may be inserted to ensure that an array of the struct type has the intended alignment. In such
 *   a case, the struct is padded such that it meets m-byte boundary where m is the largest size boundary among its
 *   members.
 */

// Align members of BitmapInfoHeader struct to 2-byte boundaries or natural alignment boundary, whichever is less
// For BitmapInfoHeader struct in particular, there will be no padding.
#pragma pack(push, 2)
struct BitmapInfoHeader {
    int32_t headerSize{40};
    int32_t width;
    int32_t height;
    int16_t planes{1};
    int16_t bitsPerPixel{24};
    int32_t compression{0};
    int32_t dataSize{0};
    int32_t horizontalResolution{2400};
    int32_t verticalResolution{2400};
    int32_t colors{0};
    int32_t importantColors{0};
};
#pragma pack(pop)


