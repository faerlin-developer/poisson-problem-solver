#pragma once

#include <cstdint>

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
