#include <fstream>
#include "bitmap.h"
#include "bitmapinfoheader.h"
#include "bitmapfileheader.h"

Bitmap::Bitmap(int width, int height) : m_width(width), m_height(height), m_pPixels(new uint8_t[width * height * 3]{}) {
    // the () or {} at the end of new uint8_t[width * height * 3] initializes each cell to zero
}

bool Bitmap::write(std::string filename) {
    BitmapFileHeader fileHeader;
    BitmapInfoHeader infoHeader;

    fileHeader.fileSize = sizeof(BitmapFileHeader) + sizeof(BitmapInfoHeader) + m_width * m_height * 3;
    fileHeader.dataOffset = sizeof(BitmapFileHeader) + sizeof(BitmapInfoHeader);

    infoHeader.width = m_width;
    infoHeader.height = m_height;

    std::ofstream file;
    file.open(filename.c_str(), std::ios::out | std::ios::binary);
    if (!file) {
        return false;
    }

    file.write((char *) &fileHeader, sizeof(fileHeader));
    file.write((char *) &infoHeader, sizeof(infoHeader));
    file.write((char *) m_pPixels.get(), m_width * m_height * 3);

    file.close();
    if (!file) {
        return false;
    }

    return true;
}

// Bitmap coordinates is cartesian coordinates
void Bitmap::setPixel(int x, int y, uint8_t red, uint8_t green, uint8_t blue) {

    uint8_t *pPixel = m_pPixels.get();
    pPixel += (y * 3) * m_width + x * 3;

    // Bitmap expects little endian format
    pPixel[2] = red;
    pPixel[1] = green;
    pPixel[0] = blue;
}
