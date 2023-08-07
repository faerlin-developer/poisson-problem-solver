#pragma once

#include <cstdint>

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

