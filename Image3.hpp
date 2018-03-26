#pragma once

#include "ImageConvert.hpp"
#include "PNGRead.hpp"

#include <vector>
#include <string>
#include <limits>
#include <exception>
#include <memory>

#include <string.h>
#include <stdint.h>
#include <assert.h>

#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>

#include <png.h>

template <class _Type, class _RT=double> // TODO: Kernel::FT
class Image3
{
public:
    typedef _Type Type;
    typedef _RT RT;
    typedef std::shared_ptr<Type> Data;

private:
    const int w, h, d;
    RT dx, dy, dz;
    Data im;

public:
    inline Image3(Data im, int w, int h, int d, RT dx=1, RT dy=1, RT dz=1)
        : im(im), w(w), h(h), d(d), dx(dx), dy(dy), dz(dz) { }

    inline int xdim() const { return w; }
    inline int ydim() const { return h; }
    inline int zdim() const { return d; }
    inline RT& vx() { return dx; }
    inline RT& vy() { return dy; }
    inline RT& vz() { return dz; }
    inline RT vx() const { return dx; }
    inline RT vy() const { return dy; }
    inline RT vz() const { return dz; }
    inline Type* data() { return im; }
    inline const Type* data() const { return im; }

    static Image3* read_png_stack(const std::vector<std::string>& filenames, RT dx=1, RT dy=1, RT dz=1)
    {
        if (filenames.size() > std::numeric_limits<int>::max())
        {
            throw std::invalid_argument("number of slices is too large");
        }
        int w = 0, h = 0, d = (int)filenames.size(), w2, h2, imsz;
        Data data;
        for (size_t i = 0; i < d; ++i)
        {
            Data im = read_png<Type>(filenames[i], w2, h2);
            if (w == 0)
            {
                // First slice
                w = w2; h = h2; imsz = w*h*sizeof(Type);
                data.reset(new Type[w*h*d], std::default_delete<Type[]>());
            }
            else if (w != w2 || h != h2) { throw std::invalid_argument("slices not the same size"); }
            memcpy(data + imsz*i, im, imsz);
        }
        return new Image3(data, w, h, d, dx, dy, dz);
    }
    static Image3* read_mrc_file(std::string& filename, RT dx=0, RT dy=0, RT dz=0)
    {
        static const int32_t IMOD = 0x444F4D49; // "IMOD"
        static const byte MAP_[4] = "MAP_";
        static const byte LE[3][4] = { "\x44\x41\x00\x00", "\x44\x00\x00\x00", "\x44\x44\x00\x00" }; // little endian stamps
        static const byte BE[2][4] = { "\x17\x17\x00\x00", "\x17\x00\x00\x00" }; // big endian stamps
        static const int32_t MODE_UINT8 = 0; // sometimes INT8
        static const int32_t MODE_INT16 = 1, MODE_UINT16 = 6, MODE_FLOAT32 = 2, MODE_RGB24 = 16;
        static const int32_t MODE_INT16_2 = 3, MODE_FLOAT32_2 = 4, MODE_UINT4 = 101; // unsupported modes
        static const int32_t FLAG_SIGNED_BYTES = 0x01;
        static const int32_t FLAG_PX_SPACING_EXT_HDR = 0x02;
        static const int32_t FLAG_ORIGIN_INVERTED = 0x04;
        static const int32_t FLAG_RMS_NEG_NA = 0x08;
        static const int32_t FLAG_BYTE_NYBBLE_SWAP = 0x10;
        static const int32_t FLAG_MASK = FLAG_SIGNED_BYTES | FLAG_PX_SPACING_EXT_HDR | FLAG_ORIGIN_INVERTED | FLAG_RMS_NEG_NA | FLAG_BYTE_NYBBLE_SWAP;
        assert(sizeof(float) == 4);

        // Open the file and memory map it
        const int fd = open(filename.c_str(), O_RDONLY);
        if (fd == -1) { throw std::invalid_argument("cannot open " + filename); }
        Destructor _close([fd] () { close(fd); });
        struct stat sb;
        if (fstat(fd, &sb) == -1) { throw std::invalid_argument("cannot stat " + filename); }
        const size_t len = sb.st_size;
        if (len < 1024) { throw std::invalid_argument(filename + " too small"); }
        void* ptr = mmap(NULL, len, PROT_READ, MAP_PRIVATE, fd, 0);
        Destructor _munmap([ptr, len] () { munmap(ptr, len); });

        // Get the header values
        int32_t* int32s = (int32_t*)ptr;
        float* floats = (float*)ptr;
        int32_t nx = int32s[0], ny = int32s[1], nz = int32s[2], mode = int32s[3];
        int32_t mx = int32s[7], my = int32s[8], mz = int32s[9];
        float xlen = floats[10], ylen = floats[11], zlen = floats[12];
        int32_t mapc = int32s[16], mapr = int32s[17], maps = int32s[18];
        int32_t next = int32s[23];
        int32_t imodStamp = int32s[38], imodFlags = int32s[39];
        bytes cmap = (bytes)&int32s[52], stamp = (bytes)&int32s[53];

        bool big_endian = false;
        if (memcmp(cmap, MAP_, 4) == 0)
        {
            big_endian = memcmp(stamp, BE[0], 4) == 0 || memcmp(stamp, BE[1], 4) == 0;
            if (!big_endian && !(memcmp(stamp, LE[0], 4) == 0 || memcmp(stamp, LE[1], 4) == 0 || memcmp(stamp, LE[2], 4) == 0))
            {
                throw std::invalid_argument(filename + " has invalid MRC header");
            }
        }
        if (big_endian)
        {
            // TODO: swap ordering of all values...
            // Would also need image converter that supported byte-swapping
            throw std::invalid_argument(filename + " is big-endian which is not supported");
        }

        // Check the header values
        if (nx <= 0 || ny <= 0 || nz <= 0 || mx <= 0 || my <= 0 || mz <= 0 ||
            xlen <= 0 || ylen <= 0 || zlen <= 0 || next < 0)
        {
            throw std::invalid_argument(filename + " has invalid MRC header");
        }
        if (mode == MODE_INT16_2 || mode == MODE_FLOAT32_2 || mode == MODE_UINT4)
        {
            throw std::invalid_argument(filename + " has unsupported mode");
        }
        if (mapc != 1 && mapr != 2 && maps != 3)
        {
            throw std::invalid_argument(filename + " has unsupported value ordering");
        }
        int32_t flags = imodStamp == IMOD ? imodFlags : 0;
        if (flags & FLAG_MASK) { throw std::invalid_argument(filename + " has unknown flags"); }

        // Get the data offset
        size_t off = 1024 + next, data_len = len - off, n_px = nx * ny * nz;
        void* data_raw = (bytes)ptr + off;

        // Check mode and get data
        typedef Data (*convert) (void* row, size_t n);
        convert conv = NULL;
        size_t bpp;
        if (mode == MODE_UINT8)
        {
            if (flags & FLAG_BYTE_NYBBLE_SWAP) { throw std::invalid_argument(filename + " has unsupported flag"); }
            bpp = sizeof(uint8_t);
            conv = (flags & FLAG_SIGNED_BYTES) ? convert_im<int8_t, Type> : convert_im<uint8_t, Type>;
        }
        else if (mode == MODE_UINT16)  { bpp = sizeof(uint16_t); conv = convert_im<uint16_t, Type>; }
        else if (mode == MODE_INT16)   { bpp = sizeof(int16_t);  conv = convert_im<int16_t, Type>; }
        else if (mode == MODE_FLOAT32) { bpp = sizeof(float);    conv = convert_im<float, Type>; }
        else if (mode == MODE_RGB24)   { bpp = sizeof(RGB24);    conv = convert_im<RGB24, Type>; }
        else { throw std::invalid_argument(filename + " has unknown mode"); }
        if (data_len != bpp*n_px) { throw std::invalid_argument(filename + " has wrong data size"); }

        // Get image
        Data data = conv(data_raw, n_px);
        if (!dx) { dx = RT(xlen) / mx; }
        if (!dy) { dy = RT(ylen) / my; }
        if (!dz) { dz = RT(zlen) / mz; }
        return new Image3(data, nx, ny, nz, dx, dy, dz);
    }
};
