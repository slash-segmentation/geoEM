#include "Image3.hpp"

#include <string.h>
#include <stdio.h>
#include <stdint.h>

#include <CGAL/ImageIO.h>

typedef unsigned char byte;

static const int32_t MRC_IMOD = 0x444F4D49; // "IMOD"
static const byte MRC_MAP_[4] = {'M', 'A', 'P', ' '};
static const byte MRC_LE[3][4] = {{'\x44','\x41','\x00','\x00'}, {'\x44','\x00','\x00','\x00'}, {'\x44','\x44','\x00','\x00'}}; // little endian stamps
static const byte MRC_BE[2][4] = {{'\x17','\x17','\x00','\x00'}, {'\x17','\x00','\x00','\x00'}}; // big endian stamps
static const int32_t MRC_MODE_BYTE = 0, MRC_MODE_INT16 = 1, MRC_MODE_FLOAT32 = 2;
static const int32_t MRC_MODE_INT16_2 = 3, MRC_MODE_FLOAT32_2 = 4, MRC_MODE_UINT16 = 6, MRC_MODE_RGB24 = 16;
static const int32_t MRC_MODE_UINT4 = 101; // unsupported mode
static const int32_t MRC_FLAG_SIGNED_BYTES = 0x01;
static const int32_t MRC_FLAG_PX_SPACING_EXT_HDR = 0x02;
static const int32_t MRC_FLAG_ORIGIN_INVERTED = 0x04;
static const int32_t MRC_FLAG_RMS_NEG_NA = 0x08;
static const int32_t MRC_FLAG_BYTE_NYBBLE_SWAP = 0x10;
static const int32_t MRC_FLAG_MASK = MRC_FLAG_SIGNED_BYTES | MRC_FLAG_PX_SPACING_EXT_HDR | MRC_FLAG_ORIGIN_INVERTED | MRC_FLAG_RMS_NEG_NA | MRC_FLAG_BYTE_NYBBLE_SWAP;

static ENDIANNESS getMrcEndianness(byte* header, bool& new_header_style)
{
    byte *cmap = &header[208], *stamp = &header[212];
    if (memcmp(cmap, MRC_MAP_, 4) == 0)
    {
        new_header_style = true;
        if (memcmp(stamp, MRC_BE[0], 4) == 0 || memcmp(stamp, MRC_BE[1], 4) == 0) { return END_BIG; }
        if (memcmp(stamp, MRC_LE[0], 4) == 0 || memcmp(stamp, MRC_LE[1], 4) == 0 || memcmp(stamp, MRC_LE[2], 4) == 0) { return END_LITTLE; }
        return END_UNKNOWN;
    }
    else { new_header_style = false; }
    return END_LITTLE;
}
static void fixMrcEndianness(byte* header)
{
    // This universal byte swapping will mess up some headers but we don't care about those
    // Mostly related to the extended or old-style headers
    int32_t* xs = (int32_t*)header;
    for (int i = 0; i < 56; ++i)
    {
        int32_t x = xs[i]; // with GCC the below code will become a single CPU instruction if available
        xs[i] = ((x & 0xff000000) >> 24) | ((x & 0xff0000) >> 8) | ((x & 0xff00) << 8) | (x << 24);
    }
}
static bool readMrcHeader(FILE *f, _image* im, bool include_labels)
{
    byte* header = (byte*)ImageIO_alloc(1024);
    if (!header) { return false; }
    if (fread(header, 1, 1024, f) != 1024) { ImageIO_free(header); return false; }

    // Check if data is big or little endian
    bool new_header_style;
    ENDIANNESS endianness = getMrcEndianness(header, new_header_style);
    if (endianness == END_UNKNOWN) { ImageIO_free(header); return false; }
    if (::_getEndianness() != endianness) { fixMrcEndianness(header); }

    // Get the header values
    int32_t* int32s = (int32_t*)header;
    float* floats = (float*)header;
    int32_t nx = int32s[0], ny = int32s[1], nz = int32s[2]; // image size
    int32_t mode = int32s[3]; // pixel type
    // ignoring: int32_t nxstart = int32s[4], nystart = int32s[5], mnzstartz = int32s[6]; // starting point of sub image
    int32_t mx = int32s[7], my = int32s[8], mz = int32s[9]; // grid size
    float xlen = floats[10], ylen = floats[11], zlen = floats[12]; // cell size
    // ignoring: float alpha = floats[13], beta = floats[14], gamma = floats[15]; // cell angles
    int32_t mapc = int32s[16], mapr = int32s[17], maps = int32s[18]; // data ordering
    // ignoring: float amin = floats[19], amax = floats[20], amean = floats[21];
    // ignoring: int32_t ispg = int32s[22]; // 0 for image stack, 1 for volume
    int32_t next = int32s[23]; // extended header size in bytes
    // ignoring the next fields as they aren't used or pertain to interpretation of extended header
    // TODO: should be checking version and only use everything after this point if imodStamp is IMOD
    int32_t imodStamp = int32s[38], imodFlags = int32s[39]; // file creator and flags
    // ignoring the next fields as they are acquisition information
    float xorg, yorg, zorg; // the origin of the image
    if (new_header_style) { xorg = floats[49]; yorg = floats[50]; zorg = floats[51]; }
    else { xorg = floats[52]; yorg = floats[53]; zorg = floats[54]; }
    // The only other new-style header field is rms which we will ignore
    int32_t nlabl = int32s[55];

    // Check the header values
    if (nx <= 0 || ny <= 0 || nz <= 0 || mx <= 0 || my <= 0 || mz <= 0 ||
        xlen <= 0 || ylen <= 0 || zlen <= 0 || next < 0 || nlabl < 0 || nlabl > 10 ||
        mapc != 1 || mapr != 2 || maps != 3) { ImageIO_free(header); return false; }
    int32_t flags = imodStamp == MRC_IMOD ? imodFlags : 0;
    if (flags & MRC_FLAG_MASK) { ImageIO_free(header); return false; }

    // Fill in basic image information
    im->endianness = endianness;
    im->xdim = nx; im->ydim = ny; im->zdim = nz; im->vdim = 1;
    im->vx = (double)xlen / mx; im->vy = (double)ylen / mx; im->vz = (double)zlen / mx;
    im->tx = -xorg, im->ty = -yorg; im->tz = -zorg; // TODO: negative or positve? use MRC_FLAG_ORIGIN_INVERTED
    im->rx = im->ry = im->rz = 0; // TODO: image rotation
    im->cx = im->cy = im->cz = 0; // TODO: image center

    // Fill in the pixel mode/kind/sign/size
    // Default mode/kind/size
    im->vectMode = VM_SCALAR;
    im->wordKind = WK_FIXED;
    im->sign = SGN_SIGNED;
    if (mode == MRC_MODE_BYTE)
    {
        if (flags & MRC_FLAG_BYTE_NYBBLE_SWAP) { ImageIO_free(header); return false; }
        im->wdim = sizeof(uint8_t);
        if (!(flags & MRC_FLAG_SIGNED_BYTES)) { im->sign = SGN_UNSIGNED; }
    }
    else if (mode == MRC_MODE_INT16)     { im->wdim = sizeof(int16_t); }
    else if (mode == MRC_MODE_FLOAT32)   { im->wdim = sizeof(float); im->wordKind = WK_FLOAT; }
    else if (mode == MRC_MODE_INT16_2)   { im->wdim = sizeof(int16_t); im->vectMode = VM_INTERLACED; im->vdim = 2; }
    else if (mode == MRC_MODE_FLOAT32_2) { im->wdim = sizeof(float); im->vectMode = VM_INTERLACED; im->vdim = 2; im->wordKind = WK_FLOAT; }
    else if (mode == MRC_MODE_UINT16)    { im->wdim = sizeof(uint16_t); im->sign = SGN_UNSIGNED; }
    else if (mode == MRC_MODE_RGB24)     { im->wdim = sizeof(uint8_t); im->vectMode = VM_INTERLACED; im->sign = SGN_UNSIGNED; im->vdim = 3; }
    else { ImageIO_free(header); return false; }

    // Fill in the labels
    if (include_labels)
    {
        im->nuser = nlabl;
        if (nlabl)
        {
            im->user = (char**)ImageIO_alloc(nlabl * sizeof(char*));
            if (!im->user) { ImageIO_free(header); return false; }
            memset(im->user, 0, nlabl * sizeof(char*));
        }
    	for (int i = 0; i < nlabl; ++i)
        {
            im->user[i] = (char*)ImageIO_alloc(80 * sizeof(char));
            if (!im->user[i]) { ImageIO_free(header); return false; }
            memcpy(im->user[i], header+224+i*80*sizeof(char), 80*sizeof(char));
    	}
    }

    // Everything is in im now
    ImageIO_free(header);

    // Skip extended header
    if (next && !fseek(f, next, SEEK_CUR)) { return false; }

    // Success!
    return !feof(f);
}

static int testMrcHeader(char *magic, const char *name)
{
    // This simply uses readMrcHeader (without labels) to check if it is a valid MRC file that we
    // can read. There is no way to tell if the file is valid from the first 4 bytes.
    FILE* f = fopen(name, "rb");
    if (!f) { return ImageIO_OPENING; }
    _image im;
    int result = readMrcHeader(f, &im, false) ? ImageIO_NO_ERROR : ImageIO_UNKNOWN_TYPE;
    fclose(f);
    return result;
}

static int _readMrcImage(const char *name, _image *im, FILE* f)
{
    // Read the header with labels
    if (!readMrcHeader(f, im, true))
    {
        fprintf(stderr, "readMrcImage: %s has invalid or unsupported header\n", name);
        return ImageIO_READING_HEADER;
    }

    // Get the image data
    size_t nbytes = ((size_t)im->xdim) * im->ydim * im->zdim * im->vdim * im->wdim;
    im->data = ImageIO_alloc(nbytes);
    if (!im->data)
    {
        fprintf(stderr, "readMrcImage: failed to allocate image memory\n");
        return ImageIO_READING_IMAGE;
    }
    if (fread(im->data, 1, nbytes, f) != nbytes)
    {
        fprintf(stderr, "readMrcImage: failed to read image data in %s\n", name);
        return ImageIO_READING_IMAGE;
    }

    // TODO: check eof

    return 1;
}

static int readMrcImage(const char *name, _image *im)
{
    // Open the file
    FILE* f = fopen(name, "rb");
    if (!f)
    {
        fprintf(stderr, "readMrcImage: failed to open %s\n", name);
        return ImageIO_OPENING;
    }

    // Read the data
    im->data = im->user = NULL;
    int result = _readMrcImage(name, im, f);

    // Cleanup
    fclose(f);
    if (result <= 0 && im->data) { ImageIO_free(im->data); im->data = NULL; }
    if (result < 0 && im->user)
    {
        for (unsigned int i = 0; i < im->nuser && im->user[i]; ++i) { ImageIO_free(im->user[i]); }
        ImageIO_free(im->user);
        im->user = NULL;
    }

    return result;
}

void registerMRCFormat()
{
    IMAGE_FORMAT* format = (IMAGE_FORMAT*)ImageIO_alloc(sizeof(IMAGE_FORMAT));
    format->testImageFormat = &testMrcHeader;
    format->readImageHeader = &readMrcImage;
    format->writeImage = NULL; // TODO
    strncpy(format->fileExtension, ".mrc,.ali,.rec,.st", IMAGE_FORMAT_NAME_LENGTH);
    strncpy(format->realName, "MRC", IMAGE_FORMAT_NAME_LENGTH);
    addImageFormat(format);
}
