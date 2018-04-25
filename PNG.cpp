#include "Image3.hpp"

#include <string.h>
#include <stdio.h>

#include <CGAL/ImageIO.h>

#include <png.h>

static int testPngHeader(char *magic, const char *name)
{
    return !png_sig_cmp((png_const_bytep)magic, 0, 4) ? ImageIO_NO_ERROR : ImageIO_UNKNOWN_TYPE;
}

static void _readPngImage_off_scale(png_structp png_ptr, png_infop info_ptr, _image *im)
{
    int unit; // PNG_SCALE_METER or PNG_SCALE_RADIAN if has_scale or
              // PNG_RESOLUTION_METER or PNG_RESOLUTION_UNKNOWN if has_phys
    double width, height; // sCAL
    png_uint_32 res_x, res_y; // pHYs
    bool has_scale, has_phys = false;
    if ((has_scale = png_get_sCAL(png_ptr, info_ptr, &unit, &width, &height)) ||
        (has_phys = png_get_pHYs(png_ptr, info_ptr, &res_x, &res_y, &unit)))
    {
        im->vx = has_scale ? width : res_x;
        im->vy = has_scale ? height : res_y;
        im->vz = 1;
    }
    else { im->vx = im->vy = im->vz = 1; }
    
    int off_unit; // PNG_OFFSET_PIXEL or PNG_OFFSET_MICROMETER
    png_int_32 x_off, y_off;
    bool has_off = png_get_oFFs(png_ptr, info_ptr, &x_off, &y_off, &off_unit);
    if (has_off)
    {
        im->tx = x_off;
        im->ty = y_off;
        im->tz = 0;
        // Try to make sure the units are consistent with the scale/physical units
        if (off_unit == PNG_OFFSET_PIXEL) { im->tx *= im->vx; im->ty *= im->vy; }
        else if (off_unit == PNG_OFFSET_MICROMETER || ((has_scale && unit == PNG_SCALE_METER) || (has_phys && unit == PNG_RESOLUTION_METER)))
        {
            im->tx *= 1000000; im->ty *= 1000000; // um to m
        }
    }
    else { im->tx = im->ty = im->tz = 0; }
}

static bool _readPngImage_labels(png_structp png_ptr, png_infop info_ptr, _image *im)
{
    png_textp text_ptr;
    im->nuser = png_get_text(png_ptr, info_ptr, &text_ptr, NULL);
    if (im->nuser)
    {
        im->user = (char**)ImageIO_alloc(im->nuser * sizeof(char*));
        if (!im->user) { return false; }
        memset(im->user, 0, im->nuser * sizeof(char*));
    }
    for (unsigned int i = 0; i < im->nuser; ++i)
    {
        int len = (strlen(text_ptr[i].key) + strlen(text_ptr[i].text) + 3) * sizeof(char);
        im->user[i] = (char*)ImageIO_alloc(len);
        if (!im->user[i]) { return false; }
        snprintf(im->user[i], len, "%s: %s", text_ptr[i].key, text_ptr[i].text);
    }
    return true;
}

static int _readPngImage(const char *name, _image *im, FILE* f)
{
    // Setup PNG structures
    png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png_ptr) { fprintf(stderr, "readPngImage: failed to initialize libpng structures\n"); return ImageIO_READING_HEADER; }
    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) { png_destroy_read_struct(&png_ptr, NULL, NULL); fprintf(stderr, "readPngImage: failed to initialize libpng structures\n"); return ImageIO_READING_HEADER; }
    if (setjmp(png_jmpbuf(png_ptr)))
    {
        png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
        fprintf(stderr, "readPngImage: failed to read PNG image\n");
        return im->data ? ImageIO_READING_IMAGE : ImageIO_READING_HEADER;
    }
    png_init_io(png_ptr, f);
    png_set_sig_bytes(png_ptr, 7);

    // Read the PNG data into the internal data structure
    // PACKING causes 1/2/4-bit color samples to become 8-bit samples
    // EXPAND causes palette to become RGB and grayscale image 1/2/4-bit images to 8-bit, and resolves tRNS chunks
    // SWAP_ENDIAN causes 16-bit samples to be low endian
    png_read_png(png_ptr, info_ptr, PNG_TRANSFORM_PACKING | PNG_TRANSFORM_EXPAND |
        (::_getEndianness() == END_LITTLE ? PNG_TRANSFORM_SWAP_ENDIAN : 0), NULL);

    // Get information about the image
    png_uint_32 h, w; //, stride = png_get_rowbytes(png_ptr, info_ptr);
    int bit_depth, channels = png_get_channels(png_ptr, info_ptr);
    png_get_IHDR(png_ptr, info_ptr, &w, &h, &bit_depth, NULL, NULL, NULL, NULL);
   
    // Fill in image information
    im->endianness = ::_getEndianness();
    im->vectMode = channels == 1 ? VM_SCALAR : VM_INTERLACED;
    im->wordKind = WK_FIXED; im->sign = SGN_UNSIGNED;
    im->wdim = bit_depth / 8; // 1 or 2
    im->xdim = w; im->ydim = h; im->zdim = 1; im->vdim = channels;
    im->rx = im->ry = im->rz = 0;
    im->cx = im->cy = im->cz = 0;

    // Possibly get scale and offset
    _readPngImage_off_scale(png_ptr, info_ptr, im);
    
    // Fill in labels
    if (!_readPngImage_labels(png_ptr, info_ptr, im))
    {
        png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
        fprintf(stderr, "readPngImage: failed to allocate memory\n");
        return ImageIO_READING_HEADER;
    }

    // Get the image data
    size_t nbytes = ((size_t)im->xdim) * im->ydim * im->zdim * im->vdim * im->wdim;
    im->data = ImageIO_alloc(nbytes);
    if (!im->data)
    {
        png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
        fprintf(stderr, "readPngImage: failed to allocate image memory\n");
        return ImageIO_READING_IMAGE;
    }
    png_bytepp rows = png_get_rows(png_ptr, info_ptr);
    size_t row_sz = im->xdim * im->vdim * im->wdim;
    for (size_t i = 0; i < im->ydim; ++i)
    {
        memcpy(((unsigned char*)im->data) + i*row_sz, rows[i], row_sz);
    }

    return 1;
}

static int readPngImage(const char *name, _image *im)
{
    // Open the file
    FILE *f = fopen(name, "rb");
    if (!f)
    {
        fprintf(stderr, "readPngImage: failed to open %s\n", name);
        return ImageIO_OPENING;
    }

    // Check the signature
    unsigned char buf[7];
    if (fread(buf, 1, 7, f) != 7 || png_sig_cmp(buf, 0, 7))
    {
        fprintf(stderr, "readPngImage: %s is not a PNG\n", name);
        fclose(f);
        return ImageIO_UNKNOWN_TYPE;
    }

    // Read the data
    im->data = im->user = NULL;
    int result = _readPngImage(name, im, f);

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

void registerPNGFormat()
{
    IMAGE_FORMAT* format = (IMAGE_FORMAT*)ImageIO_alloc(sizeof(IMAGE_FORMAT));
    format->testImageFormat = &testPngHeader;
    format->readImageHeader = &readPngImage;
    format->writeImage = NULL; // TODO
    strncpy(format->fileExtension, ".png", IMAGE_FORMAT_NAME_LENGTH);
    strncpy(format->realName, "PNG", IMAGE_FORMAT_NAME_LENGTH);
    addImageFormat(format);
}
