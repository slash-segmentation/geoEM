#include "Image3.hpp"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits>

_image* copy_image_header(const _image* src)
{
    // Allocate the new image and copy attributes
    _image* im = ::_initImage();
    if (!im) { return NULL; }
    im->wdim = src->wdim; im->vectMode = src->vectMode; im->wordKind = src->wordKind;
    im->sign = src->sign; im->endianness = src->endianness;
    im->xdim = src->xdim; im->ydim = src->ydim; im->zdim = src->zdim; im->vdim = src->vdim;
    im->vx = src->vx; im->vy = src->vy; im->vz = src->vz;
    im->tx = src->tx; im->ty = src->ty; im->tz = src->tz;
    im->rx = src->rx; im->ry = src->ry; im->rz = src->rz;
    im->cx = src->cx; im->cy = src->cy; im->cz = src->cz;
    
    // Copy user strings
    im->nuser = src->nuser;
    im->user = (char**)ImageIO_alloc(im->nuser*sizeof(char*));
    if (!im->user) { ::_freeImage(im); throw std::bad_alloc(); }
    memset(im->user, 0, im->nuser*sizeof(char*));
    for (size_t i = 0; i < im->nuser; ++i)
    {
        im->user[i] = (char*)ImageIO_alloc(strlen(src->user[i])*sizeof(char));
        if (!im->user[i]) { ::_freeImage(im); return NULL; }
        strcpy(im->user[i], src->user[i]);
    }
    
    return im;
}

void copy_image(const Image3& image, Image3& output)
{
    const _image* src = image.image();
    
    // Copy the old image's header
    _image* im = copy_image_header(src);
    if (!im) { throw std::bad_alloc(); }
    
    // Copy image data
    size_t nbytes = ((size_t)im->xdim) * im->ydim * im->zdim * im->vdim * im->wdim;
    im->data = ImageIO_alloc(nbytes);
    if (!im->data) { ::_freeImage(im); throw std::bad_alloc(); }
    memcpy(im->data, src->data, nbytes);

    // Done!
    output = Image3(im);
}

bool fp_eq(double a, double b)
{
    double max = std::max(fabs(a), fabs(a));
    return fabs(a-b) <= (max * 4 * std::numeric_limits<double>::epsilon());
}

bool stack_images(const std::vector<Image3>& images, Image3& output)
{
    if (images.size() == 0) { return false; }
    const _image* first = images[0].image();
    
    // Check the data
    size_t zdim = 0, i = 0;
    bool consistent = true;
    for (const Image3& slc : images)
    {
        const _image* im = slc.image();
        if (im->wdim != first->wdim || im->vectMode != first->vectMode || im->wordKind != first->wordKind ||
            im->sign != first->sign || im->endianness != first->endianness ||
            im->xdim != first->xdim || im->ydim != first->ydim || im->vdim != first->vdim)
        {
            fprintf(stderr, "cannot stack images of different types or dimensions\n");
            return false;
        }
        if (consistent &&
            (im->vx != first->vx || im->vy != first->vy || im->vz != first->vz ||
             im->tx != first->tx || im->ty != first->ty || !fp_eq(im->tz, first->tz + i*first->vz)))
        {
            fprintf(stderr, "Warning: stacking images of different voxel sizes or different translations\n");
            consistent = false;
        }
        zdim += im->zdim;
    }
    
    // Allocate the new image and copy attributes
    _image* im = copy_image_header(first);
    if (!im) { throw std::bad_alloc(); }
    im->zdim = zdim;
    if (!consistent) { im->vx = im->vy = im->vz = 1; im->tx = im->ty = im->tz = 0; }
    
    // TODO: these are ignored
    im->rx = im->ry = im->rz = 0;
    im->cx = im->cy = im->cz = 0;
    
    // TODO: deal with nuser and user[...] (first is copied but the rest are ignored)
    
    // Copy image data
    size_t bpp = im->vdim * im->wdim;
    size_t nbytes = ((size_t)im->xdim) * im->ydim * im->zdim * bpp;
    unsigned char* data = (unsigned char*)(im->data = ImageIO_alloc(nbytes));
    if (!data) { ::_freeImage(im); throw std::bad_alloc(); }
    for (const Image3& slc : images)
    {
        size_t slc_nbytes = slc.size() * bpp;
        memcpy(data, slc.data(), slc_nbytes);
        data += slc_nbytes;
    }

    // Done!
    output = Image3(im);
    return true;
}

void fix_endian(const Image3& image, Image3& output)
{
    copy_image(image, output);
    fix_endian(output);
}

#if defined(_MSC_VER)
	uint16_t inline byte_swap(uint16_t x) { return _byteswap_ushort(x); }
	uint32_t inline byte_swap(uint32_t x) { return _byteswap_ulong(x);  }
	uint64_t inline byte_swap(uint64_t x) { return _byteswap_uint64(x); }
#elif defined(__GNUC__) // GCC and Clang
	uint16_t inline byte_swap(uint16_t x) { return __builtin_bswap16(x); }
	uint32_t inline byte_swap(uint32_t x) { return __builtin_bswap32(x); }
	uint64_t inline byte_swap(uint64_t x) { return __builtin_bswap64(x); }
#else
	uint16_t inline byte_swap(uint16_t x) { return (x<<8)|(x>>8); }
	uint32_t inline byte_swap(uint32_t x) { return (x<<24)|((x<<8)&0x00FF0000)|((x>>8)&0x0000FF00)|(x>>24); }
	uint64_t inline byte_swap(uint64_t x) { return (x<<56)|((x<<40)&0x00FF000000000000ull)|((x<<24)&0x0000FF0000000000ull)|((x<<8)&0x000000FF00000000ull)|((x>>8)&0x00000000FF000000ull)|((x>>24)&0x0000000000FF0000)|((x>>40)&0x000000000000FF00)|(x>>56); }
#endif

void fix_endian(Image3& image)
{
    _image* im = image.image();
    if (im->endianness == ::_getEndianness()) { return; } // nothing to fix
    im->endianness = ::_getEndianness();
    size_t wdim = im->wdim, nvals = ((size_t)im->xdim) * im->ydim * im->zdim * im->vdim;
    if (wdim <= 1) { return; } // nothing to fix
    // Highly optimized endian conversions
    else if (wdim == 2)
    {
        uint16_t* data = (uint16_t*)im->data;
        for (size_t i = 0; i < nvals; ++i) { data[i] = byte_swap(data[i]); }
    }
    else if (wdim == 4)
    {
        uint32_t* data = (uint32_t*)im->data;
        for (size_t i = 0; i < nvals; ++i) { data[i] = byte_swap(data[i]); }
    }
    else if (wdim == 8)
    {
        uint64_t* data = (uint64_t*)im->data;
        for (size_t i = 0; i < nvals; ++i) { data[i] = byte_swap(data[i]); }
    }
    else
    {
        // Generic endian conversions
        unsigned char* data = (unsigned char*)im->data;
        size_t nbytes = wdim * nvals;
        for (size_t i = 0; i < nbytes; i += wdim)
        {
            for (size_t j = 0; j < wdim / 2; ++j)
            {
                unsigned char temp = data[i+j];
                data[i+j] = data[i+wdim-j-1];
                data[i+wdim-j-1] = temp;
            }
        }
    }
}
