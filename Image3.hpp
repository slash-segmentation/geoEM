#pragma once

#include <CGAL/ImageIO.h>
#include <CGAL/Image_3.h>

#include "ImageConvert.hpp"

typedef CGAL::Image_3 Image3;

void registerPNGFormat();
void registerMRCFormat();

_image* copy_image_header(const _image* src);

void copy_image(const Image3& image, Image3& output);
bool stack_images(const std::vector<Image3>& images, Image3& output);
void fix_endian(const Image3& image, Image3& output);
void fix_endian(Image3& image);

template <class WordType>
inline void convert_image(const Image3& image, Image3& output)
{
    const _image* src = image.image();
    
    // Fix the endianness before converting
    if (src->endianness != ::_getEndianness())
    {
        Image3 fixed;
        fix_endian(image, fixed);
        convert_image<WordType>(fixed, output);
        return;
    }

    // Copy the old image's header
    _image* im = copy_image_header(src);
    if (!im) { throw std::bad_alloc(); }
    
    // Set values based on type
    set_image_word_type<WordType>(im);
    
    // Convert image data
    size_t nbytes = ((size_t)im->xdim) * im->ydim * im->zdim * im->vdim * im->wdim;
    im->data = ImageIO_alloc(nbytes);
    if (!im->data) { ::_freeImage(im); throw std::bad_alloc(); }
    convert_image_data(src, (WordType*)im->data);

    // Done!
    output = Image3(im);
}
