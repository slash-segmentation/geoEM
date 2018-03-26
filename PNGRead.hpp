#pragma once

#include "Destructor.hpp"
#include "ImageConvert.hpp"

#include <png.h>

#include <string>
#include <memory>
#include <limits>
#include <exception>

bool check_if_png(FILE *fp);
void png_cpp_error_fn(png_structp png_ptr, png_const_charp error_msg);

template <class Type>
std::shared_ptr<Type> read_png(const std::string& filename, int& width, int& height)
{
    // Open the file
    FILE *fp;
    if (!(fp = fopen(filename.c_str(), "rb"))) { throw std::invalid_argument(filename + " cannot be opened"); }
    Destructor _close_fp([fp] () { fclose(fp); });
    if (!check_if_png(fp)) { throw std::invalid_argument(filename + " is nto a PNG"); }

    // Setup PNG structures
    png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, &png_cpp_error_fn, NULL, NULL);
    if (!png_ptr) { throw std::invalid_argument("Failed to create PNG read struct"); }
    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) { png_destroy_read_struct(&png_ptr, NULL, NULL); return throw std::invalid_argument("failed to create PNG info struct"); }
    Destructor _destroy_png([png_ptr, info_ptr] () { png_destroy_read_struct(&png_ptr, &info_ptr, NULL); });
    png_init_io(png_ptr, fp);
    png_set_sig_bytes(png_ptr, 7);

   // Read the PNG data into the internal data structure
   // PACKING causes 1/2/4-bit color samples to become 8-bit samples
   // EXPAND causes palette to become RGB and grayscale image 1/2/4-bit images to 8-bit, and resolves tRNS chunks
   png_read_png(png_ptr, info_ptr, PNG_TRANSFORM_PACKING | PNG_TRANSFORM_EXPAND, NULL);

   // Get information about the image data
   png_uint_32 h, w, stride = png_get_rowbytes(png_ptr, info_ptr);
   int bit_depth, color_type;
   png_get_IHDR(png_ptr, info_ptr, &w, &h, &bit_depth, &color_type, NULL, NULL, NULL);
   if (h > std::numeric_limits<int>::max() || w > std::numeric_limits<int>::max())
   {
       throw std::invalid_argument(filename + " is too large");
   }
   if (color_type == PNG_COLOR_TYPE_GRAY_ALPHA || color_type == PNG_COLOR_TYPE_RGB_ALPHA)
   {
       throw std::invalid_argument(filename + " uses alpha which is not supported");
   }
   width = w; height = h;

   // Get the conversion function
   typedef void (*convert) (png_bytep row, Type* out, size_t n);
   convert conv = NULL;
   if (color_type == PNG_COLOR_TYPE_GRAY)
   {
       if (bit_depth == 8) { conv = convert_im<uint8_t, Type>; }
       else if (bit_depth == 16) { conv = convert_im<uint16_t, Type>; }
   }
   else if (color_type == PNG_COLOR_TYPE_RGB && bit_depth == 8) { conv = convert_im<RGB24, Type>; }
   if (!conv) { throw std::invalid_argument(filename + " has unknown/unsupported color type & bit-depth"); }

   // Get the image data as an array of rows
   png_bytepp rows = png_get_rows(png_ptr, info_ptr);

   // Copy the data to a new chunk of memory
   std::shared_ptr<Type> out(new Type[h*w], std::default_delete<Type[]>());
   for (png_uint_32 i = 0; i < height; ++i) { conv(rows[i], out, w); }

   // Done (cleanup handled by RAII)
   return out;
}
