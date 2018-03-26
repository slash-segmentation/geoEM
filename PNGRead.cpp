#include "PNGRead.hpp"

bool check_if_png(FILE *fp)
{
   unsigned char buf[7];
   if (fread(buf, 1, 7, fp) != 7) { return false; }
   return !png_sig_cmp(buf, 0, 7);
}

void png_cpp_error_fn(png_structp png_ptr, png_const_charp error_msg) { throw std::invalid_argument(error_msg); }
