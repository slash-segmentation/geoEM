#pragma once

///////////////////////////////////////////////////////////////////////////////
// Low-level image conversion functions                                      //
///////////////////////////////////////////////////////////////////////////////

#include <exception>    // std::invalid_argument
#include <limits>       // std::numeric_limits
#include <type_traits>  // std::enable_if and std::is_* 
#include <string.h>     // memcpy
#include <stdint.h>     // (u)int(#)_t


// Utilities to detect type of template argument (floating point, signed integer, unsigned integer)
#define FP(T)   std::is_floating_point<T>::value
#define UINT(T) std::is_integral<T>::value && std::is_unsigned<T>::value
#define INT(T)  std::is_integral<T>::value && std::is_signed<T>::value


///////////////////////////////////////////////////////////////////////////////
// Convert single pixels between types                                       //
///////////////////////////////////////////////////////////////////////////////

template<class Type> constexpr Type uint_factor() { return std::numeric_limits<Type>::max(); }
template<class Type> constexpr Type int_factor() { return -std::numeric_limits<Type>::min(); }
template<class T1, class T2> constexpr T2 uint_uint_factor() { return uint_factor<T2>() / (T2)uint_factor<T1>(); }
template<class T1, class T2> constexpr T2 uint_uint_round() { return (T1)(uint_factor<T1>() >> 1) + (uint_factor<T2>() >> 1); }


// Same to same
template <class Src, class Dest>
typename std::enable_if<std::is_same<Src, Dest>::value, Dest>::type
inline convert_px(const Src x) { return x; }

// Floating-point to floating-point (just a cast)
template <class Src, class Dest>
typename std::enable_if<FP(Src) && FP(Dest), Dest>::type
inline convert_px(const Src x) { return (Dest)x; }

// Unsigned integer to floating-point (divide by max, result is 0.0 to 1.0)
template<class Src, class Dest>
typename std::enable_if<UINT(Src) && FP(Dest), Dest>::type
inline convert_px(const Src x) { return x / (Dest)uint_factor<Src>(); }

// Signed integer to floating-point (divide by -min, result is -1.0 to ~1.0)
template<class Src, class Dest>
typename std::enable_if<INT(Src) && FP(Dest), Dest>::type
inline convert_px(const Src x) { return x / (Dest)int_factor<Src>(); }

// Floating-point to unsigned integer (multiply by max)
template<class Src, class Dest>
typename std::enable_if<FP(Src) && UINT(Dest), Dest>::type
inline convert_px(const Src x) { return (Dest)(x * uint_factor<Src>()); }

// Floating-point to signed integer (multiply by -min)
template<class Src, class Dest>
typename std::enable_if<FP(Src) && INT(Dest), Dest>::type
inline convert_px(const Src x) { return (Dest)(x * int_factor<Src>()); }

// Unsigned integer to unsigned integer (uses shift)
template<class Src, class Dest>
typename std::enable_if<UINT(Src) && UINT(Dest), Dest>::type
inline convert_px(const Src x)
{
    if (sizeof(Src) > sizeof(Dest))
    {
        // First check if we can do it the fast way with a larger integer intermediate
        if ((sizeof(Src) + sizeof(Dest)) <= sizeof(uint32_t))
        {
            return (Dest)((((uint32_t)x) * uint_factor<Dest>() + uint_uint_round<Src, Dest>()) >> (8*sizeof(Src)));
        }
        else if ((sizeof(Src) + sizeof(Dest)) <= sizeof(uint64_t))
        {
            return (Dest)((((uint64_t)x) * uint_factor<Dest>() + uint_uint_round<Src, Dest>()) >> (8*sizeof(Src)));
        }
        // We can't - we have to do it the slow way that doesn't require a larger integer type
        else
        {
            // The high bits - most important, frequently these will be correct already
            Dest hi = (Dest)(x >> (sizeof(Src) - sizeof(Dest)));
            Dest hi_max = uint_factor<Dest>(); // max hi value
            Dest hi_half = ((Dest)1) << (sizeof(Dest) - 1); // half-way hi value
            // The low bits - used to determine rounding
            Src lo_max = (((Src)1) << (sizeof(Src) - sizeof(Dest))) - 1; // max lo value
            Src lo = x & lo_max;
            // The cutoff for rounding (takes advantage of integer overflow)
            Src cutoff = ((Src)(hi + hi_half)) * (lo_max / (Src)hi_max);
            // Return the hi plus any adjustments for rounding
            return hi + ((hi > hi_half) ? ((lo < cutoff) ? -1 : 0) : ((lo > cutoff) ? 1 : 0));
       }
    }
    else { return (Dest)(x * uint_uint_factor<Src, Dest>()); }
}

// TODO: Signed integer to signed integer
template<class Src, class Dest>
typename std::enable_if<INT(Src) && INT(Dest), Dest>::type
inline convert_px(const Src x) { throw std::invalid_argument("not implemented"); }

// TODO: Signed integer to unsigned integer
template<class Src, class Dest>
typename std::enable_if<INT(Src) && UINT(Dest), Dest>::type
inline convert_px(const Src x) { throw std::invalid_argument("not implemented"); }

// TODO: Unsigned integer to signed integer
template<class Src, class Dest>
typename std::enable_if<UINT(Src) && INT(Dest), Dest>::type
convert_px(const Src x) { throw std::invalid_argument("not implemented"); }


///////////////////////////////////////////////////////////////////////////////
// Convert multiple pixels between types                                     //
///////////////////////////////////////////////////////////////////////////////

// Specialization for same-to-same
template <class Src, class Dest>
typename std::enable_if<std::is_same<Src, Dest>::value>::type
inline convert_image_data(const Src* data, Dest* out, size_t n) { memcpy(out, data, n*sizeof(Dest)); }

// Fallback
template <class Src, class Dest>
typename std::enable_if<!std::is_same<Src, Dest>::value>::type
inline convert_image_data(const Src* data, Dest* out, size_t n)
{
    for (size_t i = 0; i < n; ++i) { out[i] = convert_px<Src, Dest>(data[i]); }
}

// Determine type from ImageIO raw image itself
template <class Dest>
inline void convert_image_data(const _image* im, Dest* out, size_t off, size_t n)
{
    if (im->wordKind == WK_FLOAT)
    {
        if (im->sign == SGN_SIGNED)
        {
            if (im->wdim == sizeof(float)) { convert_image_data(((float*)im->data)+off, out, n); return; }
            if (im->wdim == sizeof(double)) { convert_image_data(((double*)im->data)+off, out, n); return; }
            if (im->wdim == sizeof(long double)) { convert_image_data(((long double*)im->data)+off, out, n); return; }
        }
    }
    else if (im->wordKind == WK_FIXED)
    {
        if (im->sign == SGN_SIGNED)
        {
            if (im->wdim == sizeof(int8_t))  { convert_image_data(((int8_t *)im->data)+off, out, n); return; }
            if (im->wdim == sizeof(int16_t)) { convert_image_data(((int16_t*)im->data)+off, out, n); return; }
            if (im->wdim == sizeof(int32_t)) { convert_image_data(((int32_t*)im->data)+off, out, n); return; }
            if (im->wdim == sizeof(int64_t)) { convert_image_data(((int64_t*)im->data)+off, out, n); return; }
        }
        else if (im->sign == SGN_UNSIGNED)
        {
            if (im->wdim == sizeof(uint16_t)) { convert_image_data(((uint16_t*)im->data)+off, out, n); return; }
            if (im->wdim == sizeof(uint8_t))  { convert_image_data(((uint8_t *)im->data)+off, out, n); return; }
            if (im->wdim == sizeof(uint32_t)) { convert_image_data(((uint32_t*)im->data)+off, out, n); return; }
            if (im->wdim == sizeof(uint64_t)) { convert_image_data(((uint64_t*)im->data)+off, out, n); return; }
        }
    }
    throw std::invalid_argument("unknown image data type");
}

// Without offset or size
template <class Dest>
inline void convert_image_data(const _image* im, Dest* out)
{
    size_t n = ((size_t)im->xdim) * im->ydim * im->zdim * im->vdim;
    convert_image_data(im, out, 0, n);
}


///////////////////////////////////////////////////////////////////////////////
// Fill in image word type from C type                                       //
///////////////////////////////////////////////////////////////////////////////
template <class WordType>
void set_image_word_type(_image* im)
{
    im->wordKind = std::numeric_limits<WordType>::is_integer ? WK_FIXED : WK_FLOAT;
    im->sign = std::numeric_limits<WordType>::is_signed ? SGN_SIGNED : SGN_UNSIGNED;
    im->wdim = sizeof(WordType);
}


#undef FP
#undef UINT
#undef INT
