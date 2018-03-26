#pragma once

#include <exception>
#include <memory>
#include <limits>

typedef unsigned char byte;
typedef byte* bytes;

#pragma pack(push, 1)
struct RGB24 { byte R, G, B; };
#pragma pack(pop)


///////////////////////////////////////////////////////////////////////////////
// Convert single pixels between types                                       //
///////////////////////////////////////////////////////////////////////////////
#define FP(T)   !std::numeric_limits<T>::is_integer::value
#define UINT(T) std::numeric_limits<T>::is_integer::value && !std::numeric_limits<T>::is_signed::value
#define INT(T)  std::numeric_limits<T>::is_integer::value && std::numeric_limits<T>::is_signed::value
#define RGB24_(T) std::is_same<T, RGB24>::value

template <class Src, class Dest>
Dest convert_px(const Src x) { throw std::invalid_argument("not implemented"); }

// Same to same
template <class Src, class Dest>
typename std::enable_if<std::is_same<Src, Dest>::value, Dest>::type
convert_px(const Src x) { return x; }

// Floating-point to floating-point (just a cast)
template <class Src, class Dest>
typename std::enable_if<FP(Src) && FP(Dest), Dest>::type
convert_px(const Src x) { return (Dest)x; }

// Unsigned integer to floating-point (divide by max, result is 0.0 to 1.0)
template<class Src, class Dest>
typename std::enable_if<UINT(Src) && FP(Dest), Dest>::type
convert_px(const Src x) { return x / (Dest)std::numeric_limits<Src>::max(); }

// Signed integer to floating-point (divide by -min, result is -1.0 to ~1.0)
template<class Src, class Dest>
typename std::enable_if<INT(Src) && FP(Dest), Dest>::type
convert_px(const Src x) { return x / (Dest)-std::numeric_limits<Src>::min(); }

// Floating-point to unsigned integer (multiply by max)
template<class Src, class Dest>
typename std::enable_if<FP(Src) && UINT(Dest), Dest>::type
convert_px(const Src x) { return (Dest)(x * std::numeric_limits<Src>::max()); }

// Floating-point to signed integer (multiply by -min)
template<class Src, class Dest>
typename std::enable_if<FP(Src) && INT(Dest), Dest>::type
convert_px(const Src x) { return (Dest)(x * -std::numeric_limits<Src>::min()); }

// Unsigned integer to unsigned integer (uses shift)
template<class Src, class Dest>
typename std::enable_if<UINT(Src) && UINT(Dest), Dest>::type
convert_px(const Src x)
{
    if (sizeof(Src) > sizeof(Dest))
    {
        // TODO: rounding? (below code is not right though)
        //Dest round = (Dest)(x >> (8*(sizeof(Src) - sizeof(Dest)) - 1) & 0x1);
        return (Dest)(x >> (8*(sizeof(Src) - sizeof(Dest)))); // + round;
    }
    else
    {
        // TODO: fill in zero bits?
        return (Dest)(x << (8*(sizeof(Dest) - sizeof(Src))));
    }
}

// TODO: Signed integer to signed integer
template<class Src, class Dest>
typename std::enable_if<INT(Src) && INT(Dest), Dest>::type
convert_px(const Src x) { throw std::invalid_argument("not implemented"); }

// TODO: Signed integer to unsigned integer
template<class Src, class Dest>
typename std::enable_if<INT(Src) && UINT(Dest), Dest>::type
convert_px(const Src x) { throw std::invalid_argument("not implemented"); }

// TODO: Unsigned integer to signed integer
template<class Src, class Dest>
typename std::enable_if<UINT(Src) && INT(Dest), Dest>::type
convert_px(const Src x) { throw std::invalid_argument("not implemented"); }

// TODO: RGB24 to floating-point
template<class Src, class Dest>
typename std::enable_if<RGB24_(Src) && FP(Dest), Dest>::type
convert_px(const Src x) { throw std::invalid_argument("not implemented"); }

// TODO: RGB24 to unsigned integer
template<class Src, class Dest>
typename std::enable_if<RGB24_(Src) && UINT(Dest), Dest>::type
convert_px(const Src x) { throw std::invalid_argument("not implemented"); }

// TODO: RGB24 to signed integer
template<class Src, class Dest>
typename std::enable_if<RGB24_(Src) && INT(Dest), Dest>::type
convert_px(const Src x) { throw std::invalid_argument("not implemented"); }

// TODO: Floating-point to RGB24
template<class Src, class Dest>
typename std::enable_if<FP(Src) && RGB24_(Dest), Dest>::type
convert_px(const Src x) { throw std::invalid_argument("not implemented"); }

// TODO: Unsigned integer to RGB24
template<class Src, class Dest>
typename std::enable_if<UINT(Src) && RGB24_(Dest), Dest>::type
convert_px(const Src x) { throw std::invalid_argument("not implemented"); }

// TODO: Signed integer to RGB24
template<class Src, class Dest>
typename std::enable_if<INT(Src) && RGB24_(Dest), Dest>::type
convert_px(const Src x) { throw std::invalid_argument("not implemented"); }

#undef FP
#undef UINT
#undef INT
#undef RGB24_


///////////////////////////////////////////////////////////////////////////////
// Convert images between types                                              //
///////////////////////////////////////////////////////////////////////////////

template <class Src, class Dest>
void convert_im(const Src* data, Dest* out, size_t n)
{
    for (size_t i = 0; i < n; ++i) { out[i] = convert_px<Src, Dest>(data[i]); }
}

// Specialization for same-to-same
template <class Src, class Dest>
typename std::enable_if<std::is_same<Src, Dest>::value>::type
convert_im(const Src* data, Dest* out, size_t n) { memcpy(out, data, n*sizeof(Dest)); }

// With allocation
template <class Src, class Dest>
std::shared_ptr<Dest> convert_im(const Src* data, size_t n)
{
    std::shared_ptr<Dest> out(new Dest[n], std::default_delete<Dest[]>());
    convert_im(data, out, n);
    return out;
}
