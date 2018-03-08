#pragma once

#include <boost/config.hpp>

#include <stdexcept>
#include <iostream>
#include <utility>
#include <cstdint>
#include <cstring>
#include <cctype>
#include <string>
#include <array>

// TODO: one day make all flags names and None constexpr, but this can't be done until the C++ standard is clarified

///////////////////////////////////////////////////////////////////////////////
// Type-safe, bit-wise flags, for representing many true/false options in a
// compact variable. Uses template and macro magic to make them easy to use and
// safe. This uses many C++11 features including constexpr, variadic templates,
// explicit cast operators, and initializer lists. If any of these are missing,
// they are emulated with some side-effect. See the end for details.
//
// You create a flags type by using the FLAGS and FLAGS_WITH_VALUES macros:
//
//   FLAGS(TYPE_NAME, FLAG_NAME_1, FLAG_NAME_2, ...)
//   FLAGS_WITH_VALUES(TYPE_NAME, BASE_TYPE, NAME_1, VAL_1, NAME_2, VAL_2, ...)
//
// The FLAGS macro creates a new class named TYPE_NAME with static members with
// the given flag names, each having a unique power of two value. The base
// integral type is determined by the number of flags given.
//
// The FLAGS_WITH_VALUES macro allows you to specify the value for each flag.
// Like FLAGS, it creates a new class named TYPE_NAME. The base integral type
// is given by BASE_TYPE. The names and values for the flags are given
// alternating. Every flag must have a value after it.
//
// FLAGS() and FLAGS_WITH_VALUES() cannot be used inside a class directly or
// within a header without multiply-defined symbols. To do either of these
// you need to seperate declarations and defininitions, an example is:
//
//   class MyClass { FLAGS_DECL(MyFlags, A, B, C); }; // in header file
//   FLAGS_DEFN(MyClass::MyFlags, A, B, C);           // in source file
//
// The created flags class can then be explicitly constructed from a value of
// the integral base type, from another instance of the flags, or
// default-constructed with no flags set. All given flag names along with
// a special "None" name are available as static members of the class as well.
// Flags instances can also be can also be obtained using strings with their
// names through the from_name static methods.

// Flags instances behave much like you would expect them to, supporting the
// basic bit-wise operators. Additionally the flags can be converted to
// strings and directly output to ostreams.
//
// Flags have the following static members (all are const):
//   Flags *     - all of the named flags
//   Flags None  - Flags instance representing no flags set (a value of 0)
//
// Flags have the following member functions (all are const):
//   const char* const rawname() - the flag name, "None", or NULL if not a named flag
//   std::string       name()    - the flag name, "None", or a "|" joined list of named flags (e.g. "A|B")
//
// Flags have the following std::ostream adapter:
//  std::ostream& operator<< (std::ostream& s, const Flags& b) - writes the result of name() to the stream
//
// Flags have the following static functions:
//   Flags from_name(std::string name) - Gets a Flags from its string representation, supporting all return values from name() in addition to whitespace padding around each flag-name (e.g. "A | B")
//   Flags from_name(const char* name) - Same as above but with a null-terminated C-style string
//   Flags from_name(const char* name, size_t len) - Same as above but with a C-style string of a given length len
//   std::array<Flags,*> members()     - a collection of all named flags
//
// Flags support many operators (all non-compound are const and constexpr):
//   explicit operator base_type()       - Explicitly cast to the integral base type
//   base_type operator +()              - Get the integral base type
//   explicit operator bool()            - Explicitly cast to bool, true when at least one flag is set
//   bool operator!()                    - True when no flags are set
//   bool operator==(MyFlags)            - True when both MyFlags have the same set of flags set
//   bool operator!=(MyFlags)            - True when the MyFlags do not have the same flags set
//   MyFlags operator~()                 - Flips every bit that is covered by a flag
//   MyFlags operator&(MyFlags) [and &=] - Only keeps the flags set that both MyFlags have set
//   MyFlags operator|(MyFlags) [and |=] - Keeps the flags set that either MyFlags have set
//   MyFlags operator^(MyFlags) [and ^=] - Keeps the flags that are set in one but not the other
//
// The overhead of individual flags is minimal, each flags only contains a
// single variable of the integral base type. However, the class as a whole
// keeps a lot of static data but still it should be not so large and since
// it is not copied anywhere.
//
// Flags cannot have the following names (or any C++ keyword):
//   _mask, _count, _names, _values, _x,
//   super_type, base_type, type, name, rawname, from_name, members, None
//   (and _dummy/_safe_bool if explicit cast operators are not supported)
//
// The MSVC version also does not allow:
//   _fill_names, _fill_values, FT_, _count_ F#_, V#_
//
// These lists are not necessarily complete but should be close.
//
// Emulated C++11 Features Side-Effects:
//  variadic templates:      takes more run-time memory (up to 2 kb), does not
//                           detect illegal flag values (=0) during creation
//  explicit cast operators: cannot explicitly cast to base type (use +
//                           instead) and 'slow' and problematic bool casting
//  constexpr:               lacking some optimizations, minor though
///////////////////////////////////////////////////////////////////////////////

// Gets around a bug in MSVC when passing __VA_ARGS__ to another macro
#ifdef _MSC_VER
#define FLAGS___EXPAND(x) x
#else
#define FLAGS___EXPAND(...) __VA_ARGS__
#endif

#ifdef BOOST_NO_CXX11_CONSTEXPR
#define FLAGS___INIT_IF_NO_CONSTEXPR(...) = { __VA_ARGS__ }
#define FLAGS___INIT_IF_HAS_CONSTEXPR(...)
#else
#define FLAGS___INIT_IF_NO_CONSTEXPR(...)
#define FLAGS___INIT_IF_HAS_CONSTEXPR(...) = { __VA_ARGS__ }
#endif

namespace detail
{
#ifdef BOOST_NO_CXX11_VARIADIC_TEMPLATES
#define Values(B, A) \
    B V1  A, B V2  A, B V3  A, B V4  A, B V5  A, B V6  A, B V7  A, B V8  A, \
    B V9  A, B V10 A, B V11 A, B V12 A, B V13 A, B V14 A, B V15 A, B V16 A, \
    B V17 A, B V18 A, B V19 A, B V20 A, B V21 A, B V22 A, B V23 A, B V24 A, \
    B V25 A, B V26 A, B V27 A, B V28 A, B V29 A, B V30 A, B V31 A, B V32 A
#endif // BOOST_NO_CXX11_VARIADIC_TEMPLATES

#ifdef BOOST_NO_CXX11_VARIADIC_TEMPLATES
    template <class BaseType, class FlagsType, size_t Count, Values(BaseType, = 0)>
#else // BOOST_NO_CXX11_VARIADIC_TEMPLATES
    template <class BaseType, class FlagsType, size_t      , BaseType... Values>
#endif // BOOST_NO_CXX11_VARIADIC_TEMPLATES
    class flags_base
    {
    public:
        typedef BaseType base_type;
        typedef FlagsType type;
        static const type None; // no flags set

    protected:
        // Flag value (only non-static variable)
        base_type _x;

#ifdef BOOST_NO_CXX11_VARIADIC_TEMPLATES
        // Number of named flags
        static BOOST_CONSTEXPR const size_t _count = Count;

        // OR flag values together to get the mask of all flags
        static BOOST_CONSTEXPR const base_type _mask =
            V1  | V2  | V3  | V4  | V5  | V6  | V7  | V8  |
            V9  | V10 | V11 | V12 | V13 | V14 | V15 | V16 |
            V17 | V18 | V19 | V20 | V21 | V22 | V23 | V24 |
            V25 | V26 | V27 | V28 | V29 | V30 | V31 | V32 ;

        // Keep track of all flag values
        static BOOST_CONSTEXPR const base_type _values[32]
            FLAGS___INIT_IF_HAS_CONSTEXPR(Values(,));
#else // BOOST_NO_CXX11_VARIADIC_TEMPLATES
        // Number of named flags
        static BOOST_CONSTEXPR const size_t _count = sizeof...(Values);

        // Recursively OR flag values together to get the mask of all flags
        template<base_type... Vs> struct _or;
        template<base_type V>
        struct _or<V>
        {
            static_assert(V != 0, "No flag may have a value of 0");
            static BOOST_CONSTEXPR const base_type x = V;
        };
        template<base_type V, base_type... Vs>
        struct _or<V, Vs...>
        {
            static_assert(V != 0, "No flag may have a value of 0");
            static BOOST_CONSTEXPR const base_type x = V | _or<Vs...>::x;
        };
        static BOOST_CONSTEXPR const base_type _mask = _or<Values...>::x;

        // Keep track of all flag values
        static BOOST_CONSTEXPR const base_type _values[_count]
            FLAGS___INIT_IF_HAS_CONSTEXPR(Values...);
#endif // BOOST_NO_CXX11_VARIADIC_TEMPLATES
    public:
        static const std::array<type, _count>& members()
        {
#ifdef BOOST_NO_CXX11_VARIADIC_TEMPLATES
            static std::array<type, _count> m;
            if (m[0] == type::None)
            {
                for (size_t i = 0; i < type::_count; ++i)
                {
                    m[i] = type(type::_values[i]);
                }
            }
#else // BOOST_NO_CXX11_VARIADIC_TEMPLATES
            static const std::array<type, _count> m = { type(Values)... };
#endif // BOOST_NO_CXX11_VARIADIC_TEMPLATES
            return m;
        }

        ///// Constructors and casts /////
        inline BOOST_CONSTEXPR flags_base() : _x(0) { }
        inline BOOST_CONSTEXPR flags_base(const type& b) : _x(b._x) { }
        inline explicit flags_base(base_type b) : _x(b)
        {
            if ((b | type::_mask) != type::_mask)
            {
                throw std::domain_error("Flag value outside mask");
            }
        }
        inline BOOST_CONSTEXPR base_type operator+() const { return this->_x; }

#ifdef BOOST_NO_CXX11_EXPLICIT_CONVERSION_OPERATORS
    private:
        struct _dummy { void nonnull() {}; };
        typedef void (_dummy::*_safe_bool)();
    public:
        inline BOOST_CONSTEXPR operator _safe_bool() const
        {
            return this->_x != 0 ? &_dummy::nonnull : 0;
        }
        // no cast to base_type since it can't be explicit, use unary + operator instead
#else // BOOST_NO_CXX11_EXPLICIT_CONVERSION_OPERATORS
        inline BOOST_CONSTEXPR explicit operator bool() const { return _x != 0; }
        inline BOOST_CONSTEXPR explicit operator base_type() const { return _x; }
#endif // BOOST_NO_CXX11_EXPLICIT_CONVERSION_OPERATORS
        
        ///// Bit-wise operators /////
        inline BOOST_CONSTEXPR bool operator!() const { return this->_x == 0; }
        inline BOOST_CONSTEXPR bool operator==(const type b) const { return _x == b._x; }
        inline BOOST_CONSTEXPR bool operator!=(const type b) const { return _x != b._x; }
        inline BOOST_CONSTEXPR type operator ~() const { return type(~_x & type::_mask); }
        inline BOOST_CONSTEXPR type operator & (const type b) const { return type(_x & b._x); }
        inline BOOST_CONSTEXPR type operator | (const type b) const { return type(_x | b._x); }
        inline BOOST_CONSTEXPR type operator ^ (const type b) const { return type(_x ^ b._x); }
        inline type& operator &=(const type b) { _x &= b._x; return *static_cast<type*>(this); }
        inline type& operator |=(const type b) { _x |= b._x; return *static_cast<type*>(this); }
        inline type& operator ^=(const type b) { _x ^= b._x; return *static_cast<type*>(this); }

        ///// String/flag conversions /////

        // Returns "None", a predefined name, or NULL
        inline const char* const rawname() const
        {
            if (this->_x == 0) { return "None"; }
            for (size_t i = 0; i < _count; ++i)
            {
                if (this->_x == type::_values[i]) { return type::_names[i]; }
            }
            return NULL;
        }

        // Returns "None", a predefined name, or a string composed of predefined
        // names separated by a pipe "|".
        inline const std::string name() const
        {
            const char* const n = this->rawname();
            if (n) { return n; }
            std::string s = "";
            for (size_t i = 0; i < type::_count; ++i)
            {
                if ((this->_x & type::_values[i]) == type::_values[i])
                {
                    if (s.size()) { s += '|'; }
                    s += type::_names[i];
                }
            }
            return s;
        }
        friend inline std::ostream& operator<< (std::ostream& s, const type& b)
        {
            return s << b.name();
        }

        // Reverses "name()" accepting predefined names and "None" separated by
        // a pipe "|" while ignoring whitespace before and after each name -
        // throws std::invalid_argument if a name is not valid.
        static inline type from_name(const char *name, size_t len)
        {
            // Remove leading and trailing whitespace
            while (*name && std::isspace(*name)) { ++name; --len; }
            while (len>0 && std::isspace(name[len - 1])) { --len; }
            const char *pipe = (const char*)std::memchr(name, '|', len);
            if (pipe == NULL)
            {
                // No "pipe", we have a single name
                if (len == 4 && !std::strncmp(name, "None", 4)) { return type(0); }
                for (size_t i = 0; i < type::_count; ++i)
                {
                    const char* n = type::_names[i];
                    if (len == std::strlen(n) && !std::strncmp(name, n, len))
                    {
                        return type(type::_values[i]);
                    }
                }
                throw std::invalid_argument(
                    "Flags name '" + std::string(name, len) + "' not found");
            }
            // There is a pipe, split string up at pipes and process each part
            type b = type::from_name(name, pipe - name);
            len -= pipe - name;
            name = pipe + 1;
            while ((pipe = (const char*)std::memchr(name, '|', len)) != NULL)
            {
                b |= type::from_name(name, pipe - name);
                len -= pipe - name;
                name = pipe + 1;
            }
            b |= type::from_name(name, len - 1);
            return b;
        }
        static inline type from_name(const char * name)
        {
            return type::from_name(name, std::strlen(name));
        }
        static inline type from_name(const std::string name)
        {
            return type::from_name(name.c_str(), name.size());
        }
    };


#ifdef BOOST_NO_CXX11_VARIADIC_TEMPLATES
    // Define values
    template <class B, class F, size_t C, Values(B,)>
    BOOST_CONSTEXPR const B flags_base<B,F,C,Values(,)>::_values[32]
        FLAGS___INIT_IF_NO_CONSTEXPR(Values(,));

    // Define the None flag
    template <class B, class F, size_t C, Values(B, )> const F flags_base<B,F,C,Values(,)>::None;

#undef Values
#else // BOOST_NO_CXX11_VARIADIC_TEMPLATES
    // Define values
    template <class B, class F, size_t C, B... V>
    BOOST_CONSTEXPR const B flags_base<B,F,C,V...>::_values[flags_base<B,F,C,V...>::_count]
        FLAGS___INIT_IF_NO_CONSTEXPR(V...);

    // Define the None flag
    template <class B, class F, size_t C, B... V> const F flags_base<B,F,C,V...>::None;
#endif // BOOST_NO_CXX11_VARIADIC_TEMPLATES
}

////////// Helper macros //////////
// GET_ARG_32(0, __VA_ARGS__, ...)  Gets the 32nd argument from __VA_ARGS__+...
#define FLAGS___GET_ARG_32(X,X32,X31,X30,X29,X28,X27,X26,X25,X24,X23,X22,X21,X20,X19,X18,X17,X16,X15,X14,X13,X12,X11,X10,X9,X8,X7,X6,X5,X4,X3,X2,X1,N,...) N
// GET_ARG_64(0, __VA_ARGS__, ...)  Gets the 64th argument from __VA_ARGS__+...
#define FLAGS___GET_ARG_64(X,X64,X63,X62,X61,X60,X59,X58,X57,X56,X55,X54,X53,X52,X51,X50,X49,X48,X47,X46,X45,X44,X43,X42,X41,X40,X39,X38,X37,X36,X35,X34,X33,X32,X31,X30,X29,X28,X27,X26,X25,X24,X23,X22,X21,X20,X19,X18,X17,X16,X15,X14,X13,X12,X11,X10,X9,X8,X7,X6,X5,X4,X3,X2,X1,N,...) N
// NUM_ARGS(...) evaluates to the number of the passed-in arguments
#define FLAGS___NUM_ARGS(...) \
    FLAGS___EXPAND(FLAGS___GET_ARG_64(0,__VA_ARGS__, \
        64,63,62,61,60,59,58,57,56,55,54,53,52,51,50,49, \
        48,47,46,45,44,43,42,41,40,39,38,37,36,35,34,33, \
        32,31,30,29,28,27,26,25,24,23,22,21,20,19,18,17, \
        16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0))
// BASE_TYPE(...) evaluates to the base type that will fit the number of flags given by the var args
#define FLAGS___BASE_TYPE(...) \
    FLAGS___EXPAND(FLAGS___GET_ARG_32(0,__VA_ARGS__, \
        ::std::uint32_t,::std::uint32_t,::std::uint32_t,::std::uint32_t, \
        ::std::uint32_t,::std::uint32_t,::std::uint32_t,::std::uint32_t, \
        ::std::uint32_t,::std::uint32_t,::std::uint32_t,::std::uint32_t, \
        ::std::uint32_t,::std::uint32_t,::std::uint32_t,::std::uint32_t, \
        ::std::uint16_t,::std::uint16_t,::std::uint16_t,::std::uint16_t, \
        ::std::uint16_t,::std::uint16_t,::std::uint16_t,::std::uint16_t, \
        ::std::uint8_t, ::std::uint8_t, ::std::uint8_t, ::std::uint8_t, \
        ::std::uint8_t, ::std::uint8_t, ::std::uint8_t, ::std::uint8_t)) /* TODO: last one as bool? */
// POW2(...) evaluates to the maximum power of 2 that will be given to a flag for the given var args
#define FLAGS___POW2(...) \
    FLAGS___EXPAND(FLAGS___GET_ARG_32(0,__VA_ARGS__, \
        0x80000000, 0x40000000, 0x20000000, 0x10000000, \
        0x08000000, 0x04000000, 0x02000000, 0x01000000, \
        0x00800000, 0x00400000, 0x00200000, 0x00100000, \
        0x00080000, 0x00040000, 0x00020000, 0x00010000, \
        0x00008000, 0x00004000, 0x00002000, 0x00001000, \
        0x00000800, 0x00000400, 0x00000200, 0x00000100, \
        0x00000080, 0x00000040, 0x00000020, 0x00000010, \
        0x00000008, 0x00000004, 0x00000002, 0x00000001))

////////// NVLIST: expands a list of name-values using an X macro //////////
// X is a macro that takes the NAME of the class, and a flag name-value pair
// T is the type of the flags class
// S is a separator between items (once evaluated with ())
// N is the length of the list
#define FLAGS___NVLIST(X,T,S,N,...) FLAGS___EXPAND(FLAGS___NVL##N(X,T,S,__VA_ARGS__))
#define FLAGS___NVL2(X,T,S,N,V) X(T,N,V)
#define FLAGS___NVL4(X,T,S,N,V,...) X(T,N,V)S()FLAGS___EXPAND(FLAGS___NVL2(X,T,S,__VA_ARGS__))
#define FLAGS___NVL6(X,T,S,N,V,...) X(T,N,V)S()FLAGS___EXPAND(FLAGS___NVL4(X,T,S,__VA_ARGS__))
#define FLAGS___NVL8(X,T,S,N,V,...) X(T,N,V)S()FLAGS___EXPAND(FLAGS___NVL6(X,T,S,__VA_ARGS__))
#define FLAGS___NVL10(X,T,S,N,V,...) X(T,N,V)S()FLAGS___EXPAND(FLAGS___NVL8(X,T,S,__VA_ARGS__))
#define FLAGS___NVL12(X,T,S,N,V,...) X(T,N,V)S()FLAGS___EXPAND(FLAGS___NVL10(X,T,S,__VA_ARGS__))
#define FLAGS___NVL14(X,T,S,N,V,...) X(T,N,V)S()FLAGS___EXPAND(FLAGS___NVL12(X,T,S,__VA_ARGS__))
#define FLAGS___NVL16(X,T,S,N,V,...) X(T,N,V)S()FLAGS___EXPAND(FLAGS___NVL14(X,T,S,__VA_ARGS__))
#define FLAGS___NVL18(X,T,S,N,V,...) X(T,N,V)S()FLAGS___EXPAND(FLAGS___NVL16(X,T,S,__VA_ARGS__))
#define FLAGS___NVL20(X,T,S,N,V,...) X(T,N,V)S()FLAGS___EXPAND(FLAGS___NVL18(X,T,S,__VA_ARGS__))
#define FLAGS___NVL22(X,T,S,N,V,...) X(T,N,V)S()FLAGS___EXPAND(FLAGS___NVL20(X,T,S,__VA_ARGS__))
#define FLAGS___NVL24(X,T,S,N,V,...) X(T,N,V)S()FLAGS___EXPAND(FLAGS___NVL22(X,T,S,__VA_ARGS__))
#define FLAGS___NVL26(X,T,S,N,V,...) X(T,N,V)S()FLAGS___EXPAND(FLAGS___NVL24(X,T,S,__VA_ARGS__))
#define FLAGS___NVL28(X,T,S,N,V,...) X(T,N,V)S()FLAGS___EXPAND(FLAGS___NVL26(X,T,S,__VA_ARGS__))
#define FLAGS___NVL30(X,T,S,N,V,...) X(T,N,V)S()FLAGS___EXPAND(FLAGS___NVL28(X,T,S,__VA_ARGS__))
#define FLAGS___NVL32(X,T,S,N,V,...) X(T,N,V)S()FLAGS___EXPAND(FLAGS___NVL30(X,T,S,__VA_ARGS__))
#define FLAGS___NVL34(X,T,S,N,V,...) X(T,N,V)S()FLAGS___EXPAND(FLAGS___NVL32(X,T,S,__VA_ARGS__))
#define FLAGS___NVL36(X,T,S,N,V,...) X(T,N,V)S()FLAGS___EXPAND(FLAGS___NVL34(X,T,S,__VA_ARGS__))
#define FLAGS___NVL38(X,T,S,N,V,...) X(T,N,V)S()FLAGS___EXPAND(FLAGS___NVL36(X,T,S,__VA_ARGS__))
#define FLAGS___NVL40(X,T,S,N,V,...) X(T,N,V)S()FLAGS___EXPAND(FLAGS___NVL38(X,T,S,__VA_ARGS__))
#define FLAGS___NVL42(X,T,S,N,V,...) X(T,N,V)S()FLAGS___EXPAND(FLAGS___NVL40(X,T,S,__VA_ARGS__))
#define FLAGS___NVL44(X,T,S,N,V,...) X(T,N,V)S()FLAGS___EXPAND(FLAGS___NVL42(X,T,S,__VA_ARGS__))
#define FLAGS___NVL46(X,T,S,N,V,...) X(T,N,V)S()FLAGS___EXPAND(FLAGS___NVL44(X,T,S,__VA_ARGS__))
#define FLAGS___NVL48(X,T,S,N,V,...) X(T,N,V)S()FLAGS___EXPAND(FLAGS___NVL46(X,T,S,__VA_ARGS__))
#define FLAGS___NVL50(X,T,S,N,V,...) X(T,N,V)S()FLAGS___EXPAND(FLAGS___NVL48(X,T,S,__VA_ARGS__))
#define FLAGS___NVL52(X,T,S,N,V,...) X(T,N,V)S()FLAGS___EXPAND(FLAGS___NVL50(X,T,S,__VA_ARGS__))
#define FLAGS___NVL54(X,T,S,N,V,...) X(T,N,V)S()FLAGS___EXPAND(FLAGS___NVL52(X,T,S,__VA_ARGS__))
#define FLAGS___NVL56(X,T,S,N,V,...) X(T,N,V)S()FLAGS___EXPAND(FLAGS___NVL54(X,T,S,__VA_ARGS__))
#define FLAGS___NVL58(X,T,S,N,V,...) X(T,N,V)S()FLAGS___EXPAND(FLAGS___NVL56(X,T,S,__VA_ARGS__))
#define FLAGS___NVL60(X,T,S,N,V,...) X(T,N,V)S()FLAGS___EXPAND(FLAGS___NVL58(X,T,S,__VA_ARGS__))
#define FLAGS___NVL62(X,T,S,N,V,...) X(T,N,V)S()FLAGS___EXPAND(FLAGS___NVL60(X,T,S,__VA_ARGS__))
#define FLAGS___NVL64(X,T,S,N,V,...) X(T,N,V)S()FLAGS___EXPAND(FLAGS___NVL62(X,T,S,__VA_ARGS__))

// Separators for NVLIST
#define FLAGS___COMMA() ,
#define FLAGS___SEMICOLON() ;

// X-macros for NVLIST
#define FLAGS___XNAME(T, N, V) #N
#define FLAGS___XVALUE(T, N, V) V
#define FLAGS___XDECLARE(T, N, V) static const type N
#define FLAGS___XDEFINE(T, N, V) template<> const T T::N(V)

////////// CREATE: Create the actual flags class //////////
// T is the typename of the new flags-class
// BT is the base typename for the flags
// N is the number of flags (== NUM_ARGS(__VA_ARGS__))
// First stage of creation expands __VA_ARGS__ into all necessary forms
#define FLAGS___CREATE_DECL(T, BT, N, ...) \
    FLAGS___CREATE2_DECL(T, BT, N, \
        FLAGS___NVLIST(FLAGS___XNAME, T, FLAGS___COMMA, N, __VA_ARGS__), \
        FLAGS___NVLIST(FLAGS___XVALUE, T, FLAGS___COMMA, N, __VA_ARGS__), \
        FLAGS___NVLIST(FLAGS___XDECLARE, T, FLAGS___SEMICOLON, N, __VA_ARGS__) \
    )
#define FLAGS___CREATE_DEFN(T, BT, N, ...) \
    FLAGS___CREATE2_DEFN(T, BT, N, \
        FLAGS___NVLIST(FLAGS___XNAME, T, FLAGS___COMMA, N, __VA_ARGS__), \
        FLAGS___NVLIST(FLAGS___XDEFINE, T, FLAGS___SEMICOLON, N, __VA_ARGS__) \
    )
// Second stage no longer needs __VA_ARGS__
#define FLAGS___CREATE2_DECL(T, BT, N, NAMES, VALS, DECS) \
    class T : public ::detail::flags_base<BT, T, N / 2, VALS> \
    { \
        typedef ::detail::flags_base<BT, T, N / 2, VALS> super_type; \
        friend super_type; \
        static_assert(std::is_integral<BT>::value, \
            "Flags base type must be integral"); \
        static_assert(N != 0, \
            "Must have at least one flag name"); \
        static_assert(((N&1)==0) && (super_type::_count*2==N), \
            "Number of flag names and values must be the same"); \
        static BOOST_CONSTEXPR const char* const _names[super_type::_count] \
            FLAGS___INIT_IF_HAS_CONSTEXPR(NAMES); \
    public: \
        typedef BT base_type; \
        typedef T type; \
        inline BOOST_CONSTEXPR T() : super_type() { } \
        inline BOOST_CONSTEXPR T(const type& b) : super_type(b) { } \
        inline explicit T(base_type b) : super_type(b) { } \
        inline type& operator=(const type& b) { this->_x = b._x; return *this; } \
        DECS; \
    };
#define FLAGS___CREATE2_DEFN(T, BT, N, NAMES, DEFS) \
    template<> \
    BOOST_CONSTEXPR const char* const T::_names [T::_count] \
        FLAGS___INIT_IF_NO_CONSTEXPR(NAMES); \
    DEFS


////////// ADDVAL: Adds powers of two between each flag name //////////
#define FLAGS___ADDVAL(MAX,N,...) FLAGS___EXPAND(FLAGS___ADDVAL_(MAX,N,__VA_ARGS__))
#define FLAGS___ADDVAL_(M,N,...) FLAGS___EXPAND(FLAGS___ADDVAL##N(M,__VA_ARGS__)) /* indirection layer */
#define FLAGS___ADDVAL1(M,N) N,M
#define FLAGS___ADDVAL2(M,N,...) N,(M>>1),FLAGS___EXPAND(FLAGS___ADDVAL1(M,__VA_ARGS__))
#define FLAGS___ADDVAL3(M,N,...) N,(M>>2),FLAGS___EXPAND(FLAGS___ADDVAL2(M,__VA_ARGS__))
#define FLAGS___ADDVAL4(M,N,...) N,(M>>3),FLAGS___EXPAND(FLAGS___ADDVAL3(M,__VA_ARGS__))
#define FLAGS___ADDVAL5(M,N,...) N,(M>>4),FLAGS___EXPAND(FLAGS___ADDVAL4(M,__VA_ARGS__))
#define FLAGS___ADDVAL6(M,N,...) N,(M>>5),FLAGS___EXPAND(FLAGS___ADDVAL5(M,__VA_ARGS__))
#define FLAGS___ADDVAL7(M,N,...) N,(M>>6),FLAGS___EXPAND(FLAGS___ADDVAL6(M,__VA_ARGS__))
#define FLAGS___ADDVAL8(M,N,...) N,(M>>7),FLAGS___EXPAND(FLAGS___ADDVAL7(M,__VA_ARGS__))
#define FLAGS___ADDVAL9(M,N,...) N,(M>>8),FLAGS___EXPAND(FLAGS___ADDVAL8(M,__VA_ARGS__))
#define FLAGS___ADDVAL10(M,N,...) N,(M>>9),FLAGS___EXPAND(FLAGS___ADDVAL9(M,__VA_ARGS__))
#define FLAGS___ADDVAL11(M,N,...) N,(M>>10),FLAGS___EXPAND(FLAGS___ADDVAL10(M,__VA_ARGS__))
#define FLAGS___ADDVAL12(M,N,...) N,(M>>11),FLAGS___EXPAND(FLAGS___ADDVAL11(M,__VA_ARGS__))
#define FLAGS___ADDVAL13(M,N,...) N,(M>>12),FLAGS___EXPAND(FLAGS___ADDVAL12(M,__VA_ARGS__))
#define FLAGS___ADDVAL14(M,N,...) N,(M>>13),FLAGS___EXPAND(FLAGS___ADDVAL13(M,__VA_ARGS__))
#define FLAGS___ADDVAL15(M,N,...) N,(M>>14),FLAGS___EXPAND(FLAGS___ADDVAL14(M,__VA_ARGS__))
#define FLAGS___ADDVAL16(M,N,...) N,(M>>15),FLAGS___EXPAND(FLAGS___ADDVAL15(M,__VA_ARGS__))
#define FLAGS___ADDVAL17(M,N,...) N,(M>>16),FLAGS___EXPAND(FLAGS___ADDVAL16(M,__VA_ARGS__))
#define FLAGS___ADDVAL18(M,N,...) N,(M>>17),FLAGS___EXPAND(FLAGS___ADDVAL17(M,__VA_ARGS__))
#define FLAGS___ADDVAL19(M,N,...) N,(M>>18),FLAGS___EXPAND(FLAGS___ADDVAL18(M,__VA_ARGS__))
#define FLAGS___ADDVAL20(M,N,...) N,(M>>19),FLAGS___EXPAND(FLAGS___ADDVAL19(M,__VA_ARGS__))
#define FLAGS___ADDVAL21(M,N,...) N,(M>>20),FLAGS___EXPAND(FLAGS___ADDVAL20(M,__VA_ARGS__))
#define FLAGS___ADDVAL22(M,N,...) N,(M>>21),FLAGS___EXPAND(FLAGS___ADDVAL21(M,__VA_ARGS__))
#define FLAGS___ADDVAL23(M,N,...) N,(M>>22),FLAGS___EXPAND(FLAGS___ADDVAL22(M,__VA_ARGS__))
#define FLAGS___ADDVAL24(M,N,...) N,(M>>23),FLAGS___EXPAND(FLAGS___ADDVAL23(M,__VA_ARGS__))
#define FLAGS___ADDVAL25(M,N,...) N,(M>>24),FLAGS___EXPAND(FLAGS___ADDVAL24(M,__VA_ARGS__))
#define FLAGS___ADDVAL26(M,N,...) N,(M>>25),FLAGS___EXPAND(FLAGS___ADDVAL25(M,__VA_ARGS__))
#define FLAGS___ADDVAL27(M,N,...) N,(M>>26),FLAGS___EXPAND(FLAGS___ADDVAL26(M,__VA_ARGS__))
#define FLAGS___ADDVAL28(M,N,...) N,(M>>27),FLAGS___EXPAND(FLAGS___ADDVAL27(M,__VA_ARGS__))
#define FLAGS___ADDVAL29(M,N,...) N,(M>>28),FLAGS___EXPAND(FLAGS___ADDVAL28(M,__VA_ARGS__))
#define FLAGS___ADDVAL30(M,N,...) N,(M>>29),FLAGS___EXPAND(FLAGS___ADDVAL29(M,__VA_ARGS__))
#define FLAGS___ADDVAL31(M,N,...) N,(M>>30),FLAGS___EXPAND(FLAGS___ADDVAL30(M,__VA_ARGS__))
#define FLAGS___ADDVAL32(M,N,...) N,(M>>31),FLAGS___EXPAND(FLAGS___ADDVAL31(M,__VA_ARGS__))
#define FLAGS___ADDVALS_AND_CREATE(X, TYPE_NAME, ...) \
    FLAGS___EXPAND(X( \
        TYPE_NAME, \
        FLAGS___BASE_TYPE(__VA_ARGS__), \
        FLAGS___ADDVAL(\
            FLAGS___POW2(__VA_ARGS__), \
            FLAGS___NUM_ARGS(__VA_ARGS__), \
            __VA_ARGS__ \
        ) \
    ))


////////// Create a set of flags with the given values //////////
#define FLAGS_WITH_VALUES(TYPE_NAME, BASE_TYPE, ...) \
    FLAGS___CREATE_DECL( \
        TYPE_NAME, BASE_TYPE, FLAGS___NUM_ARGS(__VA_ARGS__), __VA_ARGS__ \
    ) \
    FLAGS___CREATE_DEFN( \
        TYPE_NAME, BASE_TYPE, FLAGS___NUM_ARGS(__VA_ARGS__), __VA_ARGS__ \
    )
#define FLAGS_WITH_VALUES_DECL(TYPE_NAME, BASE_TYPE, ...) \
    FLAGS___CREATE_DECL( \
        TYPE_NAME, BASE_TYPE, FLAGS___NUM_ARGS(__VA_ARGS__), __VA_ARGS__ \
    )
#define FLAGS_WITH_VALUES_DEFN(TYPE_NAME, BASE_TYPE, ...) \
    FLAGS___CREATE_DEFN( \
        TYPE_NAME, BASE_TYPE, FLAGS___NUM_ARGS(__VA_ARGS__), __VA_ARGS__ \
    )


////////// Create a set of flags with default powers-of-two values //////////
#define FLAGS(TYPE_NAME, ...) \
    FLAGS___ADDVALS_AND_CREATE(FLAGS_WITH_VALUES, TYPE_NAME, __VA_ARGS__)
#define FLAGS_DECL(TYPE_NAME, ...) \
    FLAGS___ADDVALS_AND_CREATE(FLAGS_WITH_VALUES_DECL, TYPE_NAME, __VA_ARGS__)
#define FLAGS_DEFN(TYPE_NAME, ...) \
    FLAGS___ADDVALS_AND_CREATE(FLAGS_WITH_VALUES_DEFN, TYPE_NAME, __VA_ARGS__)
