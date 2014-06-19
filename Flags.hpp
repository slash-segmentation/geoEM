#pragma once

#include <type_traits>
#include <stdexcept>
#include <iostream>
#include <utility>
#include <cstdint>
#include <cstring>
#include <cctype>
#include <string>
#include <array>

///////////////////////////////////////////////////////////////////////////////////////////////////
// Type-safe, bit-wise flags, for representing many true/false options in a compact variable.
// Uses template and macro magic to make them easy to use and safe. The full version requires many
// C++11 features including constexpr, variadic templates, and initializer lists. However there is 
// a version explicitly for MSVC which is API-compatible but only supports 8 flags (could be
// expanded). It has some other caveats as well, see the bottom of this description for more info.
//
// The easiest way to create a set of flags is to use the FLAGS() and FLAGS_WITH_VALUES() at the
// end of this document. However, you can easily do this on your own as well.
//
// To create a set of flags you first need to create flag names using the macro FLAG_NAME(Name).
// This macro creates a template class called Flag_Name in the current namespace. This class
// should not be used except as a parameter to Flags. You then typedef the Flags template class
// with a base type and the flag name classes and then choose a set of values for the flags.
//
// Example creation:
//   FLAG_NAME(A);
//   FLAG_NAME(B);
//   FLAG_NAME(C);
//   typedef Flags<std::uint8_t, Flag_A, Flag_B, Flag_C>::WithValues<0x01, 0x02, 0x04> MyFlags;
// In this case we could have used WithDDefaultValues instead of WithValues<...>.
//
// FLAG_NAME cannot be used inside a class. To define flag names inside of a class you must split
// the FLAG_NAME macro into FLAG_NAME_CORE and FLAG_NAME_INIT_NS, doing something like the
// following:
//   class MyClass
//   {
//     FLAG_NAME_CORE(A);
//     FLAG_NAME_CORE(B);
//     FLAG_NAME_CORE(C);
//   public:
//     typedef Flags<std::uint8_t, Flag_A, Flag_B, Flag_C>::WithDefaultValues MyFlags;
//   };
//   FLAG_NAME_INIT_NS(MyClass, A);
//   FLAG_NAME_INIT_NS(MyClass, B);
//   FLAG_NAME_INIT_NS(MyClass, C);
//
// MyFlags object can then be explicitly constructed from a value of the integral base type, from
// another MyFlags, or default-constructed with no flags set. All pre-defined names along with
// a special "None" name are available as static members of the class as well. They can also be
// queried with strings representing them through the from_name static methods. They behave much
// like you would expect a flag object to, supporting the basic bot-wise operators. Additionally
// the flags can be directly output to ostreams.
//
// MyFlags has the following member variables (all are const and constexpr):
//   MyFlags A, MyFlags B, MyFlags C  - all of the named flags as MyFlags objects
//   MyFlags None                     - MyFlags object representing no flags set (a value of 0)
//
// MyFlags has the following member functions (all are const):
//   const char* const rawname()    - the pre-defined flag name (e.g. "A"), "None", or NULL if not a pre-defined flag
//   std::string       name()       - the flag name, if not a pre-defined flag gives a "|" joined list of pre-defined flag names (e.g. "A|B")
//
// MyFlags has the following std::ostream adapter:
//  std::ostream& operator<< (std::ostream& s, const type& b) - writes the result of name() to the stream
//
// MyFlags has the following static functions:
//   MyFlags from_name(std::string name) - Gets the MyFlags from its string representation, supporting all return values from name() in addition to whitespace padding around each flag-name (e.g. "A | B")
//   MyFlags from_name(const char* name) - Same as above but with a null-terminated C-style string.
//   MyFlags from_name(const char* name, size_t n) - Same as above but with a C-style string of a given length n.
//   std::array<MyFlags,3> members()     - a collection of all MyFlags pre-defined objects (A, B, and C)
//
// MyFlags operators (all non-compound are const and constexpr):
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
// The overhead of individual flags objects should be minimal. Every flags object only contains a
// single variable of the integral base type. However, the class as a whole keeps a lot of static
// data but still it should be not so large and since it is not copied anywhere, it doesn't really
// matter (and in fact, I hope the compiler actually removes a lot of it - in fact if you never use
// name()/from_name()/members it should be able to remove all but one integral base type and the
// static versions of the pre-defined flags).
//
// 
// Flags cannot have the following names:
//   F_, Fs_, BT_, V_, Vs_, _name, _or, _mask, _count, _names, _values, _x, base_type, type, name, rawname, from_name, None, members
//   (and _dummy/_safe_bool if the compiler does not support explicit cast operators)
//
// The MSVC version also does not allow:
//   _fill_names, _fill_values, FT_, _count_ F#_, V#_
//
// These lists are not necessarily complete but should be close.
//
// The MSCV version is limited to 8 flags (could be expanded with copying and pasting) and slightly
// worse type safety. Additionally, it is really ugly code.
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef _MSC_VER
#include "Flags-MSVC.hpp"
#else
// Creates a name to use as a flag
#define FLAG_NAME(NAME) FLAG_NAME_CORE(NAME) FLAG_NAME_INIT(NAME)
#define FLAG_NAME_CORE(NAME) \
	template <class BT_, class F_, BT_ V_> struct Flag_##NAME \
	{ \
		static_assert(V_ != 0, "No flag may have a value of 0 (use the special None flag for that)"); \
		static constexpr const F_ NAME = F_(std::integral_constant<BT_,V_>()); \
		protected: static constexpr const char* const _name = #NAME; \
	};
#define FLAG_NAME_INIT(NAME) \
	template<class BT_, class F_, BT_ V_> constexpr const F_ Flag_##NAME<BT_,F_,V_>::NAME /* make available at runtime (only needed when flag names are used inline) */
#define FLAG_NAME_INIT_NS(NS, NAME) \
	template<class BT_, class F_, BT_ V_> constexpr const F_ NS::Flag_##NAME<BT_,F_,V_>::NAME /* make available at runtime (only needed when flag names are used inline) */

template<class BT_, template <class,class,BT_> class... Fs_>
class Flags
{
	typedef BT_ base_type;

	// F_ is a dummy class which can be given to a flag-name template
	// This is useful in cases where we don't have values to give the flag-name template or [for MSVC] when we want to simplify parameter packing
	struct F_ { typedef BT_ base_type; template <base_type V_> inline constexpr explicit F_(std::integral_constant<base_type, V_>) { } };
	// Count the number of flags
#ifdef __clang__
	static constexpr const size_t _count = sizeof...(Fs_); //(Fs_<base_type, F_, 1>);
#else
	static constexpr const size_t _count = sizeof...(Fs_<base_type, F_, 1>);
#endif
	static_assert(_count != 0, "Must have at least one flag name");

public:
	template<BT_... Vs_>
	class WithValues : public Fs_<BT_, Flags<BT_, Fs_...>::WithValues<Vs_...>, Vs_>...
	{
	public:
		typedef BT_ base_type;
		typedef Flags<base_type, Fs_...>::WithValues<Vs_...> type;

		static_assert(std::is_integral<base_type>::value, "Flags base type must be an integral type");
#ifdef __clang__
		static_assert(sizeof...(Fs_) == sizeof...(Vs_), "Number of flag names and values must be the same");
#else
		static_assert(sizeof...(Fs_<base_type, type, Vs_>) == sizeof...(Vs_), "Number of flag names and values must be the same");
#endif

		// Special flags
		static constexpr const type None = type(std::integral_constant<base_type, 0          >()); // no flags set
		//static constexpr const type All  = type(std::integral_constant<base_type, type::_mask>()); // all flags set (could do ~None as well)

	private:
		// Flag value (only non-static variable)
		base_type _x;

		// Recursively OR flag values together to get the mask of all flags
		template<base_type... Vs> struct _or;
		template<base_type V>                  struct _or<V>       { static constexpr const base_type value = V; };
		template<base_type V, base_type... Vs> struct _or<V,Vs...> { static constexpr const base_type value = V | _or<Vs...>::value; };
		static constexpr const base_type _mask = _or<Vs_...>::value;

		// Other basic properties about the collection of flags
		static constexpr const char* const _names [_count] = { Fs_<base_type, type, Vs_>::_name... };
		static constexpr const base_type   _values[_count] = { Vs_... };

	public:
		static /*constexpr*/ const std::array<type, _count>& members() { static const std::array<type, _count> m = { Flags<BT_,Fs_...>::WithValues<Vs_...>(Vs_)... }; return m; }

		// Constructors, assignment operators, and casts
		inline constexpr WithValues() : _x(0) { }
		inline constexpr WithValues(const type& b) : _x(b._x) { }
		inline type& operator=(const type& b) { this->_x = b._x; return *this; }
	
		inline explicit WithValues(base_type b) : _x(b) { if ((b|type::_mask)!=type::_mask) { throw std::domain_error("Flag value outside mask"); } } // this is not constexpr to support throwing an exception, use Flags(std::integral_constant<base_type, b>()) instead for constexpr support
		template <base_type V_>
		inline constexpr explicit WithValues(std::integral_constant<base_type, V_>) : _x(V_) { static_assert((V_ | type::_mask) == type::_mask, "Flag value outside mask"); }
		//inline type& operator=(base_type b) { if ((b|type::_mask)!=type::_mask) { throw std::domain_error("Flag value outside mask"); } this->_x = b; return *this; }
		inline constexpr base_type operator+() const { return this->_x; }
	
		inline constexpr explicit operator bool()      const { return this->_x != 0; }
		inline constexpr explicit operator base_type() const { return this->_x; }

		// Bit-wise operators
		inline constexpr bool operator!() const { return this->_x == 0; }
		inline constexpr bool operator==(const type b) const { return this->_x == b._x; }
		inline constexpr bool operator!=(const type b) const { return this->_x != b._x; }
		inline constexpr type operator ~() const { return type(~this->_x & type::_mask); }
		inline constexpr type operator & (const type b) const { return type(this->_x & b._x); }
		inline constexpr type operator | (const type b) const { return type(this->_x | b._x); }
		inline constexpr type operator ^ (const type b) const { return type(this->_x ^ b._x); }
		inline type& operator &=(const type b) { this->_x &= b._x; return *this; }
		inline type& operator |=(const type b) { this->_x |= b._x; return *this; }
		inline type& operator ^=(const type b) { this->_x ^= b._x; return *this; }

		// String/flag conversions
		inline const char* const rawname() const // returns "None", a predefined name, or NULL
		{
			if (this->_x == 0)           { return "None"; }
			//if (this->_x == type::_mask) { return "All";  }
			for (size_t i = 0; i < _count; ++i) { if (this->_x == type::_values[i]) { return type::_names[i]; } }
			return NULL;
		}
		inline const std::string name() const // returns "None", a predefined name or a string composed of predefined names separated by |
		{
			const char* const n = this->rawname();
			if (n) { return n; }
			std::string s = "";
			for (size_t i = 0; i < _count; ++i) { if ((this->_x & type::_values[i]) == type::_values[i]) { if (s.size()) { s += '|'; } s += type::_names[i]; } }
			return s;
		}
		friend inline std::ostream& operator<< (std::ostream& s, const type& b) { return s << b.name(); }
		static inline type from_name(const char *name, size_t n) // reverses "name()" accepting predefined names and "None" separated by | while ignoring whitespace before and after each name - throws std::invalid_argument if a name is not valid
		{
			while (*name && std::isspace(*name)) { ++name; --n; } // remove leading whitespace
			while (n > 0 && std::isspace(name[n-1])) { --n; } // remove training whitespace
			const char *pipe = (const char*)std::memchr(name, '|', n);
			if (pipe == NULL)
			{
				if (n == 4 && std::strncmp(name, "None", 4) == 0) { return type(0); }
				//if (n == 3 && std::strncmp(name, "All" , n)) { return type(type::_mask); }
				for (size_t i = 0; i < _count; ++i) { if (n == std::strlen(type::_names[i]) && std::strncmp(name, type::_names[i], n) == 0) { return type(type::_values[i]); } }
				throw std::invalid_argument("Flags name '"+std::string(name,n)+"' not found");
			}
			type b = type::from_name(name, pipe-name); n -= pipe - name; name = pipe + 1;
			while ((pipe = (const char*)std::memchr(name, '|', n)) != NULL)
			{
				b |= type::from_name(name, pipe-name); n -= pipe - name; name = pipe + 1;
			}
			b |= type::from_name(name, n-1);
			return b;
		}
		static inline type from_name(const char * name)      { return type::from_name(name,         std::strlen(name)); }
		static inline type from_name(const std::string name) { return type::from_name(name.c_str(), name.size());       }
	};
	
private:
	static inline constexpr base_type pow2(unsigned exp, base_type base = 2, base_type result = 1) { return exp < 1 ? result : pow2(exp/2, base*base, (exp % 2) ? result*base : result); }
	template <class dummy, unsigned... is>  struct indices                 { using next = indices<dummy, is..., sizeof...(is)>; using WithPow2Values = WithValues<pow2(is)...>; };
	template <class dummy, unsigned N>      struct build_indices           { using type = typename build_indices<dummy, N-1>::type::next; };
	template <class dummy>                  struct build_indices<dummy, 0> { using type = indices<dummy>; };

public:
	using WithDefaultValues = typename build_indices<void, _count>::type::WithPow2Values;
};

// make constexpr variables available at run-time
template <class BT_, template <class, class, BT_> class... Fs_> template <BT_... Vs_> constexpr const Flags<BT_, Fs_...>::WithValues<Vs_...> Flags<BT_, Fs_...>::WithValues<Vs_...>::None;
//template <class BT_, template <class, class, BT_> class... Fs_> template <BT_... Vs_> constexpr const Flags<BT_, Fs_...>::WithValues<Vs_...> Flags<BT_, Fs_...>::WithValues<Vs_...>::All;
template <class BT_, template <class,class,BT_> class... Fs_> template <BT_... Vs_> constexpr const char* const Flags<BT_,Fs_...>::WithValues<Vs_...>::_names [Flags<BT_,Fs_...>::_count];
template <class BT_, template <class,class,BT_> class... Fs_> template <BT_... Vs_> constexpr const BT_         Flags<BT_,Fs_...>::WithValues<Vs_...>::_values[Flags<BT_,Fs_...>::_count];

//template <class BT_, template <class,class,BT_> class... Fs_> template <BT_... Vs_> const std::vector<Flags<BT_,Fs_...>::WithValues<Vs_...>> Flags<BT_,Fs_...>::WithValues<Vs_...>::members = { Flags<BT_,Fs_...>::WithValues<Vs_...>(Vs_)... };


///////////////////////////////////////////////////////////////////////////////////////////////////
// Macro magic for making flags with default or defined values.
// Examples:
//   FLAGS(NameOfClass, FlagName1, FlagName2, FlagName3);
//   FLAGS_WITH_VALUES(NameOfClass, std::uint8_t, FlagName1, FlagName2, FlagName3, 1, 2, 4);
// When using flags created with default values, they are given the powers of two and the base
// type is automatically selected as the unsigned type that can minimally contain the number of
// flags. When providing your own values for the flags, you must specify a base type. Currently
// supports up to 32 flags in a single call but this could be expanded with copying and pasting.
//
// Note that this does "pollute" the containing namespace with a class for each flag name. One way
// around this is to do something like:
//   namespace detail { FLAGS(NameOfClass, FlagName1, FlagName2, FlagName3); }
//   typedef detail::NameOfClass NameOfClass;
// Also note that since FLAG_NAME cannot be used inside a class, these cannot be used inside a
// class either.
///////////////////////////////////////////////////////////////////////////////////////////////////

#define _xMAKE_FLAG_NAMES_1(F, ...)		FLAG_NAME(F)
#define _xMAKE_FLAG_NAMES_2(F, ...)		FLAG_NAME(F); _xMAKE_FLAG_NAMES_1(__VA_ARGS__)
#define _xMAKE_FLAG_NAMES_3(F, ...)		FLAG_NAME(F); _xMAKE_FLAG_NAMES_2(__VA_ARGS__)
#define _xMAKE_FLAG_NAMES_4(F, ...)		FLAG_NAME(F); _xMAKE_FLAG_NAMES_3(__VA_ARGS__)
#define _xMAKE_FLAG_NAMES_5(F, ...)		FLAG_NAME(F); _xMAKE_FLAG_NAMES_4(__VA_ARGS__)
#define _xMAKE_FLAG_NAMES_6(F, ...)		FLAG_NAME(F); _xMAKE_FLAG_NAMES_5(__VA_ARGS__)
#define _xMAKE_FLAG_NAMES_7(F, ...)		FLAG_NAME(F); _xMAKE_FLAG_NAMES_6(__VA_ARGS__)
#define _xMAKE_FLAG_NAMES_8(F, ...)		FLAG_NAME(F); _xMAKE_FLAG_NAMES_7(__VA_ARGS__)
#define _xMAKE_FLAG_NAMES_9(F, ...)		FLAG_NAME(F); _xMAKE_FLAG_NAMES_8(__VA_ARGS__)
#define _xMAKE_FLAG_NAMES_10(F, ...)	FLAG_NAME(F); _xMAKE_FLAG_NAMES_9(__VA_ARGS__)
#define _xMAKE_FLAG_NAMES_11(F, ...)	FLAG_NAME(F); _xMAKE_FLAG_NAMES_10(__VA_ARGS__)
#define _xMAKE_FLAG_NAMES_12(F, ...)	FLAG_NAME(F); _xMAKE_FLAG_NAMES_11(__VA_ARGS__)
#define _xMAKE_FLAG_NAMES_13(F, ...)	FLAG_NAME(F); _xMAKE_FLAG_NAMES_12(__VA_ARGS__)
#define _xMAKE_FLAG_NAMES_14(F, ...)	FLAG_NAME(F); _xMAKE_FLAG_NAMES_13(__VA_ARGS__)
#define _xMAKE_FLAG_NAMES_15(F, ...)	FLAG_NAME(F); _xMAKE_FLAG_NAMES_14(__VA_ARGS__)
#define _xMAKE_FLAG_NAMES_16(F, ...)	FLAG_NAME(F); _xMAKE_FLAG_NAMES_15(__VA_ARGS__)
#define _xMAKE_FLAG_NAMES_17(F, ...)	FLAG_NAME(F); _xMAKE_FLAG_NAMES_16(__VA_ARGS__)
#define _xMAKE_FLAG_NAMES_18(F, ...)	FLAG_NAME(F); _xMAKE_FLAG_NAMES_17(__VA_ARGS__)
#define _xMAKE_FLAG_NAMES_19(F, ...)	FLAG_NAME(F); _xMAKE_FLAG_NAMES_18(__VA_ARGS__)
#define _xMAKE_FLAG_NAMES_20(F, ...)	FLAG_NAME(F); _xMAKE_FLAG_NAMES_19(__VA_ARGS__)
#define _xMAKE_FLAG_NAMES_21(F, ...)	FLAG_NAME(F); _xMAKE_FLAG_NAMES_20(__VA_ARGS__)
#define _xMAKE_FLAG_NAMES_22(F, ...)	FLAG_NAME(F); _xMAKE_FLAG_NAMES_21(__VA_ARGS__)
#define _xMAKE_FLAG_NAMES_23(F, ...)	FLAG_NAME(F); _xMAKE_FLAG_NAMES_22(__VA_ARGS__)
#define _xMAKE_FLAG_NAMES_24(F, ...)	FLAG_NAME(F); _xMAKE_FLAG_NAMES_23(__VA_ARGS__)
#define _xMAKE_FLAG_NAMES_25(F, ...)	FLAG_NAME(F); _xMAKE_FLAG_NAMES_24(__VA_ARGS__)
#define _xMAKE_FLAG_NAMES_26(F, ...)	FLAG_NAME(F); _xMAKE_FLAG_NAMES_25(__VA_ARGS__)
#define _xMAKE_FLAG_NAMES_27(F, ...)	FLAG_NAME(F); _xMAKE_FLAG_NAMES_26(__VA_ARGS__)
#define _xMAKE_FLAG_NAMES_28(F, ...)	FLAG_NAME(F); _xMAKE_FLAG_NAMES_27(__VA_ARGS__)
#define _xMAKE_FLAG_NAMES_29(F, ...)	FLAG_NAME(F); _xMAKE_FLAG_NAMES_28(__VA_ARGS__)
#define _xMAKE_FLAG_NAMES_30(F, ...)	FLAG_NAME(F); _xMAKE_FLAG_NAMES_29(__VA_ARGS__)
#define _xMAKE_FLAG_NAMES_31(F, ...)	FLAG_NAME(F); _xMAKE_FLAG_NAMES_30(__VA_ARGS__)
#define _xMAKE_FLAG_NAMES_32(F, ...)	FLAG_NAME(F); _xMAKE_FLAG_NAMES_31(__VA_ARGS__)

#define _xFLAG_NAMES_1(F, ...)		Flag_##F
#define _xFLAG_NAMES_2(F, ...)		Flag_##F, _xFLAG_NAMES_1(__VA_ARGS__)
#define _xFLAG_NAMES_3(F, ...)		Flag_##F, _xFLAG_NAMES_2(__VA_ARGS__)
#define _xFLAG_NAMES_4(F, ...)		Flag_##F, _xFLAG_NAMES_3(__VA_ARGS__)
#define _xFLAG_NAMES_5(F, ...)		Flag_##F, _xFLAG_NAMES_4(__VA_ARGS__)
#define _xFLAG_NAMES_6(F, ...)		Flag_##F, _xFLAG_NAMES_5(__VA_ARGS__)
#define _xFLAG_NAMES_7(F, ...)		Flag_##F, _xFLAG_NAMES_6(__VA_ARGS__)
#define _xFLAG_NAMES_8(F, ...)		Flag_##F, _xFLAG_NAMES_7(__VA_ARGS__)
#define _xFLAG_NAMES_9(F, ...)		Flag_##F, _xFLAG_NAMES_8(__VA_ARGS__)
#define _xFLAG_NAMES_10(F, ...)		Flag_##F, _xFLAG_NAMES_9(__VA_ARGS__)
#define _xFLAG_NAMES_11(F, ...)		Flag_##F, _xFLAG_NAMES_10(__VA_ARGS__)
#define _xFLAG_NAMES_12(F, ...)		Flag_##F, _xFLAG_NAMES_11(__VA_ARGS__)
#define _xFLAG_NAMES_13(F, ...)		Flag_##F, _xFLAG_NAMES_12(__VA_ARGS__)
#define _xFLAG_NAMES_14(F, ...)		Flag_##F, _xFLAG_NAMES_13(__VA_ARGS__)
#define _xFLAG_NAMES_15(F, ...)		Flag_##F, _xFLAG_NAMES_14(__VA_ARGS__)
#define _xFLAG_NAMES_16(F, ...)		Flag_##F, _xFLAG_NAMES_15(__VA_ARGS__)
#define _xFLAG_NAMES_17(F, ...)		Flag_##F, _xFLAG_NAMES_16(__VA_ARGS__)
#define _xFLAG_NAMES_18(F, ...)		Flag_##F, _xFLAG_NAMES_17(__VA_ARGS__)
#define _xFLAG_NAMES_19(F, ...)		Flag_##F, _xFLAG_NAMES_18(__VA_ARGS__)
#define _xFLAG_NAMES_20(F, ...)		Flag_##F, _xFLAG_NAMES_19(__VA_ARGS__)
#define _xFLAG_NAMES_21(F, ...)		Flag_##F, _xFLAG_NAMES_20(__VA_ARGS__)
#define _xFLAG_NAMES_22(F, ...)		Flag_##F, _xFLAG_NAMES_21(__VA_ARGS__)
#define _xFLAG_NAMES_23(F, ...)		Flag_##F, _xFLAG_NAMES_22(__VA_ARGS__)
#define _xFLAG_NAMES_24(F, ...)		Flag_##F, _xFLAG_NAMES_23(__VA_ARGS__)
#define _xFLAG_NAMES_25(F, ...)		Flag_##F, _xFLAG_NAMES_24(__VA_ARGS__)
#define _xFLAG_NAMES_26(F, ...)		Flag_##F, _xFLAG_NAMES_25(__VA_ARGS__)
#define _xFLAG_NAMES_27(F, ...)		Flag_##F, _xFLAG_NAMES_26(__VA_ARGS__)
#define _xFLAG_NAMES_28(F, ...)		Flag_##F, _xFLAG_NAMES_27(__VA_ARGS__)
#define _xFLAG_NAMES_29(F, ...)		Flag_##F, _xFLAG_NAMES_28(__VA_ARGS__)
#define _xFLAG_NAMES_30(F, ...)		Flag_##F, _xFLAG_NAMES_29(__VA_ARGS__)
#define _xFLAG_NAMES_31(F, ...)		Flag_##F, _xFLAG_NAMES_30(__VA_ARGS__)
#define _xFLAG_NAMES_32(F, ...)		Flag_##F, _xFLAG_NAMES_31(__VA_ARGS__)

#define _xFLAG_VALUES_1(V)			V
#define _xFLAG_VALUES_2(V, ...)		V, _xFLAG_VALUES_1(__VA_ARGS__)
#define _xFLAG_VALUES_3(V, ...)		V, _xFLAG_VALUES_2(__VA_ARGS__)
#define _xFLAG_VALUES_4(V, ...)		V, _xFLAG_VALUES_3(__VA_ARGS__)
#define _xFLAG_VALUES_5(V, ...)		V, _xFLAG_VALUES_4(__VA_ARGS__)
#define _xFLAG_VALUES_6(V, ...)		V, _xFLAG_VALUES_5(__VA_ARGS__)
#define _xFLAG_VALUES_7(V, ...)		V, _xFLAG_VALUES_6(__VA_ARGS__)
#define _xFLAG_VALUES_8(V, ...)		V, _xFLAG_VALUES_7(__VA_ARGS__)
#define _xFLAG_VALUES_9(V, ...)		V, _xFLAG_VALUES_8(__VA_ARGS__)
#define _xFLAG_VALUES_10(V, ...)	V, _xFLAG_VALUES_9(__VA_ARGS__)
#define _xFLAG_VALUES_11(V, ...)	V, _xFLAG_VALUES_10(__VA_ARGS__)
#define _xFLAG_VALUES_12(V, ...)	V, _xFLAG_VALUES_11(__VA_ARGS__)
#define _xFLAG_VALUES_13(V, ...)	V, _xFLAG_VALUES_12(__VA_ARGS__)
#define _xFLAG_VALUES_14(V, ...)	V, _xFLAG_VALUES_13(__VA_ARGS__)
#define _xFLAG_VALUES_15(V, ...)	V, _xFLAG_VALUES_14(__VA_ARGS__)
#define _xFLAG_VALUES_16(V, ...)	V, _xFLAG_VALUES_15(__VA_ARGS__)
#define _xFLAG_VALUES_17(V, ...)	V, _xFLAG_VALUES_16(__VA_ARGS__)
#define _xFLAG_VALUES_18(V, ...)	V, _xFLAG_VALUES_17(__VA_ARGS__)
#define _xFLAG_VALUES_19(V, ...)	V, _xFLAG_VALUES_18(__VA_ARGS__)
#define _xFLAG_VALUES_20(V, ...)	V, _xFLAG_VALUES_19(__VA_ARGS__)
#define _xFLAG_VALUES_21(V, ...)	V, _xFLAG_VALUES_20(__VA_ARGS__)
#define _xFLAG_VALUES_22(V, ...)	V, _xFLAG_VALUES_21(__VA_ARGS__)
#define _xFLAG_VALUES_23(V, ...)	V, _xFLAG_VALUES_22(__VA_ARGS__)
#define _xFLAG_VALUES_24(V, ...)	V, _xFLAG_VALUES_23(__VA_ARGS__)
#define _xFLAG_VALUES_25(V, ...)	V, _xFLAG_VALUES_24(__VA_ARGS__)
#define _xFLAG_VALUES_26(V, ...)	V, _xFLAG_VALUES_25(__VA_ARGS__)
#define _xFLAG_VALUES_27(V, ...)	V, _xFLAG_VALUES_26(__VA_ARGS__)
#define _xFLAG_VALUES_28(V, ...)	V, _xFLAG_VALUES_27(__VA_ARGS__)
#define _xFLAG_VALUES_29(V, ...)	V, _xFLAG_VALUES_28(__VA_ARGS__)
#define _xFLAG_VALUES_30(V, ...)	V, _xFLAG_VALUES_29(__VA_ARGS__)
#define _xFLAG_VALUES_31(V, ...)	V, _xFLAG_VALUES_30(__VA_ARGS__)
#define _xFLAG_VALUES_32(V, ...)	V, _xFLAG_VALUES_31(__VA_ARGS__)

#define _xSKIP_FLAG_NAMES_1(N, F, ...)	_xFLAG_VALUES_##N(__VA_ARGS__)
#define _xSKIP_FLAG_NAMES_2(N, F, ...)	_xSKIP_FLAG_NAMES_1(N, __VA_ARGS__)
#define _xSKIP_FLAG_NAMES_3(N, F, ...)	_xSKIP_FLAG_NAMES_2(N, __VA_ARGS__)
#define _xSKIP_FLAG_NAMES_4(N, F, ...)	_xSKIP_FLAG_NAMES_3(N, __VA_ARGS__)
#define _xSKIP_FLAG_NAMES_5(N, F, ...)	_xSKIP_FLAG_NAMES_4(N, __VA_ARGS__)
#define _xSKIP_FLAG_NAMES_6(N, F, ...)	_xSKIP_FLAG_NAMES_5(N, __VA_ARGS__)
#define _xSKIP_FLAG_NAMES_7(N, F, ...)	_xSKIP_FLAG_NAMES_6(N, __VA_ARGS__)
#define _xSKIP_FLAG_NAMES_8(N, F, ...)	_xSKIP_FLAG_NAMES_7(N, __VA_ARGS__)
#define _xSKIP_FLAG_NAMES_9(N, F, ...)	_xSKIP_FLAG_NAMES_8(N, __VA_ARGS__)
#define _xSKIP_FLAG_NAMES_10(N, F, ...)	_xSKIP_FLAG_NAMES_9(N, __VA_ARGS__)
#define _xSKIP_FLAG_NAMES_11(N, F, ...)	_xSKIP_FLAG_NAMES_10(N, __VA_ARGS__)
#define _xSKIP_FLAG_NAMES_12(N, F, ...)	_xSKIP_FLAG_NAMES_11(N, __VA_ARGS__)
#define _xSKIP_FLAG_NAMES_13(N, F, ...)	_xSKIP_FLAG_NAMES_12(N, __VA_ARGS__)
#define _xSKIP_FLAG_NAMES_14(N, F, ...)	_xSKIP_FLAG_NAMES_13(N, __VA_ARGS__)
#define _xSKIP_FLAG_NAMES_15(N, F, ...)	_xSKIP_FLAG_NAMES_14(N, __VA_ARGS__)
#define _xSKIP_FLAG_NAMES_16(N, F, ...)	_xSKIP_FLAG_NAMES_15(N, __VA_ARGS__)
#define _xSKIP_FLAG_NAMES_17(N, F, ...)	_xSKIP_FLAG_NAMES_16(N, __VA_ARGS__)
#define _xSKIP_FLAG_NAMES_18(N, F, ...)	_xSKIP_FLAG_NAMES_17(N, __VA_ARGS__)
#define _xSKIP_FLAG_NAMES_19(N, F, ...)	_xSKIP_FLAG_NAMES_18(N, __VA_ARGS__)
#define _xSKIP_FLAG_NAMES_20(N, F, ...)	_xSKIP_FLAG_NAMES_19(N, __VA_ARGS__)
#define _xSKIP_FLAG_NAMES_21(N, F, ...)	_xSKIP_FLAG_NAMES_20(N, __VA_ARGS__)
#define _xSKIP_FLAG_NAMES_22(N, F, ...)	_xSKIP_FLAG_NAMES_21(N, __VA_ARGS__)
#define _xSKIP_FLAG_NAMES_23(N, F, ...)	_xSKIP_FLAG_NAMES_22(N, __VA_ARGS__)
#define _xSKIP_FLAG_NAMES_24(N, F, ...)	_xSKIP_FLAG_NAMES_23(N, __VA_ARGS__)
#define _xSKIP_FLAG_NAMES_25(N, F, ...)	_xSKIP_FLAG_NAMES_24(N, __VA_ARGS__)
#define _xSKIP_FLAG_NAMES_26(N, F, ...)	_xSKIP_FLAG_NAMES_25(N, __VA_ARGS__)
#define _xSKIP_FLAG_NAMES_27(N, F, ...)	_xSKIP_FLAG_NAMES_26(N, __VA_ARGS__)
#define _xSKIP_FLAG_NAMES_28(N, F, ...)	_xSKIP_FLAG_NAMES_27(N, __VA_ARGS__)
#define _xSKIP_FLAG_NAMES_29(N, F, ...)	_xSKIP_FLAG_NAMES_28(N, __VA_ARGS__)
#define _xSKIP_FLAG_NAMES_30(N, F, ...)	_xSKIP_FLAG_NAMES_29(N, __VA_ARGS__)
#define _xSKIP_FLAG_NAMES_31(N, F, ...)	_xSKIP_FLAG_NAMES_30(N, __VA_ARGS__)
#define _xSKIP_FLAG_NAMES_32(N, F, ...)	_xSKIP_FLAG_NAMES_31(N, __VA_ARGS__)

#define _xFLAGS_1(NAME, F)		_xMAKE_FLAG_NAMES_1(F);            typedef ::Flags<std::uint8_t,_xFLAG_NAMES_1(F)          >::WithDefaultValues NAME
#define _xFLAGS_2(NAME, ...)	_xMAKE_FLAG_NAMES_2(__VA_ARGS__);  typedef ::Flags<std::uint8_t,_xFLAG_NAMES_2(__VA_ARGS__)>::WithDefaultValues NAME
#define _xFLAGS_3(NAME, ...)	_xMAKE_FLAG_NAMES_3(__VA_ARGS__);  typedef ::Flags<std::uint8_t,_xFLAG_NAMES_3(__VA_ARGS__)>::WithDefaultValues NAME
#define _xFLAGS_4(NAME, ...)	_xMAKE_FLAG_NAMES_4(__VA_ARGS__);  typedef ::Flags<std::uint8_t,_xFLAG_NAMES_4(__VA_ARGS__)>::WithDefaultValues NAME
#define _xFLAGS_5(NAME, ...)	_xMAKE_FLAG_NAMES_5(__VA_ARGS__);  typedef ::Flags<std::uint8_t,_xFLAG_NAMES_5(__VA_ARGS__)>::WithDefaultValues NAME
#define _xFLAGS_6(NAME, ...)	_xMAKE_FLAG_NAMES_6(__VA_ARGS__);  typedef ::Flags<std::uint8_t,_xFLAG_NAMES_6(__VA_ARGS__)>::WithDefaultValues NAME
#define _xFLAGS_7(NAME, ...)	_xMAKE_FLAG_NAMES_7(__VA_ARGS__);  typedef ::Flags<std::uint8_t,_xFLAG_NAMES_7(__VA_ARGS__)>::WithDefaultValues NAME
#define _xFLAGS_8(NAME, ...)	_xMAKE_FLAG_NAMES_8(__VA_ARGS__);  typedef ::Flags<std::uint8_t,_xFLAG_NAMES_8(__VA_ARGS__)>::WithDefaultValues NAME
#define _xFLAGS_9(NAME, ...)	_xMAKE_FLAG_NAMES_9(__VA_ARGS__);  typedef ::Flags<std::uint16_t,_xFLAG_NAMES_9(__VA_ARGS__)>::WithDefaultValues NAME
#define _xFLAGS_10(NAME, ...)	_xMAKE_FLAG_NAMES_10(__VA_ARGS__); typedef ::Flags<std::uint16_t,_xFLAG_NAMES_10(__VA_ARGS__)>::WithDefaultValues NAME
#define _xFLAGS_11(NAME, ...)	_xMAKE_FLAG_NAMES_11(__VA_ARGS__); typedef ::Flags<std::uint16_t,_xFLAG_NAMES_11(__VA_ARGS__)>::WithDefaultValues NAME
#define _xFLAGS_12(NAME, ...)	_xMAKE_FLAG_NAMES_12(__VA_ARGS__); typedef ::Flags<std::uint16_t,_xFLAG_NAMES_12(__VA_ARGS__)>::WithDefaultValues NAME
#define _xFLAGS_13(NAME, ...)	_xMAKE_FLAG_NAMES_13(__VA_ARGS__); typedef ::Flags<std::uint16_t,_xFLAG_NAMES_13(__VA_ARGS__)>::WithDefaultValues NAME
#define _xFLAGS_14(NAME, ...)	_xMAKE_FLAG_NAMES_14(__VA_ARGS__); typedef ::Flags<std::uint16_t,_xFLAG_NAMES_14(__VA_ARGS__)>::WithDefaultValues NAME
#define _xFLAGS_15(NAME, ...)	_xMAKE_FLAG_NAMES_15(__VA_ARGS__); typedef ::Flags<std::uint16_t,_xFLAG_NAMES_15(__VA_ARGS__)>::WithDefaultValues NAME
#define _xFLAGS_16(NAME, ...)	_xMAKE_FLAG_NAMES_16(__VA_ARGS__); typedef ::Flags<std::uint16_t,_xFLAG_NAMES_16(__VA_ARGS__)>::WithDefaultValues NAME
#define _xFLAGS_17(NAME, ...)	_xMAKE_FLAG_NAMES_17(__VA_ARGS__); typedef ::Flags<std::uint32_t,_xFLAG_NAMES_17(__VA_ARGS__)>::WithDefaultValues NAME
#define _xFLAGS_18(NAME, ...)	_xMAKE_FLAG_NAMES_18(__VA_ARGS__); typedef ::Flags<std::uint32_t,_xFLAG_NAMES_18(__VA_ARGS__)>::WithDefaultValues NAME
#define _xFLAGS_19(NAME, ...)	_xMAKE_FLAG_NAMES_19(__VA_ARGS__); typedef ::Flags<std::uint32_t,_xFLAG_NAMES_19(__VA_ARGS__)>::WithDefaultValues NAME
#define _xFLAGS_20(NAME, ...)	_xMAKE_FLAG_NAMES_20(__VA_ARGS__); typedef ::Flags<std::uint32_t,_xFLAG_NAMES_20(__VA_ARGS__)>::WithDefaultValues NAME
#define _xFLAGS_21(NAME, ...)	_xMAKE_FLAG_NAMES_21(__VA_ARGS__); typedef ::Flags<std::uint32_t,_xFLAG_NAMES_21(__VA_ARGS__)>::WithDefaultValues NAME
#define _xFLAGS_22(NAME, ...)	_xMAKE_FLAG_NAMES_22(__VA_ARGS__); typedef ::Flags<std::uint32_t,_xFLAG_NAMES_22(__VA_ARGS__)>::WithDefaultValues NAME
#define _xFLAGS_23(NAME, ...)	_xMAKE_FLAG_NAMES_23(__VA_ARGS__); typedef ::Flags<std::uint32_t,_xFLAG_NAMES_23(__VA_ARGS__)>::WithDefaultValues NAME
#define _xFLAGS_24(NAME, ...)	_xMAKE_FLAG_NAMES_24(__VA_ARGS__); typedef ::Flags<std::uint32_t,_xFLAG_NAMES_24(__VA_ARGS__)>::WithDefaultValues NAME
#define _xFLAGS_25(NAME, ...)	_xMAKE_FLAG_NAMES_25(__VA_ARGS__); typedef ::Flags<std::uint32_t,_xFLAG_NAMES_25(__VA_ARGS__)>::WithDefaultValues NAME
#define _xFLAGS_26(NAME, ...)	_xMAKE_FLAG_NAMES_26(__VA_ARGS__); typedef ::Flags<std::uint32_t,_xFLAG_NAMES_26(__VA_ARGS__)>::WithDefaultValues NAME
#define _xFLAGS_27(NAME, ...)	_xMAKE_FLAG_NAMES_27(__VA_ARGS__); typedef ::Flags<std::uint32_t,_xFLAG_NAMES_27(__VA_ARGS__)>::WithDefaultValues NAME
#define _xFLAGS_28(NAME, ...)	_xMAKE_FLAG_NAMES_28(__VA_ARGS__); typedef ::Flags<std::uint32_t,_xFLAG_NAMES_28(__VA_ARGS__)>::WithDefaultValues NAME
#define _xFLAGS_29(NAME, ...)	_xMAKE_FLAG_NAMES_29(__VA_ARGS__); typedef ::Flags<std::uint32_t,_xFLAG_NAMES_29(__VA_ARGS__)>::WithDefaultValues NAME
#define _xFLAGS_30(NAME, ...)	_xMAKE_FLAG_NAMES_30(__VA_ARGS__); typedef ::Flags<std::uint32_t,_xFLAG_NAMES_30(__VA_ARGS__)>::WithDefaultValues NAME
#define _xFLAGS_31(NAME, ...)	_xMAKE_FLAG_NAMES_31(__VA_ARGS__); typedef ::Flags<std::uint32_t,_xFLAG_NAMES_31(__VA_ARGS__)>::WithDefaultValues NAME
#define _xFLAGS_32(NAME, ...)	_xMAKE_FLAG_NAMES_32(__VA_ARGS__); typedef ::Flags<std::uint32_t,_xFLAG_NAMES_32(__VA_ARGS__)>::WithDefaultValues NAME

#define _xFLAGS_WITH_VALUES_2(NAME, BT, F, V)	_xMAKE_FLAG_NAMES_1(F);            typedef ::Flags<BT,_xFLAG_NAMES_1(F)          >::WithValues<_xSKIP_FLAG_NAMES_1(1,F,V)        > NAME
#define _xFLAGS_WITH_VALUES_4(NAME, BT, ...)	_xMAKE_FLAG_NAMES_2(__VA_ARGS__);  typedef ::Flags<BT,_xFLAG_NAMES_2(__VA_ARGS__)>::WithValues<_xSKIP_FLAG_NAMES_2(2,__VA_ARGS__)> NAME
#define _xFLAGS_WITH_VALUES_6(NAME, BT, ...)	_xMAKE_FLAG_NAMES_3(__VA_ARGS__);  typedef ::Flags<BT,_xFLAG_NAMES_3(__VA_ARGS__)>::WithValues<_xSKIP_FLAG_NAMES_3(3,__VA_ARGS__)> NAME
#define _xFLAGS_WITH_VALUES_8(NAME, BT, ...)	_xMAKE_FLAG_NAMES_4(__VA_ARGS__);  typedef ::Flags<BT,_xFLAG_NAMES_4(__VA_ARGS__)>::WithValues<_xSKIP_FLAG_NAMES_4(4,__VA_ARGS__)> NAME
#define _xFLAGS_WITH_VALUES_10(NAME, BT, ...)	_xMAKE_FLAG_NAMES_5(__VA_ARGS__);  typedef ::Flags<BT,_xFLAG_NAMES_5(__VA_ARGS__)>::WithValues<_xSKIP_FLAG_NAMES_5(5,__VA_ARGS__)> NAME
#define _xFLAGS_WITH_VALUES_12(NAME, BT, ...)	_xMAKE_FLAG_NAMES_6(__VA_ARGS__);  typedef ::Flags<BT,_xFLAG_NAMES_6(__VA_ARGS__)>::WithValues<_xSKIP_FLAG_NAMES_6(6,__VA_ARGS__)> NAME
#define _xFLAGS_WITH_VALUES_14(NAME, BT, ...)	_xMAKE_FLAG_NAMES_7(__VA_ARGS__);  typedef ::Flags<BT,_xFLAG_NAMES_7(__VA_ARGS__)>::WithValues<_xSKIP_FLAG_NAMES_7(7,__VA_ARGS__)> NAME
#define _xFLAGS_WITH_VALUES_16(NAME, BT, ...)	_xMAKE_FLAG_NAMES_8(__VA_ARGS__);  typedef ::Flags<BT,_xFLAG_NAMES_8(__VA_ARGS__)>::WithValues<_xSKIP_FLAG_NAMES_8(8,__VA_ARGS__)> NAME
#define _xFLAGS_WITH_VALUES_18(NAME, BT, ...)	_xMAKE_FLAG_NAMES_9(__VA_ARGS__);  typedef ::Flags<BT,_xFLAG_NAMES_9(__VA_ARGS__)>::WithValues<_xSKIP_FLAG_NAMES_9(9,__VA_ARGS__)> NAME
#define _xFLAGS_WITH_VALUES_20(NAME, BT, ...)	_xMAKE_FLAG_NAMES_10(__VA_ARGS__); typedef ::Flags<BT,_xFLAG_NAMES_10(__VA_ARGS__)>::WithValues<_xSKIP_FLAG_NAMES_10(10,__VA_ARGS__)> NAME
#define _xFLAGS_WITH_VALUES_22(NAME, BT, ...)	_xMAKE_FLAG_NAMES_11(__VA_ARGS__); typedef ::Flags<BT,_xFLAG_NAMES_11(__VA_ARGS__)>::WithValues<_xSKIP_FLAG_NAMES_11(11,__VA_ARGS__)> NAME
#define _xFLAGS_WITH_VALUES_24(NAME, BT, ...)	_xMAKE_FLAG_NAMES_12(__VA_ARGS__); typedef ::Flags<BT,_xFLAG_NAMES_12(__VA_ARGS__)>::WithValues<_xSKIP_FLAG_NAMES_12(12,__VA_ARGS__)> NAME
#define _xFLAGS_WITH_VALUES_26(NAME, BT, ...)	_xMAKE_FLAG_NAMES_13(__VA_ARGS__); typedef ::Flags<BT,_xFLAG_NAMES_13(__VA_ARGS__)>::WithValues<_xSKIP_FLAG_NAMES_13(13,__VA_ARGS__)> NAME
#define _xFLAGS_WITH_VALUES_28(NAME, BT, ...)	_xMAKE_FLAG_NAMES_14(__VA_ARGS__); typedef ::Flags<BT,_xFLAG_NAMES_14(__VA_ARGS__)>::WithValues<_xSKIP_FLAG_NAMES_14(14,__VA_ARGS__)> NAME
#define _xFLAGS_WITH_VALUES_30(NAME, BT, ...)	_xMAKE_FLAG_NAMES_15(__VA_ARGS__); typedef ::Flags<BT,_xFLAG_NAMES_15(__VA_ARGS__)>::WithValues<_xSKIP_FLAG_NAMES_15(15,__VA_ARGS__)> NAME
#define _xFLAGS_WITH_VALUES_32(NAME, BT, ...)	_xMAKE_FLAG_NAMES_16(__VA_ARGS__); typedef ::Flags<BT,_xFLAG_NAMES_16(__VA_ARGS__)>::WithValues<_xSKIP_FLAG_NAMES_16(16,__VA_ARGS__)> NAME
#define _xFLAGS_WITH_VALUES_34(NAME, BT, ...)	_xMAKE_FLAG_NAMES_17(__VA_ARGS__); typedef ::Flags<BT,_xFLAG_NAMES_17(__VA_ARGS__)>::WithValues<_xSKIP_FLAG_NAMES_17(17,__VA_ARGS__)> NAME
#define _xFLAGS_WITH_VALUES_36(NAME, BT, ...)	_xMAKE_FLAG_NAMES_18(__VA_ARGS__); typedef ::Flags<BT,_xFLAG_NAMES_18(__VA_ARGS__)>::WithValues<_xSKIP_FLAG_NAMES_18(18,__VA_ARGS__)> NAME
#define _xFLAGS_WITH_VALUES_38(NAME, BT, ...)	_xMAKE_FLAG_NAMES_19(__VA_ARGS__); typedef ::Flags<BT,_xFLAG_NAMES_19(__VA_ARGS__)>::WithValues<_xSKIP_FLAG_NAMES_19(19,__VA_ARGS__)> NAME
#define _xFLAGS_WITH_VALUES_40(NAME, BT, ...)	_xMAKE_FLAG_NAMES_20(__VA_ARGS__); typedef ::Flags<BT,_xFLAG_NAMES_20(__VA_ARGS__)>::WithValues<_xSKIP_FLAG_NAMES_20(20,__VA_ARGS__)> NAME
#define _xFLAGS_WITH_VALUES_42(NAME, BT, ...)	_xMAKE_FLAG_NAMES_21(__VA_ARGS__); typedef ::Flags<BT,_xFLAG_NAMES_21(__VA_ARGS__)>::WithValues<_xSKIP_FLAG_NAMES_21(21,__VA_ARGS__)> NAME
#define _xFLAGS_WITH_VALUES_44(NAME, BT, ...)	_xMAKE_FLAG_NAMES_22(__VA_ARGS__); typedef ::Flags<BT,_xFLAG_NAMES_22(__VA_ARGS__)>::WithValues<_xSKIP_FLAG_NAMES_22(22,__VA_ARGS__)> NAME
#define _xFLAGS_WITH_VALUES_46(NAME, BT, ...)	_xMAKE_FLAG_NAMES_23(__VA_ARGS__); typedef ::Flags<BT,_xFLAG_NAMES_23(__VA_ARGS__)>::WithValues<_xSKIP_FLAG_NAMES_23(23,__VA_ARGS__)> NAME
#define _xFLAGS_WITH_VALUES_48(NAME, BT, ...)	_xMAKE_FLAG_NAMES_24(__VA_ARGS__); typedef ::Flags<BT,_xFLAG_NAMES_24(__VA_ARGS__)>::WithValues<_xSKIP_FLAG_NAMES_24(24,__VA_ARGS__)> NAME
#define _xFLAGS_WITH_VALUES_50(NAME, BT, ...)	_xMAKE_FLAG_NAMES_25(__VA_ARGS__); typedef ::Flags<BT,_xFLAG_NAMES_25(__VA_ARGS__)>::WithValues<_xSKIP_FLAG_NAMES_25(25,__VA_ARGS__)> NAME
#define _xFLAGS_WITH_VALUES_52(NAME, BT, ...)	_xMAKE_FLAG_NAMES_26(__VA_ARGS__); typedef ::Flags<BT,_xFLAG_NAMES_26(__VA_ARGS__)>::WithValues<_xSKIP_FLAG_NAMES_26(26,__VA_ARGS__)> NAME
#define _xFLAGS_WITH_VALUES_54(NAME, BT, ...)	_xMAKE_FLAG_NAMES_27(__VA_ARGS__); typedef ::Flags<BT,_xFLAG_NAMES_27(__VA_ARGS__)>::WithValues<_xSKIP_FLAG_NAMES_27(27,__VA_ARGS__)> NAME
#define _xFLAGS_WITH_VALUES_56(NAME, BT, ...)	_xMAKE_FLAG_NAMES_28(__VA_ARGS__); typedef ::Flags<BT,_xFLAG_NAMES_28(__VA_ARGS__)>::WithValues<_xSKIP_FLAG_NAMES_28(28,__VA_ARGS__)> NAME
#define _xFLAGS_WITH_VALUES_58(NAME, BT, ...)	_xMAKE_FLAG_NAMES_29(__VA_ARGS__); typedef ::Flags<BT,_xFLAG_NAMES_29(__VA_ARGS__)>::WithValues<_xSKIP_FLAG_NAMES_29(29,__VA_ARGS__)> NAME
#define _xFLAGS_WITH_VALUES_60(NAME, BT, ...)	_xMAKE_FLAG_NAMES_30(__VA_ARGS__); typedef ::Flags<BT,_xFLAG_NAMES_30(__VA_ARGS__)>::WithValues<_xSKIP_FLAG_NAMES_30(30,__VA_ARGS__)> NAME
#define _xFLAGS_WITH_VALUES_62(NAME, BT, ...)	_xMAKE_FLAG_NAMES_31(__VA_ARGS__); typedef ::Flags<BT,_xFLAG_NAMES_31(__VA_ARGS__)>::WithValues<_xSKIP_FLAG_NAMES_31(31,__VA_ARGS__)> NAME
#define _xFLAGS_WITH_VALUES_64(NAME, BT, ...)	_xMAKE_FLAG_NAMES_32(__VA_ARGS__); typedef ::Flags<BT,_xFLAG_NAMES_32(__VA_ARGS__)>::WithValues<_xSKIP_FLAG_NAMES_32(32,__VA_ARGS__)> NAME

// NUM_ARGS(...) evaluates to the literal number of the passed-in arguments.
#define _xNUM_ARGS2(X,X64,X63,X62,X61,X60,X59,X58,X57,X56,X55,X54,X53,X52,X51,X50,X49,X48,X47,X46,X45,X44,X43,X42,X41,X40,X39,X38,X37,X36,X35,X34,X33,X32,X31,X30,X29,X28,X27,X26,X25,X24,X23,X22,X21,X20,X19,X18,X17,X16,X15,X14,X13,X12,X11,X10,X9,X8,X7,X6,X5,X4,X3,X2,X1,N,...) N
#define NUM_ARGS(...) _xNUM_ARGS2(0,__VA_ARGS__ ,64,63,62,61,60,59,58,57,56,55,54,53,52,51,50,49,48,47,46,45,44,43,42,41,40,39,38,37,36,35,34,33,32,31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0)

#define _xFLAGS3(NAME, N, ...)	_xFLAGS_##N(NAME, __VA_ARGS__)
#define _xFLAGS2(NAME, N, ...)	_xFLAGS3(NAME, N, __VA_ARGS__)
#define FLAGS(NAME, ...)		_xFLAGS2(NAME, NUM_ARGS(__VA_ARGS__), __VA_ARGS__)

#define _xFLAGS_WITH_VALUES3(NAME, BT, N, ...)	_xFLAGS_WITH_VALUES_##N(NAME, BT, __VA_ARGS__)
#define _xFLAGS_WITH_VALUES2(NAME, BT, N, ...)	_xFLAGS_WITH_VALUES3(NAME, BT, N, __VA_ARGS__)
#define FLAGS_WITH_VALUES(NAME, BT, ...)		_xFLAGS_WITH_VALUES2(NAME, BT, NUM_ARGS(__VA_ARGS__), __VA_ARGS__)
#endif // _MSC_VER
