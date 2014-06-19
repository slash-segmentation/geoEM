#pragma once

// This is a special version for MSVC since MSVC only half-heartily supports C++11 in its newest
// version. This avoids using constexpr, variadic templates, and has a backup when explicit casts
// are not supported. Note that even though VS2013 "supports" variadic templates, there are many
// bugs with template variadic templates their usage in the Flags class is completely not
// supported. Maybe MSVC 2014 will provide full support for constexpr and variadic templates thus
// enabling us to make it so this hacked version is not used for it, until then, we need this. When
// it does happen we may need to support MSVCs interpretation of variadic macros (requiring EXPAND)
// into the default flags.

#include <type_traits>
#include <stdexcept>
#include <iostream>
#include <utility>
#include <cstdint>
#include <cstring>
#include <cctype>
#include <string>
#include <array>

#if _MSC_VER < 1600
#error Only Visual Studio 2010 and newer are supported.
#endif
#if !defined(BOOST_NO_CXX11_EXPLICIT_CONVERSION_OPERATORS) && _MSC_FULL_VER < 180020827
#define BOOST_NO_CXX11_EXPLICIT_CONVERSION_OPERATORS
#endif

// Creates a name to use as a flag
#define FLAG_NAME(NAME) FLAG_NAME_CORE(NAME) FLAG_NAME_INIT(NAME)
#define FLAG_NAME_CORE(NAME) \
	template <class BT_, class F_, BT_ V_> struct Flag_##NAME \
	{ \
		static_assert(V_ != 0, "No flag may have a value of 0 (use the special None flag for that)"); \
		static const F_ NAME; \
		protected: static const char* const _name() { return #NAME; } \
	};
#define FLAG_NAME_INIT(NAME) \
	template<class BT_, class F_, BT_ V_> const F_ Flag_##NAME<BT_,F_,V_>::NAME(V_) /* make available at runtime (only needed when flag names are used inline) */
#define FLAG_NAME_INIT_NS(NS, NAME) \
	template<class BT_, class F_, BT_ V_> const F_ NS::Flag_##NAME<BT_,F_,V_>::NAME(V_) /* make available at runtime (only needed when flag names are used inline) */


namespace detail
{
	template <class BT_, class FT_, size_t _count_>
	class _Flags_base
	{
	public:
		typedef BT_ base_type;
		typedef FT_ type;
		static_assert(std::is_integral<base_type>::value, "Flags base type must be an integral type");

		// Special flags
		static const type None; // no flags set
		//static const type All; // all flags set (could do ~None as well)

	protected:
		static const size_t _count = _count_;

		// Flag value (only non-static variable)
		base_type _x;

		static const char* const* _names()
		{
			static const char* _names[_count] = { NULL };
			if (_names[0] == NULL) { type::_fill_names(_names); }
			return _names;
		}
		static const base_type* _values()
		{
			static base_type _values[_count] = { 0 };
			if (_values[0] == 0) { type::_fill_values(_values); }
			return _values;
		}

	public:
		static const std::array<type, _count>& members()
		{
			static std::array<type, _count> m;
			if (m[0] == type::None) { const base_type* vs = type::_values(); for (size_t i = 0; i < _count; ++i) { m[i] = type(vs[i]); } }
			return m;
		}

		// Constructors, assignment operators, and casts
		inline _Flags_base() : _x(0) { }
		inline _Flags_base(const type& b) : _x(b._x) { }
		//inline type& operator=(const type& b) { this->_x = b._x; return *this; }
		inline explicit _Flags_base(base_type b) : _x(b) { if ((b | type::_mask) != type::_mask) { throw std::domain_error("Flag value outside mask"); } }
		//inline type& operator=(base_type b) { if ((b|type::_mask)!=type::_mask) { throw std::domain_error("Flag value outside mask"); } this->_x = b; return *this; }
		inline base_type operator+() const { return this->_x; }

#ifdef BOOST_NO_CXX11_EXPLICIT_CONVERSION_OPERATORS
	private:
		struct _dummy { void nonnull() {}; };
		typedef void (_dummy::*_safe_bool)();
	public:
		inline operator _safe_bool() const { return this->_x != 0 ? &_dummy::nonnull : 0; }
		// no cast to base_type since it can't be explicit, use unary + operator instead
#else
		inline explicit operator bool()      const { return this->_x != 0; }
		inline explicit operator base_type() const { return this->_x; }
#endif

		// Bit-wise operators
		inline bool operator!() const { return this->_x == 0; }
		inline bool operator==(const type b) const { return this->_x == b._x; }
		inline bool operator!=(const type b) const { return this->_x != b._x; }
		inline type operator ~() const { return type(~this->_x & type::_mask); }
		inline type operator & (const type b) const { return type(this->_x & b._x); }
		inline type operator | (const type b) const { return type(this->_x | b._x); }
		inline type operator ^ (const type b) const { return type(this->_x ^ b._x); }
		inline type& operator &=(const type b) { this->_x &= b._x; return *static_cast<type*>(this); }
		inline type& operator |=(const type b) { this->_x |= b._x; return *static_cast<type*>(this); }
		inline type& operator ^=(const type b) { this->_x ^= b._x; return *static_cast<type*>(this); }

		// String/flag conversions
		inline const char* const rawname() const // returns "None", a predefined name, or NULL
		{
			if (this->_x == 0)           { return "None"; }
			//if (this->_x == type::_mask) { return "All";  }
			const char* const* ns = type::_names();
			const base_type* vs = type::_values();
			for (size_t i = 0; i < _count; ++i) { if (this->_x == vs[i]) { return ns[i]; } }
			return NULL;
		}
		inline const std::string name() const // returns "None", a predefined name or a string composed of predefined names separated by |
		{
			const char* const n = this->rawname();
			if (n) { return n; }
			const char* const* ns = type::_names();
			const base_type* vs = type::_values();
			std::string s = "";
			for (size_t i = 0; i < _count; ++i) { if ((this->_x & vs[i]) == vs[i]) { if (s.size()) { s += '|'; } s += ns[i]; } }
			return s;
		}
		friend inline std::ostream& operator<< (std::ostream& s, const type& b) { return s << b.name(); }
		static inline type from_name(const char *name, size_t n) // reverses "name()" accepting predefined names and "None" separated by | while ignoring whitespace before and after each name - throws std::invalid_argument if a name is not valid
		{
			while (*name && std::isspace(*name)) { ++name; --n; } // remove leading whitespace
			while (n > 0 && std::isspace(name[n - 1])) { --n; } // remove training whitespace
			const char *pipe = (const char*)std::memchr(name, '|', n);
			if (pipe == NULL)
			{
				if (n == 4 && std::strncmp(name, "None", 4) == 0) { return type(/*0*/); }
				//if (n == 3 && std::strncmp(name, "All" , n)) { return type(type::_mask); }
				const char* const* ns = type::_names();
				const base_type* vs = type::_values();
				for (size_t i = 0; i < _count; ++i) { if (n == std::strlen(ns[i]) && std::strncmp(name, ns[i], n) == 0) { return type(vs[i]); } }
				throw std::invalid_argument("Flags name '" + std::string(name, n) + "' not found");
			}
			type b = type::from_name(name, pipe - name); n -= pipe - name; name = pipe + 1;
			while ((pipe = (const char*)std::memchr(name, '|', n)) != NULL)
			{
				b |= type::from_name(name, pipe - name); n -= pipe - name; name = pipe + 1;
			}
			b |= type::from_name(name, n - 1);
			return b;
		}
		static inline type from_name(const char * name)      { return type::from_name(name, std::strlen(name)); }
		static inline type from_name(const std::string name) { return type::from_name(name.c_str(), name.size()); }
	};
	template <class BT_, class FT_, size_t _count_> const FT_ _Flags_base<BT_, FT_, _count_>::None;
	//template <class BT_, class FT_, size_t _count_> const FT_ _Flags_base<BT_, FT_, _count_>::All(FT_::_mask);
	template <class BT_, class F_, BT_ V_> class NULL_FLAG_ { };
}

#pragma region Create a flags class for each size
// These macros are simply for building the various sized flags and are removed once they are made
// Currently the max number of flags is 8. To increase this:
//   Add additional LIST_#, LISTX_#, and POW_# macros
//   Increase value of MAX_FLAGS
//   Update all REM_# values and add additional REM_# values
//   Add additional FLAGS_PARTIAL(#); lines
//   Add an #undef at the end for each new define created

// LIST macro repeats an item multiple times with a separator
// ITEM is a macro that takes the 0-based index of the item, the 1-based index of the item, and the total length of the list
// SEP is a separator between items (once evaluated with ())
// N is the length of the list
#define LIST_1(ITEM, SEP, N) ITEM(0, 1, N)
#define LIST_2(ITEM, SEP, N) LIST_1(ITEM, SEP, N) SEP() ITEM(1, 2, N)
#define LIST_3(ITEM, SEP, N) LIST_2(ITEM, SEP, N) SEP() ITEM(2, 3, N)
#define LIST_4(ITEM, SEP, N) LIST_3(ITEM, SEP, N) SEP() ITEM(3, 4, N)
#define LIST_5(ITEM, SEP, N) LIST_4(ITEM, SEP, N) SEP() ITEM(4, 5, N)
#define LIST_6(ITEM, SEP, N) LIST_5(ITEM, SEP, N) SEP() ITEM(5, 6, N)
#define LIST_7(ITEM, SEP, N) LIST_6(ITEM, SEP, N) SEP() ITEM(6, 7, N)
#define LIST_8(ITEM, SEP, N) LIST_7(ITEM, SEP, N) SEP() ITEM(7, 8, N)
#define LIST_(ITEM, SEP, N) LIST_##N(ITEM, SEP, N) /* indirection layer */
#define LIST(ITEM, SEP, N) LIST_(ITEM, SEP, N)

// LISTX is very similar to LIST except that it uses a different ITEM macro (ITEM0) for the first list item and can be nested inside LIST
#define LISTX_1(ITEM0, ITEM, SEP, N) ITEM0(0, 1, N)
#define LISTX_2(ITEM0, ITEM, SEP, N) LISTX_1(ITEM0, ITEM, SEP, N) SEP() ITEM(1, 2, N)
#define LISTX_3(ITEM0, ITEM, SEP, N) LISTX_2(ITEM0, ITEM, SEP, N) SEP() ITEM(2, 3, N)
#define LISTX_4(ITEM0, ITEM, SEP, N) LISTX_3(ITEM0, ITEM, SEP, N) SEP() ITEM(3, 4, N)
#define LISTX_5(ITEM0, ITEM, SEP, N) LISTX_4(ITEM0, ITEM, SEP, N) SEP() ITEM(4, 5, N)
#define LISTX_6(ITEM0, ITEM, SEP, N) LISTX_5(ITEM0, ITEM, SEP, N) SEP() ITEM(5, 6, N)
#define LISTX_7(ITEM0, ITEM, SEP, N) LISTX_6(ITEM0, ITEM, SEP, N) SEP() ITEM(6, 7, N)
#define LISTX_8(ITEM0, ITEM, SEP, N) LISTX_7(ITEM0, ITEM, SEP, N) SEP() ITEM(7, 8, N)
#define LISTX_(ITEM0, ITEM, SEP, N) LISTX_##N(ITEM0, ITEM, SEP, N) /* indirection layer */
#define LISTX(ITEM0, ITEM, SEP, N) LISTX_(ITEM0, ITEM, SEP, N)

// Separator macros for use with LIST/LISTX
#define COMMA() ,
#define PIPE()  |
#define SEMICOLON() ;

// The powers of two
#define POW2_0 1
#define POW2_1 2
#define POW2_2 4
#define POW2_3 8
#define POW2_4 16
#define POW2_5 32
#define POW2_6 64
#define POW2_7 128

// The maximum number of flags (see above for how to increase this value)
#define MAX_FLAGS 8

// The value of MAX_FLAGS-X
#define REM_0  8
#define REM_1  7
#define REM_2  6
#define REM_3  5
#define REM_4  4
#define REM_5  3
#define REM_6  2
#define REM_7  1
#define REM_8  0
#define REM(X) REM_##X


// Items used with LIST/LISTX
#define FLG_TMPLT_ARG_ITEM(I0, I1, N)        F##I1##_
#define VAL_TMPLT_ARG_ITEM(I0, I1, N)        V##I1##_
#define FLG_TMPLT_ARG_W_TYPE_AND_DEF_ITEM(I0, I1, N) FC F##I1##_ = NF
#define FLG_TMPLT_ARG_W_TYPE_ITEM(I0, I1, N) FC  F##I1##_
#define VAL_TMPLT_ARG_W_TYPE_ITEM(I0, I1, N) BT_ V##I1##_
#define FLG_ASSIGN_ITEM(I0, I1, N)           n[I0] = F##I1##_<BT_, type, V##I1##_>::_name()
#define VAL_ASSIGN_ITEM(I0, I1, N)           v[I0] = V##I1##_
#define FLG_TMPLT_BASE_ITEM(I0, I1, N)       public F##I1##_<BT_, FLAGS_TYPE(N), V##I1##_>
#define POW2_ITEM(I0, I1, N)                 POW2_##I0
#define NF_ITEM(I0, I1, N)                   NF

#define NF       detail::NULL_FLAG_
#define FB(N)    detail::_Flags_base<BT_, FLAGS_TYPE(N), N>
#define FC       template <class,class,BT_> class
#define FLAGS_TYPE(N) typename Flags<BT_, LISTX(FLG_TMPLT_ARG_ITEM, FLG_TMPLT_ARG_ITEM, COMMA, N)>::WithValues<LISTX(VAL_TMPLT_ARG_ITEM, VAL_TMPLT_ARG_ITEM, COMMA, N)> /* uses LISTX so it can be inside LIST */
#define FLAGS_CONSTRUCTORS \
	public: \
		inline WithValues() : _Flags_base() { } \
		inline WithValues(const type& b) : _Flags_base(b) { } \
		inline explicit WithValues(base_type b) : _Flags_base(b) { } \
		inline type& operator=(const type& b) { this->_x = b._x; return *this; }
#define FLAGS_CORE(N) \
	{ \
	public: \
		template<LIST(VAL_TMPLT_ARG_W_TYPE_ITEM, COMMA, N)> \
		class WithValues : public FB(N), LIST(FLG_TMPLT_BASE_ITEM, COMMA, N) \
		{ \
			friend FB(N); \
			static const base_type _mask = LIST(VAL_TMPLT_ARG_ITEM, PIPE, N); \
			static void _fill_names(const char* n[N]) { LIST(FLG_ASSIGN_ITEM, SEMICOLON, N); } \
			static void _fill_values(base_type  v[N]) { LIST(VAL_ASSIGN_ITEM, SEMICOLON, N); } \
			FLAGS_CONSTRUCTORS \
		}; \
		typedef WithValues<LIST(POW2_ITEM, COMMA, N)> WithDefaultValues; \
	}
#define FLAGS_PARTIAL(N) \
	template<class BT_, LIST(FLG_TMPLT_ARG_W_TYPE_ITEM, COMMA, N)> \
	class Flags<BT_, LIST(FLG_TMPLT_ARG_ITEM, COMMA, N), LIST(NF_ITEM, COMMA, REM(N))> FLAGS_CORE(N)

template<class BT_, LISTX(FLG_TMPLT_ARG_W_TYPE_ITEM, FLG_TMPLT_ARG_W_TYPE_AND_DEF_ITEM, COMMA, MAX_FLAGS)> class Flags FLAGS_CORE(MAX_FLAGS);
FLAGS_PARTIAL(1);
FLAGS_PARTIAL(2);
FLAGS_PARTIAL(3);
FLAGS_PARTIAL(4);
FLAGS_PARTIAL(5);
FLAGS_PARTIAL(6);
FLAGS_PARTIAL(7);

#undef LIST_1
#undef LIST_2
#undef LIST_3
#undef LIST_4
#undef LIST_5
#undef LIST_6
#undef LIST_7
#undef LIST_8
#undef LIST_
#undef LIST
#undef LISTX_1
#undef LISTX_2
#undef LISTX_3
#undef LISTX_4
#undef LISTX_5
#undef LISTX_6
#undef LISTX_7
#undef LISTX_8
#undef LISTX_
#undef LISTX
#undef COMMA
#undef PIPE
#undef SEMICOLON
#undef POW2_0
#undef POW2_1
#undef POW2_2
#undef POW2_3
#undef POW2_4
#undef POW2_5
#undef POW2_6
#undef POW2_7
#undef REM_0
#undef REM_1
#undef REM_2
#undef REM_3
#undef REM_4
#undef REM_5
#undef REM_6
#undef REM_7
#undef REM_8
#undef REM
#undef FLG_TMPLT_ARG_ITEM
#undef VAL_TMPLT_ARG_ITEM
#undef FLG_TMPLT_ARG_W_TYPE_AND_DEF_ITEM
#undef FLG_TMPLT_ARG_W_TYPE_ITEM
#undef VAL_TMPLT_ARG_W_TYPE_ITEM
#undef FLG_ASSIGN_ITEM
#undef VAL_ASSIGN_ITEM
#undef FLG_TMPLT_BASE_ITEM
#undef POW2_ITEM
#undef NF_ITEM
#undef FC
#undef FB
#undef NF
#undef FLAGS_PARTIAL
#undef FLAGS_CORE
#undef FLAGS_CONSTRUCTORS
#undef FLAGS_TYPE
//#undef MAX_FLAGS
#pragma endregion

///////////////////////////////////////////////////////////////////////////////////////////////////
// Macro magic for making flags with default or defined values.
///////////////////////////////////////////////////////////////////////////////////////////////////

#define EXPAND(x) x

#define _xMAKE_FLAG_NAMES_1(F, ...)		FLAG_NAME(F)
#define _xMAKE_FLAG_NAMES_2(F, ...)		FLAG_NAME(F); EXPAND(_xMAKE_FLAG_NAMES_1(__VA_ARGS__))
#define _xMAKE_FLAG_NAMES_3(F, ...)		FLAG_NAME(F); EXPAND(_xMAKE_FLAG_NAMES_2(__VA_ARGS__))
#define _xMAKE_FLAG_NAMES_4(F, ...)		FLAG_NAME(F); EXPAND(_xMAKE_FLAG_NAMES_3(__VA_ARGS__))
#define _xMAKE_FLAG_NAMES_5(F, ...)		FLAG_NAME(F); EXPAND(_xMAKE_FLAG_NAMES_4(__VA_ARGS__))
#define _xMAKE_FLAG_NAMES_6(F, ...)		FLAG_NAME(F); EXPAND(_xMAKE_FLAG_NAMES_5(__VA_ARGS__))
#define _xMAKE_FLAG_NAMES_7(F, ...)		FLAG_NAME(F); EXPAND(_xMAKE_FLAG_NAMES_6(__VA_ARGS__))
#define _xMAKE_FLAG_NAMES_8(F, ...)		FLAG_NAME(F); EXPAND(_xMAKE_FLAG_NAMES_7(__VA_ARGS__))

#define _xFLAG_NAMES_1(F, ...)		Flag_##F
#define _xFLAG_NAMES_2(F, ...)		Flag_##F, EXPAND(_xFLAG_NAMES_1(__VA_ARGS__))
#define _xFLAG_NAMES_3(F, ...)		Flag_##F, EXPAND(_xFLAG_NAMES_2(__VA_ARGS__))
#define _xFLAG_NAMES_4(F, ...)		Flag_##F, EXPAND(_xFLAG_NAMES_3(__VA_ARGS__))
#define _xFLAG_NAMES_5(F, ...)		Flag_##F, EXPAND(_xFLAG_NAMES_4(__VA_ARGS__))
#define _xFLAG_NAMES_6(F, ...)		Flag_##F, EXPAND(_xFLAG_NAMES_5(__VA_ARGS__))
#define _xFLAG_NAMES_7(F, ...)		Flag_##F, EXPAND(_xFLAG_NAMES_6(__VA_ARGS__))
#define _xFLAG_NAMES_8(F, ...)		Flag_##F, EXPAND(_xFLAG_NAMES_7(__VA_ARGS__))

#define _xFLAG_VALUES_1(V)			V
#define _xFLAG_VALUES_2(V, ...)		V, EXPAND(_xFLAG_VALUES_1(__VA_ARGS__))
#define _xFLAG_VALUES_3(V, ...)		V, EXPAND(_xFLAG_VALUES_2(__VA_ARGS__))
#define _xFLAG_VALUES_4(V, ...)		V, EXPAND(_xFLAG_VALUES_3(__VA_ARGS__))
#define _xFLAG_VALUES_5(V, ...)		V, EXPAND(_xFLAG_VALUES_4(__VA_ARGS__))
#define _xFLAG_VALUES_6(V, ...)		V, EXPAND(_xFLAG_VALUES_5(__VA_ARGS__))
#define _xFLAG_VALUES_7(V, ...)		V, EXPAND(_xFLAG_VALUES_6(__VA_ARGS__))
#define _xFLAG_VALUES_8(V, ...)		V, EXPAND(_xFLAG_VALUES_7(__VA_ARGS__))

#define _xSKIP_FLAG_NAMES_1(N, F, ...)	EXPAND(_xFLAG_VALUES_##N(__VA_ARGS__))
#define _xSKIP_FLAG_NAMES_2(N, F, ...)	EXPAND(_xSKIP_FLAG_NAMES_1(N, __VA_ARGS__))
#define _xSKIP_FLAG_NAMES_3(N, F, ...)	EXPAND(_xSKIP_FLAG_NAMES_2(N, __VA_ARGS__))
#define _xSKIP_FLAG_NAMES_4(N, F, ...)	EXPAND(_xSKIP_FLAG_NAMES_3(N, __VA_ARGS__))
#define _xSKIP_FLAG_NAMES_5(N, F, ...)	EXPAND(_xSKIP_FLAG_NAMES_4(N, __VA_ARGS__))
#define _xSKIP_FLAG_NAMES_6(N, F, ...)	EXPAND(_xSKIP_FLAG_NAMES_5(N, __VA_ARGS__))
#define _xSKIP_FLAG_NAMES_7(N, F, ...)	EXPAND(_xSKIP_FLAG_NAMES_6(N, __VA_ARGS__))
#define _xSKIP_FLAG_NAMES_8(N, F, ...)	EXPAND(_xSKIP_FLAG_NAMES_7(N, __VA_ARGS__))

#define _xFLAGS_1(NAME, F)		       _xMAKE_FLAG_NAMES_1(F);             typedef ::Flags<std::uint8_t,       _xFLAG_NAMES_1(F)           >::WithDefaultValues NAME
#define _xFLAGS_2(NAME, ...)	EXPAND(_xMAKE_FLAG_NAMES_2(__VA_ARGS__));  typedef ::Flags<std::uint8_t,EXPAND(_xFLAG_NAMES_2(__VA_ARGS__))>::WithDefaultValues NAME
#define _xFLAGS_3(NAME, ...)	EXPAND(_xMAKE_FLAG_NAMES_3(__VA_ARGS__));  typedef ::Flags<std::uint8_t,EXPAND(_xFLAG_NAMES_3(__VA_ARGS__))>::WithDefaultValues NAME
#define _xFLAGS_4(NAME, ...)	EXPAND(_xMAKE_FLAG_NAMES_4(__VA_ARGS__));  typedef ::Flags<std::uint8_t,EXPAND(_xFLAG_NAMES_4(__VA_ARGS__))>::WithDefaultValues NAME
#define _xFLAGS_5(NAME, ...)	EXPAND(_xMAKE_FLAG_NAMES_5(__VA_ARGS__));  typedef ::Flags<std::uint8_t,EXPAND(_xFLAG_NAMES_5(__VA_ARGS__))>::WithDefaultValues NAME
#define _xFLAGS_6(NAME, ...)	EXPAND(_xMAKE_FLAG_NAMES_6(__VA_ARGS__));  typedef ::Flags<std::uint8_t,EXPAND(_xFLAG_NAMES_6(__VA_ARGS__))>::WithDefaultValues NAME
#define _xFLAGS_7(NAME, ...)	EXPAND(_xMAKE_FLAG_NAMES_7(__VA_ARGS__));  typedef ::Flags<std::uint8_t,EXPAND(_xFLAG_NAMES_7(__VA_ARGS__))>::WithDefaultValues NAME
#define _xFLAGS_8(NAME, ...)	EXPAND(_xMAKE_FLAG_NAMES_8(__VA_ARGS__));  typedef ::Flags<std::uint8_t,EXPAND(_xFLAG_NAMES_8(__VA_ARGS__))>::WithDefaultValues NAME

#define _xFLAGS_WITH_VALUES_2(NAME, BT, F, V)	       _xMAKE_FLAG_NAMES_1(F);             typedef ::Flags<BT,EXPAND(_xFLAG_NAMES_1(F)          >::WithValues<        _xSKIP_FLAG_NAMES_1(1,F,V)         > NAME
#define _xFLAGS_WITH_VALUES_4(NAME, BT, ...)	EXPAND(_xMAKE_FLAG_NAMES_2(__VA_ARGS__));  typedef ::Flags<BT,EXPAND(_xFLAG_NAMES_2(__VA_ARGS__))>::WithValues<EXPAND(_xSKIP_FLAG_NAMES_2(2,__VA_ARGS__))> NAME
#define _xFLAGS_WITH_VALUES_6(NAME, BT, ...)	EXPAND(_xMAKE_FLAG_NAMES_3(__VA_ARGS__));  typedef ::Flags<BT,EXPAND(_xFLAG_NAMES_3(__VA_ARGS__))>::WithValues<EXPAND(_xSKIP_FLAG_NAMES_3(3,__VA_ARGS__))> NAME
#define _xFLAGS_WITH_VALUES_8(NAME, BT, ...)	EXPAND(_xMAKE_FLAG_NAMES_4(__VA_ARGS__));  typedef ::Flags<BT,EXPAND(_xFLAG_NAMES_4(__VA_ARGS__))>::WithValues<EXPAND(_xSKIP_FLAG_NAMES_4(4,__VA_ARGS__))> NAME
#define _xFLAGS_WITH_VALUES_10(NAME, BT, ...)	EXPAND(_xMAKE_FLAG_NAMES_5(__VA_ARGS__));  typedef ::Flags<BT,EXPAND(_xFLAG_NAMES_5(__VA_ARGS__))>::WithValues<EXPAND(_xSKIP_FLAG_NAMES_5(5,__VA_ARGS__))> NAME
#define _xFLAGS_WITH_VALUES_12(NAME, BT, ...)	EXPAND(_xMAKE_FLAG_NAMES_6(__VA_ARGS__));  typedef ::Flags<BT,EXPAND(_xFLAG_NAMES_6(__VA_ARGS__))>::WithValues<EXPAND(_xSKIP_FLAG_NAMES_6(6,__VA_ARGS__))> NAME
#define _xFLAGS_WITH_VALUES_14(NAME, BT, ...)	EXPAND(_xMAKE_FLAG_NAMES_7(__VA_ARGS__));  typedef ::Flags<BT,EXPAND(_xFLAG_NAMES_7(__VA_ARGS__))>::WithValues<EXPAND(_xSKIP_FLAG_NAMES_7(7,__VA_ARGS__))> NAME
#define _xFLAGS_WITH_VALUES_16(NAME, BT, ...)	EXPAND(_xMAKE_FLAG_NAMES_8(__VA_ARGS__));  typedef ::Flags<BT,EXPAND(_xFLAG_NAMES_8(__VA_ARGS__))>::WithValues<EXPAND(_xSKIP_FLAG_NAMES_8(8,__VA_ARGS__))> NAME

// NUM_ARGS(...) evaluates to the literal number of the passed-in arguments.
#define _xNUM_ARGS2(X,X64,X63,X62,X61,X60,X59,X58,X57,X56,X55,X54,X53,X52,X51,X50,X49,X48,X47,X46,X45,X44,X43,X42,X41,X40,X39,X38,X37,X36,X35,X34,X33,X32,X31,X30,X29,X28,X27,X26,X25,X24,X23,X22,X21,X20,X19,X18,X17,X16,X15,X14,X13,X12,X11,X10,X9,X8,X7,X6,X5,X4,X3,X2,X1,N,...) N
#define NUM_ARGS(...) EXPAND(_xNUM_ARGS2(0,__VA_ARGS__ ,64,63,62,61,60,59,58,57,56,55,54,53,52,51,50,49,48,47,46,45,44,43,42,41,40,39,38,37,36,35,34,33,32,31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0))

#define _xFLAGS3(NAME, N, ...)	EXPAND(_xFLAGS_##N(NAME, __VA_ARGS__))
#define _xFLAGS2(NAME, N, ...)	_xFLAGS3(NAME, N, __VA_ARGS__)
#define FLAGS(NAME, ...)		_xFLAGS2(NAME, NUM_ARGS(__VA_ARGS__), __VA_ARGS__)

#define _xFLAGS_WITH_VALUES3(NAME, BT, N, ...)	EXPAND(_xFLAGS_WITH_VALUES_##N(NAME, BT, __VA_ARGS__))
#define _xFLAGS_WITH_VALUES2(NAME, BT, N, ...)	_xFLAGS_WITH_VALUES3(NAME, BT, N, __VA_ARGS__)
#define FLAGS_WITH_VALUES(NAME, BT, ...)		_xFLAGS_WITH_VALUES2(NAME, BT, NUM_ARGS(__VA_ARGS__), __VA_ARGS__)
