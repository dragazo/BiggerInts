#ifndef BIGGER_INTS_H
#define BIGGER_INTS_H

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <climits>
#include <cmath>
#include <type_traits>
#include <iostream>
#include <string>
#include <utility>
#include <exception>
#include <stdexcept>
#include <algorithm>

/*

Adds signed and unsigned types accessible by alias templates uint_t and int_t for any given number of bits.
These aliased types display behavior identical to the built-in integral types (moreover, uint_t<> and int_t<> for 8, 16, 32, or 64 bits returns a built-in type).
Any type returned has a full set of specialized type_traits and numeric_limits.
For a given number of bits, the signed/unsigned structures are identical (e.g. you can safely bind a reference from signed to unsigned with reinterpret_cast to avoid a (potentially expensive) copy).
Negative signed values are stored via 2's complement, just as with built-in integral types on anything made since the dark ages.

Aside from these two alias templates, everything else is subject to change/removal and should not be used directly.

Report bugs to https://github.com/dragazo/BiggerInts/issues

*/

// -- optimization switches -- //

// the types to use for zero/sign extension. zero ext must be unsigned. sign ext must be signed.
#define ZEXT_TYPE unsigned
#define SEXT_TYPE signed

// by default, a template expression is generated for each operataion, which can be significantly faster than looping.
// however, at some point the compiler will decline to inline things, at which point it becomes a normal recursive function, which is significantly slower than looping.
// the following macros control the threshholds for selecting loop algorithms vs template expression algorithms at compile time.
// at and above the specified bit count will use loops.
// if you want to squeeze out as much speed as possible, fiddling with these values can significantly impact performance.

#define BOOL_LOOP_THRESH 1024

#define INCDEC_LOOP_THRESH 1024

#define ADDSUB_LOOP_THRESH 1024

#define BITWISE_LOOP_THRESH 1024

#define SHIFT_LOOP_THRESH 1024

#define DIVMOD_BLOCK_SKIP_THRESH 1024

#define CMP_LOOP_THRESH 1024

#define NOT_LOOP_THRESH 1024

// -- utility macros -- //

// shorthand used for SFINAE overloading - returns true if the type is compact
#define COMPACT(sign) (sizeof(double_int<bits, sign>) == bits / CHAR_BIT)

// -------------------- //

namespace BiggerInts
{
	// -- helpers and types -- //

	typedef std::uint64_t u64;
	typedef std::int64_t i64;

	// returns true if val is a power of 2
	inline constexpr bool is_pow2(u64 val) noexcept { return val != 0 && (val & (val - 1)) == 0; }

	// given a size in bits, returns a power of 2 (also in bits) large enough to contain it and that is no smaller than 8
	inline constexpr u64 round_bits_up(u64 size)
	{
		if (size < 8) return 8;
		if (size > 0x8000000000000000ul) throw std::domain_error("bit size was too large");
		while (!is_pow2(size)) size = (size | (size - 1)) + 1;
		return size;
	}

	// -- #hax typedef omega synergy highjinks - don't even need to be defined, boi -- //
	
	template<typename T, u64 bits> struct masked_single_int;
	template<u64 bits, bool sign> struct double_int;

	template<u64 bits, bool sign, std::enable_if_t<(bits > 64 && is_pow2(bits)), int> = 0> double_int<bits, sign> returns_proper_type();
	template<u64 bits, bool sign, std::enable_if_t<(bits > 64 && !is_pow2(bits)), int> = 0> masked_single_int<double_int<round_bits_up(bits), sign>, bits> returns_proper_type();

	template<u64 bits, bool sign, std::enable_if_t<(bits == 64), int> = 0> std::conditional_t<sign, std::int64_t, std::uint64_t> returns_proper_type();
	template<u64 bits, bool sign, std::enable_if_t<(bits > 32 && bits < 64), int> = 0> masked_single_int<std::conditional_t<sign, std::int64_t, std::uint64_t>, bits> returns_proper_type();

	template<u64 bits, bool sign, std::enable_if_t<(bits == 32), int> = 0> std::conditional_t<sign, std::int32_t, std::uint32_t> returns_proper_type();
	template<u64 bits, bool sign, std::enable_if_t<(bits > 16 && bits < 32), int> = 0> masked_single_int<std::conditional_t<sign, std::int32_t, std::uint32_t>, bits> returns_proper_type();

	template<u64 bits, bool sign, std::enable_if_t<(bits == 16), int> = 0> std::conditional_t<sign, std::int16_t, std::uint16_t> returns_proper_type();
	template<u64 bits, bool sign, std::enable_if_t<(bits > 8 && bits < 16), int> = 0> masked_single_int<std::conditional_t<sign, std::int16_t, std::uint16_t>, bits> returns_proper_type();
	
	template<u64 bits, bool sign, std::enable_if_t<(bits == 8), int> = 0> std::conditional_t<sign, std::int8_t, std::uint8_t> returns_proper_type();
	template<u64 bits, bool sign, std::enable_if_t<(bits < 8), int> = 0> masked_single_int<std::conditional_t<sign, std::int8_t, std::uint8_t>, bits> returns_proper_type();

	template<u64 bits> using uint_t = decltype(returns_proper_type<bits, false>());
	template<u64 bits> using int_t = decltype(returns_proper_type<bits, true>());

	// -- sigh, back to normal templating :( -- //
	
	// contains an integral constant that represents the number of bits in a given (supported) (integral) type
	template<typename T> struct bit_count {};

	template<> struct bit_count<unsigned char> : std::integral_constant<u64, sizeof(unsigned char) * CHAR_BIT> {};
	template<> struct bit_count<unsigned short> : std::integral_constant<u64, sizeof(unsigned short) * CHAR_BIT> {};
	template<> struct bit_count<unsigned int> : std::integral_constant<u64, sizeof(unsigned int) * CHAR_BIT> {};
	template<> struct bit_count<unsigned long> : std::integral_constant<u64, sizeof(unsigned long) * CHAR_BIT> {};
	template<> struct bit_count<unsigned long long> : std::integral_constant<u64, sizeof(unsigned long long) * CHAR_BIT> {};

	template<> struct bit_count<signed char> : std::integral_constant<u64, sizeof(signed char) * CHAR_BIT> {};
	template<> struct bit_count<signed short> : std::integral_constant<u64, sizeof(signed short) * CHAR_BIT> {};
	template<> struct bit_count<signed int> : std::integral_constant<u64, sizeof(signed int) * CHAR_BIT> {};
	template<> struct bit_count<signed long> : std::integral_constant<u64, sizeof(signed long) * CHAR_BIT> {};
	template<> struct bit_count<signed long long> : std::integral_constant<u64, sizeof(signed long long) * CHAR_BIT> {};

	template<typename T, u64 bits> struct bit_count<masked_single_int<T, bits>> : std::integral_constant<u64, bits> {};
	template<u64 bits, bool sign> struct bit_count<double_int<bits, sign>> : std::integral_constant<u64, bits> {};
	
	// -- container impl -- //

	// holds an integral type T that, upon assignment, is masked to the specified number of bits. bits must be in the range (bit count T / 2, bit count T).
	template<typename T, u64 bits> struct masked_single_int
	{
	private: // -- data -- //

		static_assert(bits > bit_count<T>::value / 2 && bits < bit_count<T>::value, "masked_single_int: bit count out of intended range");

		static constexpr T mask = ((T)1 << bits) - 1; // mask to apply - shift can't overflow because bits < bit_count<T>::count

		T val; // the actual stored value

	public: // -- wrapper interface -- // 

		constexpr operator T() const { return val; }

		constexpr masked_single_int(T _val = 0) noexcept : val(_val & mask) {}
		constexpr masked_single_int &operator=(T _val) { val = _val & mask; return *this; }

		constexpr masked_single_int &operator++() { val = (val + 1) & mask; return *this; }
		constexpr masked_single_int &operator--() { val = (val - 1) & mask; return *this; }

		constexpr T operator++(int) { T ret = val; val = (val + 1) & mask; return ret; }
		constexpr T operator--(int) { T ret = val; val = (val - 1) & mask; return ret; }

		constexpr masked_single_int &operator+=(T _val) { val = (val + _val) & mask; return *this; }
		constexpr masked_single_int &operator-=(T _val) { val = (val + _val) & mask; return *this; }

		constexpr masked_single_int &operator*=(T _val) { val = (val * _val) & mask; return *this; }
		constexpr masked_single_int &operator/=(T _val) { val = (val / _val) & mask; return *this; }
		constexpr masked_single_int &operator%=(T _val) { val = (val % _val) & mask; return *this; }

		constexpr masked_single_int &operator&=(T _val) { val = (val & _val) & mask; return *this; }
		constexpr masked_single_int &operator|=(T _val) { val = (val | _val) & mask; return *this; }
		constexpr masked_single_int &operator^=(T _val) { val = (val ^ _val) & mask; return *this; }

		constexpr masked_single_int &operator<<=(T _val) { val = (val << _val) & mask; return *this; }
		constexpr masked_single_int &operator>>=(T _val) { val = (val >> _val) & mask; return *this; }
	};
	
	// given a desired number of bits (power of 2), simulates it with 2 ints of half that size - sign flag marks signed/unsigned
	template<u64 bits, bool sign> struct double_int
	{
	public: // -- data -- //

		static_assert(is_pow2(bits) && bits >= 128, "double_int: bits was out of valid domain");

		typedef uint_t<bits / 2> half;

		half low, high;

	public: // -- core -- //

		inline constexpr double_int() = default;

		inline constexpr double_int(const double_int&) = default;
		inline constexpr double_int &operator=(const double_int&) = default;

		inline constexpr explicit operator double_int<bits, !sign>() const noexcept { return *(double_int<bits, !sign>*)this; }

	public: // -- promotion constructors -- //

		inline constexpr double_int(unsigned short val) noexcept : high((ZEXT_TYPE)0), low(val) {}
		inline constexpr double_int(unsigned int val) noexcept : high((ZEXT_TYPE)0), low(val) {}
		inline constexpr double_int(unsigned long val) noexcept : high((ZEXT_TYPE)0), low(val) {}
		inline constexpr double_int(unsigned long long val) noexcept : high((ZEXT_TYPE)0), low(val) {}
		
		inline constexpr double_int(signed short val) noexcept : high(val < 0 ? (SEXT_TYPE)-1 : (SEXT_TYPE)0), low(val) {}
		inline constexpr double_int(signed int val) noexcept : high(val < 0 ? (SEXT_TYPE)-1 : (SEXT_TYPE)0), low(val) {}
		inline constexpr double_int(signed long val) noexcept : high(val < 0 ? (SEXT_TYPE)-1 : (SEXT_TYPE)0), low(val) {}
		inline constexpr double_int(signed long long val) noexcept : high(val < 0 ? (SEXT_TYPE)-1 : (SEXT_TYPE)0), low(val) {}

		template<u64 _bits, std::enable_if_t<(bits > _bits), int> = 0> inline constexpr double_int(const double_int<_bits, false> &val) noexcept : high((ZEXT_TYPE)0), low(val) {}
		template<u64 _bits, std::enable_if_t<(bits > _bits), int> = 0> inline constexpr double_int(const double_int<_bits, true> &val) noexcept : high(is_neg(val) ? (SEXT_TYPE)-1 : (SEXT_TYPE)0), low(val) {}

		template<u64 _bits, bool _sign, std::enable_if_t<(bits == _bits), int> = 0> inline constexpr double_int(const double_int<_bits, _sign> &val) noexcept { *this = *(double_int*)&val; }

	public: // -- promotion assignment -- //

		inline constexpr double_int &operator=(unsigned short val) noexcept { high = (ZEXT_TYPE)0; low = val; return *this; }
		inline constexpr double_int &operator=(unsigned int val) noexcept { high = (ZEXT_TYPE)0; low = val; return *this; }
		inline constexpr double_int &operator=(unsigned long val) noexcept { high = (ZEXT_TYPE)0; low = val; return *this; }
		inline constexpr double_int &operator=(unsigned long long val) noexcept { high = (ZEXT_TYPE)0; low = val; return *this; }

		inline constexpr double_int &operator=(signed short val) noexcept { high = val < 0 ? (SEXT_TYPE)-1 : (SEXT_TYPE)0; low = val; return *this; }
		inline constexpr double_int &operator=(signed int val) noexcept { high = val < 0 ? (SEXT_TYPE)-1 : (SEXT_TYPE)0; low = val; return *this; }
		inline constexpr double_int &operator=(signed long val) noexcept { high = val < 0 ? (SEXT_TYPE)-1 : (SEXT_TYPE)0; low = val; return *this; }
		inline constexpr double_int &operator=(signed long long val) noexcept { high = val < 0 ? (SEXT_TYPE)-1 : (SEXT_TYPE)0; low = val; return *this; }
		
		template<u64 _bits, std::enable_if_t<(bits > _bits), int> = 0> inline constexpr double_int &operator=(const double_int<_bits, false> &val) noexcept { high = (ZEXT_TYPE)0; low = val; return *this; }
		template<u64 _bits, std::enable_if_t<(bits > _bits), int> = 0> inline constexpr double_int &operator=(const double_int<_bits, true> &val) noexcept { high = is_neg(val) ? (SEXT_TYPE)-1 : (SEXT_TYPE)0; low = val; return *this; }

		template<u64 _bits, bool _sign, std::enable_if_t<(bits == _bits), int> = 0> inline constexpr double_int &operator=(const double_int<_bits, _sign> &val) noexcept { *this = *(double_int*)&val; return *this; }

	public: // -- demotion conversion -- //

		inline constexpr explicit operator unsigned short() const noexcept { return (unsigned short)low64(*this); }
		inline constexpr explicit operator unsigned int() const noexcept { return (unsigned int)low64(*this); }
		inline constexpr explicit operator unsigned long() const noexcept { return (unsigned long)low64(*this); }
		inline constexpr explicit operator unsigned long long() const noexcept { return (unsigned long long)low64(*this); }

		inline constexpr explicit operator signed short() const noexcept { return (signed short)low64(*this); }
		inline constexpr explicit operator signed int() const noexcept { return (signed int)low64(*this); }
		inline constexpr explicit operator signed long() const noexcept { return (signed long)low64(*this); }
		inline constexpr explicit operator signed long long() const noexcept { return (signed long long)low64(*this); }

		template<u64 _bits, bool _sign, std::enable_if_t<(_bits < bits), int> = 0>
		inline constexpr explicit operator double_int<_bits, _sign>() const noexcept
		{
			// we're casting to a smaller type and low is guaranteed to be at offset 0 for all forms, so we can just reinterpret this
			return *(double_int<_bits, _sign>*)this;
		}

	public: // -- bool conversion -- //

		inline constexpr explicit operator bool() const noexcept { return to_bool(*this); }
		inline constexpr friend bool operator!(const double_int &a) noexcept { return !to_bool(a); }
	};

	// shorthand for binary ops +, -, *, etc. will refer to a form that uses both sides of double_int
	#define SHORTHAND_BINARY_FORMATTER(op, sign, type) \
		template<u64 bits> inline constexpr double_int<bits, sign> operator op(const double_int<bits, sign> &a, type b) noexcept { return a op (double_int<bits, sign>)b; } \
		template<u64 bits> inline constexpr double_int<bits, sign> operator op(type a, const double_int<bits, sign> &b) noexcept { return (double_int<bits, sign>)a op b; }
	// does all the standard shorthands for a given op
	#define SHORTERHAND_BINARY_FORMATTER(op) \
		SHORTHAND_BINARY_FORMATTER(op, false, unsigned short) \
		SHORTHAND_BINARY_FORMATTER(op, false, unsigned int) \
		SHORTHAND_BINARY_FORMATTER(op, false, unsigned long) \
		SHORTHAND_BINARY_FORMATTER(op, false, unsigned long long) \
		\
		SHORTHAND_BINARY_FORMATTER(op, true, signed short) \
		SHORTHAND_BINARY_FORMATTER(op, true, signed int) \
		SHORTHAND_BINARY_FORMATTER(op, true, signed long) \
		SHORTHAND_BINARY_FORMATTER(op, true, signed long long)

	// -- inc/dec -- //

	template<u64 bits, bool sign, std::enable_if_t<COMPACT(sign), int> = 0> inline constexpr double_int<bits, sign> &operator++(double_int<bits, sign> &a) noexcept
	{
		// wind up 64-bit blocks - break on the first that doesn't overflow
		u64 *ptr = (u64*)&a;
		for (u64 i = 0; i < bits / 64; ++i)
			if (++ptr[i]) break;
		return a;
	}
	template<u64 bits, bool sign, std::enable_if_t<!COMPACT(sign), int> = 0> inline constexpr double_int<bits, sign> &operator++(double_int<bits, sign> &a) noexcept
	{
		if (!++a.low) ++a.high;
		return a;
	}

	template<u64 bits, bool sign, std::enable_if_t<COMPACT(sign), int> = 0> inline constexpr double_int<bits, sign> &operator--(double_int<bits, sign> &a) noexcept
	{
		// wind up 64-bit blocks - break on the first that doesn't underflow
		u64 *ptr = (u64*)&a;
		for (u64 i = 0; i < bits / 64; ++i)
			if (ptr[i]--) break;
		return a;
	}
	template<u64 bits, bool sign, std::enable_if_t<!COMPACT(sign), int> = 0> inline constexpr double_int<bits, sign> &operator--(double_int<bits, sign> &a) noexcept
	{
		if (!a.low) --a.high;
		--a.low;
		return a;
	}

	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> operator++(double_int<bits, sign> &a, int) noexcept { double_int<bits, sign> val = a; ++a; return val; }
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> operator--(double_int<bits, sign> &a, int) noexcept { double_int<bits, sign> val = a; --a; return val; }

	// -- add -- //

	template<u64 bits, bool sign1, bool sign2> inline constexpr auto operator+=(double_int<bits, sign1> &a, const double_int<bits, sign2> &b) noexcept -> std::enable_if_t<(COMPACT(sign1) && COMPACT(sign2) && bits >= ADDSUB_LOOP_THRESH), double_int<bits, sign1>&>
	{
		u64 *dest = (u64*)&a, *src = (u64*)&b;
		u64 i, j;

		for (i = 0; i < bits / 64; ++i)
		{
			dest[i] += src[i];
			if (dest[i] < src[i])
			{
				for (j = i + 1; j < bits / 64; ++j)
					if (++dest[j]) break;
			}
		}
		
		return a;
	}
	template<u64 bits, bool sign1, bool sign2> inline constexpr auto operator+=(double_int<bits, sign1> &a, const double_int<bits, sign2> &b) noexcept -> std::enable_if_t<(!COMPACT(sign1) || !COMPACT(sign2) || bits < ADDSUB_LOOP_THRESH), double_int<bits, sign1>&>
	{
		a.low += b.low;
		if (a.low < b.low) ++a.high;
		a.high += b.high;
		return a;
	}

	template<u64 bits, bool sign, typename T> inline constexpr double_int<bits, sign> &operator+=(double_int<bits, sign> &a, const T &b) noexcept { a += (double_int<bits, sign>)b; return a; }

	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> operator+(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { double_int<bits, sign> res = a; res += b; return res; }

	SHORTERHAND_BINARY_FORMATTER(+)

	// -- sub -- //

	template<u64 bits, bool sign1, bool sign2> inline constexpr auto operator-=(double_int<bits, sign1> &a, const double_int<bits, sign2> &b) noexcept -> std::enable_if_t<(COMPACT(sign1) && COMPACT(sign2) && bits >= ADDSUB_LOOP_THRESH), double_int<bits, sign1>&>
	{
		u64 *dest = (u64*)&a, *src = (u64*)&b;
		u64 i, j;

		for (i = 0; i < bits / 64; ++i)
		{
			if (dest[i] < src[i])
			{
				for (j = i + 1; j < bits / 64; ++j)
					if (dest[j]--) break;
			}
			dest[i] -= src[i];
		}

		return a;
	}
	template<u64 bits, bool sign1, bool sign2> inline constexpr auto operator-=(double_int<bits, sign1> &a, const double_int<bits, sign2> &b) noexcept -> std::enable_if_t<(!COMPACT(sign1) || !COMPACT(sign2) || bits < ADDSUB_LOOP_THRESH), double_int<bits, sign1>&>
	{
		if (a.low < b.low) --a.high;
		a.low -= b.low;
		a.high -= b.high;
		return a;
	}
	

	template<u64 bits, bool sign, typename T> inline constexpr double_int<bits, sign> &operator-=(double_int<bits, sign> &a, const T &b) noexcept { a -= (double_int<bits, sign>)b; return a; }

	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> operator-(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { double_int<bits, sign> res = a; res -= b; return res; }

	SHORTERHAND_BINARY_FORMATTER(-)

	// -- and -- //

	template<u64 bits, bool sign1, bool sign2> inline constexpr auto operator&=(double_int<bits, sign1> &a, const double_int<bits, sign2> &b) noexcept -> std::enable_if_t<(COMPACT(sign1) && COMPACT(sign2) && bits >= BITWISE_LOOP_THRESH), double_int<bits, sign1>&>
	{
		u64 *dest = (u64*)&a, *src = (u64*)&b;
		for (u64 i = 0; i < bits / 64; ++i) dest[i] &= src[i];
		return a;
	}
	template<u64 bits, bool sign1, bool sign2> inline constexpr auto operator&=(double_int<bits, sign1> &a, const double_int<bits, sign2> &b) noexcept -> std::enable_if_t<(!COMPACT(sign1) || !COMPACT(sign2) || bits < BITWISE_LOOP_THRESH), double_int<bits, sign1>&>
	{
		a.low &= b.low;
		a.high &= b.high;
		return a;
	}

	template<u64 bits, bool sign, typename T> inline constexpr double_int<bits, sign> &operator&=(double_int<bits, sign> &a, const T &b) noexcept { a &= (double_int<bits, sign>)b; return a; }

	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator&=(double_int<bits, sign> &a, unsigned short b) noexcept { a = low64(a) & b; return a; }
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator&=(double_int<bits, sign> &a, unsigned int b) noexcept { a = low64(a) & b; return a; }
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator&=(double_int<bits, sign> &a, unsigned long b) noexcept { a = low64(a) & b; return a; }
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator&=(double_int<bits, sign> &a, unsigned long long b) noexcept { a = low64(a) & b; return a; }
	
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> operator&(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { double_int<bits, sign> res = a; res &= b; return res; }

	SHORTERHAND_BINARY_FORMATTER(&)

	// -- or -- //

	template<u64 bits, bool sign1, bool sign2> inline constexpr auto operator|=(double_int<bits, sign1> &a, const double_int<bits, sign2> &b) noexcept -> std::enable_if_t<(COMPACT(sign1) && COMPACT(sign2) && bits >= BITWISE_LOOP_THRESH), double_int<bits, sign1>&>
	{
		u64 *dest = (u64*)&a, *src = (u64*)&b;
		for (u64 i = 0; i < bits / 64; ++i) dest[i] |= src[i];
		return a;
	}
	template<u64 bits, bool sign1, bool sign2> inline constexpr auto operator|=(double_int<bits, sign1> &a, const double_int<bits, sign2> &b) noexcept->std::enable_if_t<(!COMPACT(sign1) || !COMPACT(sign2) || bits < BITWISE_LOOP_THRESH), double_int<bits, sign1>&>
	{
		a.low |= b.low;
		a.high |= b.high;
		return a;
	}

	template<u64 bits, bool sign, typename T> inline constexpr double_int<bits, sign> &operator|=(double_int<bits, sign> &a, const T &b) noexcept { a |= (double_int<bits, sign>)b; return a; }

	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator|=(double_int<bits, sign> &a, unsigned short b) noexcept { ref_low64(a) |= b; return a; }
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator|=(double_int<bits, sign> &a, unsigned int b) noexcept { ref_low64(a) |= b; return a; }
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator|=(double_int<bits, sign> &a, unsigned long b) noexcept { ref_low64(a) |= b; return a; }
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator|=(double_int<bits, sign> &a, unsigned long long b) noexcept { ref_low64(a) |= b; return a; }

	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> operator|(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { double_int<bits, sign> res = a; res |= b; return res; }

	SHORTERHAND_BINARY_FORMATTER(|)

	// -- xor -- //

	template<u64 bits, bool sign1, bool sign2> inline constexpr double_int<bits, sign1> &operator^=(double_int<bits, sign1> &a, const double_int<bits, sign2> &b) noexcept
	{
		if constexpr (COMPACT(sign1) && COMPACT(sign2) && bits >= BITWISE_LOOP_THRESH)
		{
			u64 *dest = (u64*)&a, *src = (u64*)&b;
			for (u64 i = 0; i < bits / 64; ++i) dest[i] ^= src[i];
			return a;
		}
		else
		{
			a.low ^= b.low;
			a.high ^= b.high;
			return a;
		}
	}

	template<u64 bits, bool sign, typename T> inline constexpr double_int<bits, sign> &operator^=(double_int<bits, sign> &a, const T &b) noexcept { a ^= (double_int<bits, sign>)b; return a; }

	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator^=(double_int<bits, sign> &a, unsigned short b) noexcept { ref_low64(a) ^= b; return a; }
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator^=(double_int<bits, sign> &a, unsigned int b) noexcept { ref_low64(a) ^= b; return a; }
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator^=(double_int<bits, sign> &a, unsigned long b) noexcept { ref_low64(a) ^= b; return a; }
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator^=(double_int<bits, sign> &a, unsigned long long b) noexcept { ref_low64(a) ^= b; return a; }

	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> operator^(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { auto res = a; res ^= b; return res; }

	SHORTERHAND_BINARY_FORMATTER(^)

	// -- shl/sal -- //

	template<u64 bits, bool sign, std::enable_if_t<(COMPACT(sign) && bits >= SHIFT_LOOP_THRESH), int> = 0> inline constexpr double_int<bits, sign> &operator<<=(double_int<bits, sign> &val, u64 count) noexcept
	{
		u64 *ptr = (u64*)&val;

		// count masked to bits limit - no-op on zero
		if (count &= bits - 1)
		{
			// get number of blocks being traversed
			u64 blocks = count / 64;
			// if that's not zero
			if (blocks != 0)
			{
				// shuffle blocks
				for (u64 i = (bits / 64) - 1; i >= blocks; --i) ptr[i] = ptr[i - blocks];
				// fill the rest with zero
				for (u64 i = 0; i < blocks; ++i) ptr[i] = 0;
			}
			// count masked to block limit - no-op on zero
			if (count &= 0x3f)
			{
				// perform the intra-block shifting
				for (u64 i = (bits / 64) - 1; i > blocks; --i)
				{
					ptr[i] <<= count;
					ptr[i] |= ptr[i - 1] >> (64 - count);
				}
				// perform the last block separately (avoids an out of bounds error where blocks = 0)
				ptr[blocks] <<= count;
			}
		}

		return val;
	}
	template<u64 bits, bool sign, std::enable_if_t<(!COMPACT(sign) || bits < SHIFT_LOOP_THRESH), int> = 0> inline constexpr double_int<bits, sign> &operator<<=(double_int<bits, sign> &val, u64 count) noexcept
	{
		// count masked to bits limit - no-op on zero
		if (count &= bits - 1)
		{
			if (count >= bits / 2)
			{
				val.high = val.low;
				val.high <<= count - bits / 2;
				val.low = (ZEXT_TYPE)0;
			}
			else
			{
				val.high <<= count;
				val.high |= val.low >> (bits / 2 - count);
				val.low <<= count;
			}
		}
		return val;
	}
	
	template<u64 bits, bool sign, typename T> inline constexpr double_int<bits, sign> operator<<(const double_int<bits, sign> &a, const T &count) noexcept { double_int<bits, sign> res = a; res <<= (u64)count; return res; }

	// -- shr/sar -- //

	// shr
	template<u64 bits> inline constexpr auto operator>>=(double_int<bits, false> &val, u64 count) noexcept -> std::enable_if_t<(COMPACT(false) && bits >= SHIFT_LOOP_THRESH), double_int<bits, false>&>
	{
		u64 *ptr = (u64*)&val;

		// count masked to bits limit - no-op on zero
		if (count &= bits - 1)
		{
			// get number of blocks being traversed
			u64 blocks = count / 64;
			// if that's not zero
			if (blocks != 0)
			{
				// shuffle blocks
				for (u64 i = 0; i < (bits / 64) - blocks; ++i) ptr[i] = ptr[i + blocks];
				// fill the rest with zero
				for (u64 i = (bits / 64) - blocks; i < (bits / 64); ++i) ptr[i] = 0;
			}
			// count masked to block limit - no-op on zero
			if (count &= 0x3f)
			{
				// perform the intra-block shifting
				for (u64 i = 0; i < (bits / 64) - blocks - 1; ++i)
				{
					ptr[i] >>= count;
					ptr[i] |= ptr[i + 1] << (64 - count);
				}
				// perform the last block separately (avoids an out of bounds error where blocks = 0)
				ptr[(bits / 64) - blocks - 1] >>= count;
			}
		}

		return val;
	}
	template<u64 bits> inline constexpr auto operator>>=(double_int<bits, false> &val, u64 count) noexcept -> std::enable_if_t<(!COMPACT(false) || bits < SHIFT_LOOP_THRESH), double_int<bits, false>&>
	{
		// count masked to bits limit - no-op on zero
		if (count &= bits - 1)
		{
			if (count >= bits / 2)
			{
				val.low = val.high;
				val.low >>= count - bits / 2;
				val.high = (ZEXT_TYPE)0;
			}
			else
			{
				val.low >>= count;
				val.low |= val.high << (bits / 2 - count);
				val.high >>= count;
			}
		}
		return val;
	}

	// sar
	template<u64 bits> inline constexpr auto operator>>=(double_int<bits, true> &val, u64 count) noexcept -> std::enable_if_t<(COMPACT(true) && bits >= SHIFT_LOOP_THRESH), double_int<bits, true>&>
	{
		u64 *ptr = (u64*)&val;
		u64 fill = is_neg(val) ? 0xffffffffffffffff : 0; // value to fill in for sign filling

		// count masked to bits limit - no-op on zero
		if (count &= bits - 1)
		{
			// get number of blocks being traversed
			u64 blocks = count / 64;
			// if that's not zero
			if (blocks != 0)
			{
				// shuffle blocks
				for (u64 i = 0; i < (bits / 64) - blocks; ++i) ptr[i] = ptr[i + blocks];
				// fill the rest with fill
				for (u64 i = (bits / 64) - blocks; i < (bits / 64); ++i) ptr[i] = fill;
			}
			// count masked to block limit - no-op on zero
			if (count &= 0x3f)
			{
				// perform the intra-block shifting
				for (u64 i = 0; i < (bits / 64) - blocks - 1; ++i)
				{
					ptr[i] >>= count;
					ptr[i] |= ptr[i + 1] << (64 - count);
				}
				// perform the last block separately (avoids an out of bounds error where blocks = 0)
				ptr[(bits / 64) - blocks - 1] >>= count;
				ptr[(bits / 64) - blocks - 1] |= fill << (64 - count); // last block must account for fill value
			}
		}

		return val;
	}
	template<u64 bits> inline constexpr auto operator>>=(double_int<bits, true> &val, u64 count) noexcept -> std::enable_if_t<(!COMPACT(true) || bits < SHIFT_LOOP_THRESH), double_int<bits, true>&>
	{
		// count masked to bits limit - no-op on zero
		if (count &= bits - 1)
		{
			if (count >= bits / 2)
			{
				val.low = val.high;
				val.low >>= count - bits / 2;
				if (is_neg(val.high)) // fill with sign bit
				{
					val.high = (SEXT_TYPE)-1;
					if (count != bits / 2) val.low |= val.high << (bits - count);
				}
				else val.high = 0;
			}
			else
			{
				val.low >>= count;
				val.low |= val.high << (bits / 2 - count);
				if (is_neg(val.high)) // fill with sign bit
				{
					val.high >>= count;
					val.high |= (decltype(val.high))((SEXT_TYPE)-1) << (bits / 2 - count);
				}
				else val.high >>= count;
			}
		}
		return val;
	}

	template<u64 bits, bool sign, typename T> inline constexpr double_int<bits, sign> operator>>(const double_int<bits, sign> &a, const T &count) noexcept { double_int<bits, sign> res = a; res >>= (u64)count; return res; }

	// -- unary -- //

	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> operator+(const double_int<bits, sign> &a) noexcept { return a; }

	template<u64 bits> inline constexpr double_int<bits, true> operator-(const double_int<bits, true> &a) noexcept { double_int<bits, true> res = a; make_neg(res); return res; }

	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> operator~(const double_int<bits, sign> &a) noexcept { double_int<bits, sign> res = a; make_not(res); return res; }

	// -- mul -- //

	template<u64 bits> inline constexpr double_int<bits, false> operator*(const double_int<bits, false> &a, const double_int<bits, false> &b) noexcept
	{
		double_int<bits, false> res = (ZEXT_TYPE)0, _a = a;
		for (u64 bit = 0; _a; ++bit, _a <<= 1)
			if (bit_test(b, bit)) res += _a;
		return res;
	}
	template<u64 bits> inline constexpr double_int<bits, true> operator*(const double_int<bits, true> &a, const double_int<bits, true> &b) noexcept
	{
		// we'll do signed divmod in terms of unsigned

		double_int<bits, false> ua = *(double_int<bits, false>*)&a;
		double_int<bits, false> ub = *(double_int<bits, false>*)&b;

		bool neg = false;
		if (bit_test(ua, bits - 1)) { make_neg(ua); neg = !neg; }
		if (bit_test(ub, bits - 1)) { make_neg(ub); neg = !neg; }

		ua = ua * ub;

		if (neg) make_neg(ua);

		return *(double_int<bits, true>*)&ua;
	}

	SHORTERHAND_BINARY_FORMATTER(*)

	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator*=(double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { a = a * b; return a; }
	template<u64 bits, bool sign, typename T> inline constexpr double_int<bits, sign> &operator*=(double_int<bits, sign> &a, const T &b) noexcept { a *= (double_int<bits, sign>)b; return a; }

	// -- divmod -- //

	// helper for divmod_unchecked - takes the bit to start with - don't use this directly
	template<u64 bits> inline constexpr std::pair<double_int<bits, false>, double_int<bits, false>> __divmod_unchecked(const double_int<bits, false> &num, const double_int<bits, false> &den, u64 bit)
	{
		std::pair<double_int<bits, false>, double_int<bits, false>> res(0, 0);

		while (true)
		{
			res.second <<= 1;
			if (bit_test(num, bit)) ++ref_low64(res.second);

			if (res.second >= den)
			{
				res.second -= den;
				bit_set(res.first, bit);
			}

			if (bit-- == 0) break;
		}
		return res;
	}

	// utility function used internally - no divide by zero check
	template<u64 bits> inline constexpr std::pair<double_int<bits, false>, double_int<bits, false>> divmod_unchecked(const double_int<bits, false> &num, const double_int<bits, false> &den)
	{
		if constexpr (COMPACT(false) && bits >= DIVMOD_BLOCK_SKIP_THRESH)
		{
			// chop off any 64-bit blocks that are all zero
			u64 bit = bits - 1, *ptr = (u64*)&num;
			for (u64 i = bits / 64; i > 1; bit -= 64)
			{
				if (ptr[--i]) break; // if the previous byte is nonzero, we're done
			}

			return __divmod_unchecked(num, den, bit);
		}
		else
		{
			return __divmod_unchecked(num, den, bits - 1);
		}
	}

	inline constexpr std::pair<std::uint8_t, std::uint8_t> divmod(std::uint8_t num, std::uint8_t den) { return {num / den, num % den}; }
	inline constexpr std::pair<std::uint16_t, std::uint16_t> divmod(std::uint16_t num, std::uint16_t den) { return {num / den, num % den}; }
	inline constexpr std::pair<std::uint32_t, std::uint32_t> divmod(std::uint32_t num, std::uint32_t den) { return {num / den, num % den}; }
	inline constexpr std::pair<std::uint64_t, std::uint64_t> divmod(std::uint64_t num, std::uint64_t den) { return {num / den, num % den}; }
	template<u64 bits> inline constexpr std::pair<double_int<bits, false>, double_int<bits, false>> divmod(const double_int<bits, false> &num, const double_int<bits, false> &den)
	{
		if (!den) throw std::domain_error("divide by zero");
		return divmod_unchecked(num, den);
	}
	
	inline constexpr std::pair<std::int8_t, std::int8_t> divmod(std::int8_t num, std::int8_t den) { return {num / den, num % den}; }
	inline constexpr std::pair<std::int16_t, std::int16_t> divmod(std::int16_t num, std::int16_t den) { return {num / den, num % den}; }
	inline constexpr std::pair<std::int32_t, std::int32_t> divmod(std::int32_t num, std::int32_t den) { return {num / den, num % den}; }
	inline constexpr std::pair<std::int64_t, std::int64_t> divmod(std::int64_t num, std::int64_t den) { return {num / den, num % den}; }
	template<u64 bits> inline constexpr std::pair<double_int<bits, true>, double_int<bits, true>> divmod(const double_int<bits, true> &num, const double_int<bits, true> &den)
	{
		// we'll do signed divmod in terms of unsigned

		double_int<bits, false> unum = *(double_int<bits, false>*)&num;
		double_int<bits, false> uden = *(double_int<bits, false>*)&den;

		bool neg = false;
		if (bit_test(unum, bits - 1)) { make_neg(unum); neg = !neg; }
		if (bit_test(uden, bits - 1)) { make_neg(uden); neg = !neg; }

		auto res = divmod(unum, uden);
		
		if (neg)
		{
			make_neg(res.first);
			make_neg(res.second);
		}

		return { *(double_int<bits, true>*)&res.first, *(double_int<bits, true>*)&res.second };
	}

	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator/=(double_int<bits, sign> &num, const double_int<bits, sign> &den) { num = divmod(num, den).first; return num; }
	template<u64 bits, bool sign, typename T> inline constexpr double_int<bits, sign> &operator/=(double_int<bits, sign> &a, const T &b) noexcept { a /= (double_int<bits, sign>)b; return a; }
	
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator%=(double_int<bits, sign> &num, const double_int<bits, sign> &den) { num = divmod(num, den).second; return num; }
	template<u64 bits, bool sign, typename T> inline constexpr double_int<bits, sign> &operator%=(double_int<bits, sign> &a, const T &b) noexcept { a %= (double_int<bits, sign>)b; return a; }

	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> operator/(const double_int<bits, sign> &num, const double_int<bits, sign> &den) { return divmod(num, den).first; }
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> operator%(const double_int<bits, sign> &num, const double_int<bits, sign> &den) { return divmod(num, den).second; }

	SHORTERHAND_BINARY_FORMATTER(/)
	SHORTERHAND_BINARY_FORMATTER(%)

	// -- cmp -- //

	// cmp() will be equivalent to <=> but works for any version of C++
	
	// unsigned compare
	template<u64 bits> inline constexpr int cmp(const double_int<bits, false> &a, const double_int<bits, false> &b) noexcept
	{
		if constexpr (COMPACT(false) && bits >= CMP_LOOP_THRESH)
		{
			u64 *dest = (u64*)&a, *src = (u64*)&b;
			for (u64 i = bits / 64; i > 0; )
			{
				--i;
				if (dest[i] < src[i]) return -1;
				else if (dest[i] > src[i]) return 1;
			}
			return 0;
		}
		else
		{
			return a.high < b.high ? -1 : a.high > b.high ? 1 : a.low < b.low ? -1 : a.low > b.low ? 1 : 0;
		}
	}

	// signed compare
	template<u64 bits> inline constexpr int cmp(const double_int<bits, true> &a, const double_int<bits, true> &b) noexcept
	{
		// get unsigned comparison
		int ucmp = cmp(*(const double_int<bits, false>*)&a, *(const double_int<bits, false>*)&b);
		// do quick mafs to transform it into a signed comparison
		return is_neg(a) ^ is_neg(b) ? -ucmp : ucmp;
	}
	
	template<u64 bits1, u64 bits2, bool sign> inline constexpr bool operator==(const double_int<bits1, sign> &a, const double_int<bits2, sign> &b) noexcept { return cmp<std::max(bits1, bits2)>(a, b) == 0; }
	template<u64 bits1, u64 bits2, bool sign> inline constexpr bool operator!=(const double_int<bits1, sign> &a, const double_int<bits2, sign> &b) noexcept { return cmp<std::max(bits1, bits2)>(a, b) != 0; }
	template<u64 bits1, u64 bits2, bool sign> inline constexpr bool operator<(const double_int<bits1, sign> &a, const double_int<bits2, sign> &b) noexcept { return cmp<std::max(bits1, bits2)>(a, b) < 0; }
	template<u64 bits1, u64 bits2, bool sign> inline constexpr bool operator<=(const double_int<bits1, sign> &a, const double_int<bits2, sign> &b) noexcept { return cmp<std::max(bits1, bits2)>(a, b) <= 0; }
	template<u64 bits1, u64 bits2, bool sign> inline constexpr bool operator>(const double_int<bits1, sign> &a, const double_int<bits2, sign> &b) noexcept { return cmp<std::max(bits1, bits2)>(a, b) > 0; }
	template<u64 bits1, u64 bits2, bool sign> inline constexpr bool operator>=(const double_int<bits1, sign> &a, const double_int<bits2, sign> &b) noexcept { return cmp<std::max(bits1, bits2)>(a, b) >= 0; }

	// shorthand for comparing double_int with specified signage to a primitive type
	#define SHORTHAND_CMP_FORMATTER(sign, type) \
		template<u64 bits> inline constexpr bool operator==(const double_int<bits, sign> &a, type b) noexcept { return cmp(a, (double_int<bits, sign>)b) == 0; } \
		template<u64 bits> inline constexpr bool operator!=(const double_int<bits, sign> &a, type b) noexcept { return cmp(a, (double_int<bits, sign>)b) != 0; } \
		template<u64 bits> inline constexpr bool operator<(const double_int<bits, sign> &a, type b) noexcept { return cmp(a, (double_int<bits, sign>)b) < 0; } \
		template<u64 bits> inline constexpr bool operator<=(const double_int<bits, sign> &a, type b) noexcept { return cmp(a, (double_int<bits, sign>)b) <= 0; } \
		template<u64 bits> inline constexpr bool operator>(const double_int<bits, sign> &a, type b) noexcept { return cmp(a, (double_int<bits, sign>)b) > 0; } \
		template<u64 bits> inline constexpr bool operator>=(const double_int<bits, sign> &a, type b) noexcept { return cmp(a, (double_int<bits, sign>)b) >= 0; } \
		\
		template<u64 bits> inline constexpr bool operator==(type a, const double_int<bits, sign> &b) noexcept { return cmp((double_int<bits, sign>)a, b) == 0; } \
		template<u64 bits> inline constexpr bool operator!=(type a, const double_int<bits, sign> &b) noexcept { return cmp((double_int<bits, sign>)a, b) != 0; } \
		template<u64 bits> inline constexpr bool operator<(type a, const double_int<bits, sign> &b) noexcept { return cmp((double_int<bits, sign>)a, b) < 0; } \
		template<u64 bits> inline constexpr bool operator<=(type a, const double_int<bits, sign> &b) noexcept { return cmp((double_int<bits, sign>)a, b) <= 0; } \
		template<u64 bits> inline constexpr bool operator>(type a, const double_int<bits, sign> &b) noexcept { return cmp((double_int<bits, sign>)a, b) > 0; } \
		template<u64 bits> inline constexpr bool operator>=(type a, const double_int<bits, sign> &b) noexcept { return cmp((double_int<bits, sign>)a, b) >= 0; }

	SHORTHAND_CMP_FORMATTER(false, unsigned long long)
	SHORTHAND_CMP_FORMATTER(false, unsigned long)
	SHORTHAND_CMP_FORMATTER(false, unsigned int)
	SHORTHAND_CMP_FORMATTER(false, unsigned short)

	SHORTHAND_CMP_FORMATTER(true, signed long long)
	SHORTHAND_CMP_FORMATTER(true, signed long)
	SHORTHAND_CMP_FORMATTER(true, signed int)
	SHORTHAND_CMP_FORMATTER(true, signed short)
	
	// -- io -- //

	// given a hex character, converts it to an integer [0, 15] - returns true if it was a valid hex digit. ch is only meaningful on success.
	inline constexpr bool ext_hex(int &ch) noexcept
	{
		if (ch >= '0' && ch <= '9') { ch -= '0'; return true; }
		ch |= 32;
		if (ch >= 'a' && ch <= 'f') { ch = ch - 'a' + 10; return true; }
		return false;
	}

	template<u64 bits> std::ostream &operator<<(std::ostream &ostr, const double_int<bits, false> &val)
	{
		std::string str;
		int digit, dcount;
		u64 block;

		std::ostream::sentry sentry(ostr);
		if (!sentry) return ostr;

		// build the string
		if (ostr.flags() & std::ios::oct)
		{
			double_int<bits, false> cpy = val;

			while (true)
			{
				// get a block
				block = low64(cpy) & 0777777777777777777777;
				cpy >>= 63;
				dcount = 0;

				// write the block - do-while to ensure 0 gets printed
				do
				{
					str.push_back('0' + (block & 7));
					block >>= 3;
					++dcount;
				}
				while (block);

				// if there's still stuff, pad with zeroes and continue
				if (cpy) { for (; dcount < 21; ++dcount) str.push_back('0'); }
				// otherwise we're done
				else break;
			}
			if (ostr.flags() & std::ios::showbase) ostr.put('0');
		}
		else if(ostr.flags() & std::ios::hex)
		{
			double_int<bits, false> cpy = val;
			const char hex_alpha = ostr.flags() & std::ios::uppercase ? 'A' : 'a';

			while (true)
			{
				// get a block
				block = low64(cpy);
				cpy >>= 64;
				dcount = 0;

				// write the block - do-while to ensure 0 gets printed
				do
				{
					digit = block & 15;
					str.push_back(digit < 10 ? '0' + digit : hex_alpha + digit - 10);
					block >>= 4;
					++dcount;
				}
				while (block);

				// if there's still stuff, pad with zeroes and continue
				if (cpy) { for (; dcount < 16; ++dcount) str.push_back('0'); }
				// otherwise we're done
				else break;
			}
			if (ostr.flags() & std::ios::showbase)
			{
				ostr.put('0'); // uses put() to make sure we don't clobber ostr.width()
				ostr.put('x');
			}
		}
		else // default to dec mode
		{
			// divmod for double_int is really slow, so we'll extract several decimal digits in one go and then process them with built-in integers
			static constexpr double_int<bits, false> base = 10000000000000000000ul;
			
			std::pair<double_int<bits, false>, double_int<bits, false>> temp;
			temp.first = val; // temp.first takes the role of cpy in the other radicies

			while (true)
			{
				// get a block
				temp = divmod_unchecked(temp.first, base);
				block = low64(temp.second);
				dcount = 0;

				// write the block - do-while to ensure 0 gets printed
				do
				{
					str.push_back('0' + block % 10);
					block /= 10;
					++dcount;
				}
				while (block);

				// if there's still stuff, pad with zeroes and continue
				if (temp.first) { for (; dcount < 19; ++dcount) str.push_back('0'); }
				// otherwise we're done
				else break;
			}
		}

		// write the string
		if (ostr.flags() & std::ios::left)
		{
			for (std::size_t i = str.size(); i > 0;) ostr.put(str[--i]);
			for (int i = (int)ostr.width() - (int)str.size(); i > 0; --i) ostr.put(ostr.fill());
		}
		else // default to right
		{
			for (int i = (int)ostr.width() - (int)str.size(); i > 0; --i) ostr.put(ostr.fill());
			for (std::size_t i = str.size(); i > 0;) ostr.put(str[--i]);
		}

		// clobber width - makes sure the next write isn't weird
		ostr.width(0);
		return ostr;
	}
	template<u64 bits> std::ostream &operator<<(std::ostream &ostr, const double_int<bits, true> &val)
	{
		// we'll do signed io in terms of unsigned

		double_int<bits, true> cpy = val;

		// if it's negative make it positing and print the - sign
		if (bit_test(cpy, bits - 1))
		{
			make_neg(cpy);
			ostr.put('-'); // put the - sign and dec width by 1
			if (ostr.width() > 0) ostr.width(ostr.width() - 1);
		}
		// otherwise, is positive. print the + sign if showpos is set
		else if (ostr.flags() & std::ios::showpos)
		{
			ostr.put('+');// put the + sign and dec width by 1
			if (ostr.width() > 0) ostr.width(ostr.width() - 1);
		}

		// refer to cpy as unsigned
		ostr << *(double_int<bits, false>*)&cpy;

		return ostr;
	}
	
	template<u64 bits> std::istream &parse_udouble_int(std::istream &istr, double_int<bits, false> &val, bool noskipws = false)
	{
		val = 0; // start by zeroing value
		int digit, num_digits;
		u64 block;

		std::istream::sentry sentry(istr, noskipws);
		if (!sentry) return istr;

		// parse the string
		if (istr.flags() & std::ios::oct)
		{
			parse_oct:

			// first char must be an oct digit - don't extract
			if ((digit = istr.peek()) == EOF) { istr.setstate(std::ios::failbit | std::ios::eofbit); return istr; }
			if (digit < '0' || digit > '7') { istr.setstate(std::ios::failbit); return istr; }

			while (true)
			{
				block = 0; // read a block of 21 oct digits
				for (num_digits = 0; num_digits < 21; ++num_digits)
				{
					// get the digit and make sure it's in range - only extract it if it's good
					if ((digit = istr.peek()) == EOF) break;
					if (digit < '0' || digit > '7') break;
					istr.get();

					block <<= 3;
					block |= digit - '0';
				}
				// incorporate it into value
				val <<= (u64)(3 * num_digits);
				val |= block;

				if (num_digits < 21) break;
			}
		}
		else if (istr.flags() & std::ios::hex)
		{
			parse_hex:

			// first char must be an hex digit - don't extract
			if ((digit = istr.peek()) == EOF) { istr.setstate(std::ios::failbit | std::ios::eofbit); return istr; }
			if (!ext_hex(digit)) { istr.setstate(std::ios::failbit); return istr; }

			while (true)
			{
				block = 0; // read a block of 16 hex digits
				for (num_digits = 0; num_digits < 16; ++num_digits)
				{
					// get the digit and make sure it's in range - only extract it if it's good
					if ((digit = istr.peek()) == EOF) break;
					if (!ext_hex(digit)) break;
					istr.get();

					block <<= 4;
					block |= digit;
				}
				// incorporate it into value
				val <<= (u64)(4 * num_digits);
				val |= block;

				if (num_digits < 16) break;
			}
		}
		else if(istr.flags() & std::ios::dec)
		{
			parse_dec:

			// first char must be an dec digit - don't extract
			if ((digit = istr.peek()) == EOF) { istr.setstate(std::ios::failbit | std::ios::eofbit); return istr; }
			if (digit < '0' || digit > '9') { istr.setstate(std::ios::failbit); return istr; }

			while (true)
			{
				block = 0; // read a block of 19 dec digits
				for (num_digits = 0; num_digits < 19; ++num_digits)
				{
					// get the digit and make sure it's in range - only extract it if it's good
					if ((digit = istr.peek()) == EOF) break;
					if (digit < '0' || digit > '9') break;
					istr.get();

					block *= 10;
					block += digit - '0';
				}
				// incorporate it into value
				val *= (double_int<bits, false>)(u64)std::pow((u64)10, (u64)num_digits);
				val += (double_int<bits, false>)block;

				if (num_digits < 19) break;
			}
		}
		else // if no base specified, determine it from the suffix
		{
			// look at the top character - fail if none
			if ((digit = istr.peek()) == EOF) { istr.setstate(std::ios::failbit | std::ios::eofbit); return istr; }
			// if it's a zero we have a prefix
			if (digit == '0')
			{
				// get the next character - if there is none, it's a valid dec zero
				istr.get();
				if ((digit = istr.peek()) == EOF) return istr;

				// if second char is 'x', extract it and parse a hex value
				if (digit == 'x') { istr.get(); goto parse_hex; }
				// otherwise parse as oct
				else goto parse_oct;
			}
			// no prefix is dec
			else goto parse_dec;
		}

		return istr;
	}

	template<u64 bits> inline std::istream &operator>>(std::istream &istr, double_int<bits, false> &val) { parse_udouble_int(istr, val); return istr; }
	template<u64 bits> std::istream &operator>>(std::istream &istr, double_int<bits, true> &val)
	{
		// we'll do signed io in terms of unsigned

		int ch;
		bool neg = false;

		std::istream::sentry sentry(istr);
		if (!sentry) return istr;

		// look at the first char - fail if we can't get one
		if ((ch = istr.peek()) == EOF) { istr.setstate(std::ios::failbit | std::ios::eofbit); return istr; }

		// account for optional sign - extract it if present
		if (ch == '-') { istr.get(); neg = true; }
		else if (ch == '+') istr.get();

		// parse the value (and don't skip ws since we already did that)
		parse_udouble_int(istr, *(double_int<bits, false>*)&val, true);

		// account for sign in result
		if (neg) make_neg(val);

		return istr;
	}

	// -- extra utilities -- //

	// helper function for bool conversion (SFINAE wasn't working inside the class definition)
	template<u64 bits, bool sign, std::enable_if_t<(COMPACT(sign) && bits >= BOOL_LOOP_THRESH), int> = 0> inline constexpr bool to_bool(const double_int<bits, sign> &val) noexcept
	{
		u64 *ptr = (u64*)&val;
		for (u64 i = 0; i < bits / 64; ++i)
			if (ptr[i]) return true;
		return false;
	}
	template<u64 bits, bool sign, std::enable_if_t<(!COMPACT(sign) || bits < BOOL_LOOP_THRESH), int> = 0> inline constexpr bool to_bool(const double_int<bits, sign> &val) noexcept { return val.low || val.high; }

	// counts the number of set bits
	inline constexpr u64 num_set_bits(u64 val) noexcept
	{
		u64 num = 0;
		for (; val; ++num, val = val & (val - 1));
		return num;
	}
	template<u64 bits, bool sign> inline constexpr u64 num_set_bits(const double_int<bits, sign> &val) noexcept { return num_set_bits(val.low) + num_set_bits(val.high); }

	// gets the number of 64-bit blocks with significant bits (i.e. not including leading zero blocks)
	template<u64 bits, bool sign, std::enable_if_t<COMPACT(sign), int> = 0> inline constexpr u64 num_blocks(const double_int<bits, sign> &val) noexcept
	{
		// get pointer to u64 blocks (if compact, will be little-endian wrt 64-bit blocks)
		u64 *ptr = (u64*)&val;
		// go backwards through blocks (starting at high blocks) and stop on the first non-zero one
		for (int num = bits / 64; num > 0; --num)
			if (ptr[num - 1]) return num;
		// otherwise, value is bit-zero
		return 0;
	}
	template<u64 bits, bool sign, std::enable_if_t<!COMPACT(sign), int> = 0> inline constexpr u64 num_blocks(const double_int<bits, sign> &val) noexcept
	{
		// start at zero, return when we've shifted out all the significant blocks
		double_int<bits, sign> cpy = val;
		int num;
		for (num = 0; cpy; ++num, cpy >>= 64);
		return num;
	}

	// -- block access -- //

	// gets the low 64-bit word of a given value
	template<u64 bits, bool sign> inline constexpr u64 low64(const double_int<bits, sign> &val) noexcept
	{
		// because the low half always comes first, the entire recursive definition for any number of bits is little-endian wrt 64-bit blocks
		return *(u64*)&val;
	}
	// same as low64() but returns by reference
	template<u64 bits, bool sign> inline constexpr u64 &ref_low64(double_int<bits, sign> &val) noexcept
	{
		// because the low half always comes first, the entire recursive definition for any number of bits is little-endian wrt 64-bit blocks
		return *(u64*)&val;
	}

	// -- in-place unary -- //
	
	template<u64 bits, bool sign, std::enable_if_t<(COMPACT(sign) && bits >= NOT_LOOP_THRESH), int> = 0> inline constexpr void make_not(double_int<bits, sign> &val) noexcept
	{
		u64 *ptr = (u64*)&val;
		for (u64 i = 0; i < bits / 64; ++i) ptr[i] = ~ptr[i];
	}
	template<u64 bits, bool sign, std::enable_if_t<(!COMPACT(sign) || bits < NOT_LOOP_THRESH), int> = 0> inline constexpr void make_not(double_int<bits, sign> &val) noexcept
	{
		make_not(val.low);
		make_not(val.high);
	}
	template<bool sign> inline constexpr void make_not(double_int<128, sign> &val) noexcept { val.low = ~val.low; val.high = ~val.high; }

	template<u64 bits, bool sign> inline constexpr void make_neg(double_int<bits, sign> &val) noexcept { make_not(val); ++val; }

	// -- bit tests -- //

	inline constexpr bool bit_test(std::uint64_t val, u64 bit) noexcept { return (val >> bit) & 1; }
	template<u64 bits, bool sign, std::enable_if_t<COMPACT(sign), int> = 0> inline constexpr bool bit_test(const double_int<bits, sign> &val, u64 bit) noexcept
	{
		bit &= bits - 1;
		return bit_test(((u64*)&val)[bit >> 6], bit & 0x3f);
	}
	template<u64 bits, bool sign, std::enable_if_t<!COMPACT(sign), int> = 0> inline constexpr bool bit_test(const double_int<bits, sign> &val, u64 bit) noexcept
	{
		bit &= bits - 1;
		return bit >= bits / 2 ? bit_test(val.high, bit - bits / 2) : bit_test(val.low, bit);
	}
	
	inline constexpr void bit_set(std::uint64_t &val, u64 bit) noexcept { val |= (std::uint64_t)1 << bit; }
	template<u64 bits, bool sign, std::enable_if_t<COMPACT(sign), int> = 0> inline constexpr void bit_set(double_int<bits, sign> &val, u64 bit) noexcept
	{
		bit &= bits - 1;
		return bit_set(((u64*)&val)[bit >> 6], bit & 0x3f);
	}
	template<u64 bits, bool sign, std::enable_if_t<!COMPACT(sign), int> = 0> inline constexpr void bit_set(double_int<bits, sign> &val, u64 bit) noexcept
	{
		bit &= bits - 1;
		if (bit >= bits / 2) bit_set(val.high, bit - bits / 2); else bit_set(val.low, bit);
	}

	template<u64 bits, bool sign, std::enable_if_t<COMPACT(sign), int> = 0> inline constexpr bool is_neg(const double_int<bits, sign> &val) noexcept { return ((u64*)&val)[bits / 64 - 1] & 0x8000000000000000; }
	template<u64 bits, bool sign, std::enable_if_t<!COMPACT(sign), int> = 0> inline constexpr bool is_neg(const double_int<bits, sign> &val) noexcept { return bit_test(val, bits - 1); }
}

// -- std info definitions -- //

namespace std
{
	// -- double_int info -- //

	template<std::uint64_t bits, bool sign> struct is_integral<BiggerInts::double_int<bits, sign>> : std::true_type {};

	template<std::uint64_t bits, bool sign> struct is_signed<BiggerInts::double_int<bits, sign>> : std::integral_constant<bool, sign> {};
	template<std::uint64_t bits, bool sign> struct is_unsigned<BiggerInts::double_int<bits, sign>> : std::integral_constant<bool, !sign> {};

	template<std::uint64_t bits, bool sign> struct make_signed<BiggerInts::double_int<bits, sign>> { typedef BiggerInts::double_int<bits, true> type; };
	template<std::uint64_t bits, bool sign> struct make_unsigned<BiggerInts::double_int<bits, sign>> { typedef BiggerInts::double_int<bits, false> type; };
	
	template<std::uint64_t bits, bool sign> struct numeric_limits<BiggerInts::double_int<bits, sign>>
	{
		static constexpr bool is_specialized = true;
		static constexpr bool is_signed = sign;
		static constexpr bool is_integer = true;
		static constexpr bool is_exact = true;
		static constexpr bool has_infinity = false;
		static constexpr bool has_quiet_NaN = false;
		static constexpr bool has_signaling_NaN = false;
		static constexpr bool has_denorm = false;
		static constexpr bool has_denorm_loss = false;
		static constexpr std::float_round_style round_style = std::round_toward_zero;
		static constexpr bool is_iec559 = false;
		static constexpr bool is_bounded = true;
		static constexpr bool is_modulo = true;
		static constexpr int digits = sign ? bits - 1 : bits;
		static constexpr int digits10 = (int)(digits * 0.301029995663981195213738894724493026768189881462108541310); // log10(2)
		static constexpr int max_digits10 = 0;
		static constexpr int radix = 2;
		static constexpr int min_exponent = 0;
		static constexpr int min_exponent10 = 0;
		static constexpr int max_exponent = 0;
		static constexpr int max_exponent10 = 0;
		static constexpr bool traps = true;
		static constexpr bool tinyness_before = false;

		static constexpr BiggerInts::double_int<bits, sign> _min = sign ? (BiggerInts::double_int<bits, sign>)1 << (bits - 1) : (BiggerInts::double_int<bits, sign>)0;
		static constexpr BiggerInts::double_int<bits, sign> _max = sign ? ~((BiggerInts::double_int<bits, sign>)1 << (bits - 1)) : ~(BiggerInts::double_int<bits, sign>)0;
		static constexpr BiggerInts::double_int<bits, sign> _zero = (BiggerInts::double_int<bits, sign>)0;

		static constexpr const auto &min() { return _min; }
		static constexpr const auto &lowest() { return min(); }
		static constexpr const auto &max() {  return _max; }
		static constexpr const auto &epsilon() {  return _zero; }
		static constexpr const auto &round_error() { return epsilon(); }
		static constexpr const auto &infinity() { return epsilon(); }
		static constexpr const auto &quiet_NaN() { return epsilon(); }
		static constexpr const auto &signaling_NaN() { return epsilon(); }
		static constexpr const auto &denorm_min() { return epsilon(); }
	};

	// -- masked_single_int info -- //

	template<typename T, std::uint64_t bits> struct is_integral<BiggerInts::masked_single_int<T, bits>> : std::is_integral<T> {};

	template<typename T, std::uint64_t bits> struct is_signed<BiggerInts::masked_single_int<T, bits>> : std::is_signed<T> {};
	template<typename T, std::uint64_t bits> struct is_unsigned<BiggerInts::masked_single_int<T, bits>> : std::is_unsigned<T> {};

	template<typename T, std::uint64_t bits> struct make_signed<BiggerInts::masked_single_int<T, bits>> { typedef std::make_signed_t<T> type; };
	template<typename T, std::uint64_t bits> struct make_unsigned<BiggerInts::masked_single_int<T, bits>> { typedef std::make_unsigned_t<T> type; };

	template<typename T, std::uint64_t bits> struct numeric_limits<BiggerInts::masked_single_int<T, bits>>
	{
		static constexpr bool is_specialized = true;
		static constexpr bool is_signed = std::is_signed<T>::value;
		static constexpr bool is_integer = true;
		static constexpr bool is_exact = true;
		static constexpr bool has_infinity = false;
		static constexpr bool has_quiet_NaN = false;
		static constexpr bool has_signaling_NaN = false;
		static constexpr bool has_denorm = false;
		static constexpr bool has_denorm_loss = false;
		static constexpr std::float_round_style round_style = std::round_toward_zero;
		static constexpr bool is_iec559 = false;
		static constexpr bool is_bounded = true;
		static constexpr bool is_modulo = true;
		static constexpr int digits = std::is_signed<T>::value ? bits - 1 : bits;
		static constexpr int digits10 = (int)(digits * 0.301029995663981195213738894724493026768189881462108541310); // log10(2)
		static constexpr int max_digits10 = 0;
		static constexpr int radix = 2;
		static constexpr int min_exponent = 0;
		static constexpr int min_exponent10 = 0;
		static constexpr int max_exponent = 0;
		static constexpr int max_exponent10 = 0;
		static constexpr bool traps = true;
		static constexpr bool tinyness_before = false;

		static constexpr T _min = std::is_signed<T>::value ? (T)1 << (bits - 1) : 0;
		static constexpr T _max = std::is_signed<T>::value ? ~((T)1 << (bits - 1)) : ~(T)0;
		static constexpr T _zero = (T)0;

		static constexpr const auto &min() { return _min; }
		static constexpr const auto &lowest() { return min(); }
		static constexpr const auto &max() {  return _max; }
		static constexpr const auto &epsilon() { return _zero; }
		static constexpr const auto &round_error() { return epsilon(); }
		static constexpr const auto &infinity() { return epsilon(); }
		static constexpr const auto &quiet_NaN() { return epsilon(); }
		static constexpr const auto &signaling_NaN() { return epsilon(); }
		static constexpr const auto &denorm_min() { return epsilon(); }
	};
}

#endif
