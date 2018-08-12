#ifndef BIGGER_INTS_H
#define BIGGER_INTS_H

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <type_traits>
#include <iostream>
#include <string>
#include <utility>
#include <exception>
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

// if true, uses shortcut: address of struct == address of first member
#define USE_C_PTR_TO_FIRST 1
// if true, uses compactness shortcuts in double_int functions
#define USE_COMPACTNESS 1

// the types to use for zero/sign extension. zero ext must be unsigned. sign ext must be signed.
#define ZEXT_TYPE unsigned
#define SEXT_TYPE signed

// ADC is a faster addition algorithm, but has a higher constant cost. sizes below this will use the (worse) algorithm. at and above will use ADC.
#define ADC_THRESH 1024

// --------------------------- //

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
		if (size > 0x8000000000000000) throw std::domain_error("bit size was too large");
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

	template<> struct bit_count<std::uint8_t> : std::integral_constant<u64, 8> {};
	template<> struct bit_count<std::uint16_t> : std::integral_constant<u64, 16> {};
	template<> struct bit_count<std::uint32_t> : std::integral_constant<u64, 32> {};
	template<> struct bit_count<std::uint64_t> : std::integral_constant<u64, 64> {};

	template<> struct bit_count<std::int8_t> : std::integral_constant<u64, 8> {};
	template<> struct bit_count<std::int16_t> : std::integral_constant<u64, 16> {};
	template<> struct bit_count<std::int32_t> : std::integral_constant<u64, 32> {};
	template<> struct bit_count<std::int64_t> : std::integral_constant<u64, 64> {};

	template<typename T, u64 bits> struct bit_count<masked_single_int<T, bits>> : std::integral_constant<u64, bits> {};
	template<u64 bits, bool sign> struct bit_count<double_int<bits, sign>> : std::integral_constant<u64, bits> {};
	
	// -- container iml -- //

	// holds an integral type T that, upon assignment, is masked to the specified number of bits. bits must be in the range (bit count T / 2, bit count T).
	template<typename T, u64 bits> struct masked_single_int
	{
	private:

		static_assert(bits > bit_count<T>::value / 2 && bits < bit_count<T>::value, "masked_single_int: bit count out of range");

		static constexpr T mask = ((T)1 << bits) - 1;
		T val;

	public:
		inline constexpr operator T() const { return val; }

		inline constexpr masked_single_int(T _val = 0) noexcept : val(_val & mask) {}
		inline constexpr masked_single_int &operator=(T _val) { val = _val & mask; return *this; }

		inline constexpr masked_single_int &operator++() { val = (val + 1) & mask; return *this; }
		inline constexpr masked_single_int &operator--() { val = (val - 1) & mask; return *this; }

		inline constexpr T operator++(int) { T ret = val; val = (val + 1) & mask; return ret; }
		inline constexpr T operator--(int) { T ret = val; val = (val - 1) & mask; return ret; }

		inline constexpr masked_single_int &operator+=(T _val) { val = (val + _val) & mask; return *this; }
		inline constexpr masked_single_int &operator-=(T _val) { val = (val + _val) & mask; return *this; }

		inline constexpr masked_single_int &operator*=(T _val) { val = (val * _val) & mask; return *this; }
		inline constexpr masked_single_int &operator/=(T _val) { val = (val / _val) & mask; return *this; }
		inline constexpr masked_single_int &operator%=(T _val) { val = (val % _val) & mask; return *this; }

		inline constexpr masked_single_int &operator&=(T _val) { val = (val & _val) & mask; return *this; }
		inline constexpr masked_single_int &operator|=(T _val) { val = (val | _val) & mask; return *this; }
		inline constexpr masked_single_int &operator^=(T _val) { val = (val ^ _val) & mask; return *this; }

		inline constexpr masked_single_int &operator<<=(T _val) { val = (val << _val) & mask; return *this; }
		inline constexpr masked_single_int &operator>>=(T _val) { val = (val >> _val) & mask; return *this; }
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

		inline constexpr double_int(std::uint8_t val) noexcept : high((ZEXT_TYPE)0), low(val) {}
		inline constexpr double_int(std::uint16_t val) noexcept : high((ZEXT_TYPE)0), low(val) {}
		inline constexpr double_int(std::uint32_t val) noexcept : high((ZEXT_TYPE)0), low(val) {}
		inline constexpr double_int(std::uint64_t val) noexcept : high((ZEXT_TYPE)0), low(val) {}
		
		inline constexpr double_int(std::int8_t val) noexcept : high(val < 0 ? (SEXT_TYPE)-1 : (SEXT_TYPE)0), low(val) {}
		inline constexpr double_int(std::int16_t val) noexcept : high(val < 0 ? (SEXT_TYPE)-1 : (SEXT_TYPE)0), low(val) {}
		inline constexpr double_int(std::int32_t val) noexcept : high(val < 0 ? (SEXT_TYPE)-1 : (SEXT_TYPE)0), low(val) {}
		inline constexpr double_int(std::int64_t val) noexcept : high(val < 0 ? (SEXT_TYPE)-1 : (SEXT_TYPE)0), low(val) {}

		template<u64 _bits, std::enable_if_t<(bits > _bits), int> = 0> inline constexpr double_int(const double_int<_bits, false> &val) noexcept : high((ZEXT_TYPE)0), low(val) {}
		template<u64 _bits, std::enable_if_t<(bits > _bits), int> = 0> inline constexpr double_int(const double_int<_bits, true> &val) noexcept : high(is_neg(val) ? (SEXT_TYPE)-1 : (SEXT_TYPE)0), low(val) {}

		template<u64 _bits, bool _sign, std::enable_if_t<(bits == _bits), int> = 0> inline constexpr double_int(const double_int<_bits, _sign> &val) noexcept { *this = *(double_int*)&val; }

	public: // -- promotion assignment -- //

		inline constexpr double_int &operator=(std::uint8_t val) noexcept { high = (ZEXT_TYPE)0; low = val; return *this; }
		inline constexpr double_int &operator=(std::uint16_t val) noexcept { high = (ZEXT_TYPE)0; low = val; return *this; }
		inline constexpr double_int &operator=(std::uint32_t val) noexcept { high = (ZEXT_TYPE)0; low = val; return *this; }
		inline constexpr double_int &operator=(std::uint64_t val) noexcept { high = (ZEXT_TYPE)0; low = val; return *this; }

		inline constexpr double_int &operator=(std::int8_t val) noexcept { high = val < 0 ? (SEXT_TYPE)-1 : (SEXT_TYPE)0; low = val; return *this; }
		inline constexpr double_int &operator=(std::int16_t val) noexcept { high = val < 0 ? (SEXT_TYPE)-1 : (SEXT_TYPE)0; low = val; return *this; }
		inline constexpr double_int &operator=(std::int32_t val) noexcept { high = val < 0 ? (SEXT_TYPE)-1 : (SEXT_TYPE)0; low = val; return *this; }
		inline constexpr double_int &operator=(std::int64_t val) noexcept { high = val < 0 ? (SEXT_TYPE)-1 : (SEXT_TYPE)0; low = val; return *this; }
		
		template<u64 _bits, std::enable_if_t<(bits > _bits), int> = 0> inline constexpr double_int &operator=(const double_int<_bits, false> &val) noexcept { high = (ZEXT_TYPE)0; low = val; return *this; }
		template<u64 _bits, std::enable_if_t<(bits > _bits), int> = 0> inline constexpr double_int &operator=(const double_int<_bits, true> &val) noexcept { high = is_neg(val) ? (SEXT_TYPE)-1 : (SEXT_TYPE)0; low = val; return *this; }

		template<u64 _bits, bool _sign, std::enable_if_t<(bits == _bits), int> = 0> inline constexpr double_int &operator=(const double_int<_bits, _sign> &val) noexcept { *this = *(double_int*)&val; return *this; }

	public: // -- demotion conversion -- //

		inline constexpr explicit operator std::uint8_t() const noexcept { return (std::uint8_t)low64(*this); }
		inline constexpr explicit operator std::uint16_t() const noexcept { return (std::uint16_t)low64(*this); }
		inline constexpr explicit operator std::uint32_t() const noexcept { return (std::uint32_t)low64(*this); }
		inline constexpr explicit operator std::uint64_t() const noexcept { return (std::uint64_t)low64(*this); }

		inline constexpr explicit operator std::int8_t() const noexcept { return (std::int8_t)low64(*this); }
		inline constexpr explicit operator std::int16_t() const noexcept { return (std::int16_t)low64(*this); }
		inline constexpr explicit operator std::int32_t() const noexcept { return (std::int32_t)low64(*this); }
		inline constexpr explicit operator std::int64_t() const noexcept { return (std::int64_t)low64(*this); }

		template<u64 _bits, bool _sign, std::enable_if_t<(_bits < bits), int> = 0>
		inline constexpr explicit operator double_int<_bits, _sign>() const noexcept
		{
			#if USE_C_PTR_TO_FIRST
			return *(double_int<_bits, _sign>*)this;
			#else
			return (double_int<_bits, _sign>)low;
			#endif
		}

	public: // -- bool conversion -- //

		inline constexpr explicit operator bool() const noexcept { return low || high; }
		inline constexpr friend bool operator!(const double_int &a) noexcept { return !(bool)a; }
	};

	// shorthand for binary ops +, -, *, etc. will refer to a form that uses both sides of double_int
	#define SHORTHAND_BINARY_FORMATTER(op, sign, type) \
		template<u64 bits> inline constexpr double_int<bits, sign> operator op(const double_int<bits, sign> &a, type b) noexcept { return a op (double_int<bits, sign>)b; } \
		template<u64 bits> inline constexpr double_int<bits, sign> operator op(type a, const double_int<bits, sign> &b) noexcept { return (double_int<bits, sign>)a op b; }
	// does all the standard shorthands for a given op
	#define SHORTERHAND_BINARY_FORMATTER(op) \
		SHORTHAND_BINARY_FORMATTER(op, false, std::uint64_t) \
		SHORTHAND_BINARY_FORMATTER(op, false, std::uint32_t) \
		SHORTHAND_BINARY_FORMATTER(op, false, std::uint16_t) \
		SHORTHAND_BINARY_FORMATTER(op, false, std::uint8_t) \
		\
		SHORTHAND_BINARY_FORMATTER(op, true, std::int64_t) \
		SHORTHAND_BINARY_FORMATTER(op, true, std::int32_t) \
		SHORTHAND_BINARY_FORMATTER(op, true, std::int16_t) \
		SHORTHAND_BINARY_FORMATTER(op, true, std::int8_t)

	// -- inc/dec -- //

	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator++(double_int<bits, sign> &a) noexcept
	{
		#if USE_C_PTR_TO_FIRST && USE_COMPACTNESS

		// we need to make sure it's actually compact
		static_assert(sizeof(double_int<bits, sign>) == bits / CHAR_BIT, "type was not compact");

		// wind up 64-bit blocks - break on the first that doesn't overflow
		u64 *ptr = (u64*)&a;
		for (u64 i = 0; i < bits / 64; ++i)
			if (++ptr[i]) break;
		return a;

		#else

		if (!++a.low) ++a.high;
		return a;

		#endif
	}
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator--(double_int<bits, sign> &a) noexcept
	{
		#if USE_C_PTR_TO_FIRST && USE_COMPACTNESS

		// we need to make sure it's actually compact
		static_assert(sizeof(double_int<bits, sign>) == bits / CHAR_BIT, "type was not compact");

		// wind up 64-bit blocks - break on the first that doesn't underflow
		u64 *ptr = (u64*)&a;
		for (u64 i = 0; i < bits / 64; ++i)
			if (ptr[i]--) break;
		return a;

		#else

		if (!a.low) --a.high; --a.low;
		return a;

		#endif
	}

	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> operator++(double_int<bits, sign> &a, int) noexcept { double_int<bits, sign> val = a; ++a; return val; }
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> operator--(double_int<bits, sign> &a, int) noexcept { double_int<bits, sign> val = a; --a; return val; }

	// -- add -- //

	template<u64 bits, bool sign1, bool sign2, std::enable_if_t<(bits < ADC_THRESH), int> = 0> inline constexpr double_int<bits, sign1> &operator+=(double_int<bits, sign1> &a, const double_int<bits, sign2> &b) noexcept
	{
		a.low += b.low;
		if (a.low < b.low) ++a.high;
		a.high += b.high;
		return a;
	}
	template<u64 bits, bool sign1, bool sign2, std::enable_if_t<(bits >= ADC_THRESH), int> = 0> inline constexpr double_int<bits, sign1> &operator+=(double_int<bits, sign1> &a, const double_int<bits, sign2> &b) noexcept
	{
		ADC(a, *(double_int<bits, sign1>*)&b, false);
		return a;
	}

	template<u64 bits, bool sign, typename T> inline constexpr double_int<bits, sign> &operator+=(double_int<bits, sign> &a, const T &b) noexcept { a += (double_int<bits, sign>)b; return a; }

	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> operator+(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { double_int<bits, sign> res = a; res += b; return res; }

	SHORTERHAND_BINARY_FORMATTER(+)

	// -- sub -- //

	template<u64 bits, bool sign1, bool sign2, std::enable_if_t<(bits < ADC_THRESH), int> = 0> inline constexpr double_int<bits, sign1> &operator-=(double_int<bits, sign1> &a, const double_int<bits, sign2> &b) noexcept
	{
		if (a.low < b.low) --a.high;
		a.low -= b.low;
		a.high -= b.high;
		return a;
	}
	template<u64 bits, bool sign1, bool sign2, std::enable_if_t<(bits >= ADC_THRESH), int> = 0> inline constexpr double_int<bits, sign1> &operator-=(double_int<bits, sign1> &a, const double_int<bits, sign2> &b) noexcept
	{
		ADCN(a, *(double_int<bits, sign1>*)&b, true);
		return a;
	}

	template<u64 bits, bool sign, typename T> inline constexpr double_int<bits, sign> &operator-=(double_int<bits, sign> &a, const T &b) noexcept { a -= (double_int<bits, sign>)b; return a; }

	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> operator-(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { double_int<bits, sign> res = a; res -= b; return res; }

	SHORTERHAND_BINARY_FORMATTER(-)

	// -- and -- //

	template<u64 bits, bool sign1, bool sign2> inline constexpr double_int<bits, sign1> &operator&=(double_int<bits, sign1> &a, const double_int<bits, sign2> &b) noexcept { a.low &= b.low; a.high &= b.high; return a; }
	template<u64 bits, bool sign, typename T> inline constexpr double_int<bits, sign> &operator&=(double_int<bits, sign> &a, const T &b) noexcept { a &= (double_int<bits, sign>)b; return a; }
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator&=(double_int<bits, sign> &a, std::uint64_t b) noexcept { a = low64(a) & b; return a; }
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator&=(double_int<bits, sign> &a, std::uint32_t b) noexcept { a = low64(a) & b; return a; }
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator&=(double_int<bits, sign> &a, std::uint16_t b) noexcept { a = low64(a) & b; return a; }
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator&=(double_int<bits, sign> &a, std::uint8_t b) noexcept { a = low64(a) & b; return a; }
	
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> operator&(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { double_int<bits, sign> res = a; res &= b; return res; }

	SHORTERHAND_BINARY_FORMATTER(&)

	// -- or -- //

	template<u64 bits, bool sign1, bool sign2> inline constexpr double_int<bits, sign1> &operator|=(double_int<bits, sign1> &a, const double_int<bits, sign2> &b) noexcept { a.low |= b.low; a.high |= b.high; return a; }
	template<u64 bits, bool sign, typename T> inline constexpr double_int<bits, sign> &operator|=(double_int<bits, sign> &a, const T &b) noexcept { a |= (double_int<bits, sign>)b; return a; }
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator|=(double_int<bits, sign> &a, std::uint64_t b) noexcept { ref_low64(a) |= b; return a; }
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator|=(double_int<bits, sign> &a, std::uint32_t b) noexcept { ref_low64(a) |= b; return a; }
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator|=(double_int<bits, sign> &a, std::uint16_t b) noexcept { ref_low64(a) |= b; return a; }
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator|=(double_int<bits, sign> &a, std::uint8_t b) noexcept { ref_low64(a) |= b; return a; }

	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> operator|(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { double_int<bits, sign> res = a; res |= b; return res; }

	SHORTERHAND_BINARY_FORMATTER(|)

	// -- xor -- //

	template<u64 bits, bool sign1, bool sign2> inline constexpr double_int<bits, sign1> &operator^=(double_int<bits, sign1> &a, const double_int<bits, sign2> &b) noexcept { a.low ^= b.low; a.high ^= b.high; return a; }
	template<u64 bits, bool sign, typename T> inline constexpr double_int<bits, sign> &operator^=(double_int<bits, sign> &a, const T &b) noexcept { a ^= (double_int<bits, sign>)b; return a; }
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator^=(double_int<bits, sign> &a, std::uint64_t b) noexcept { ref_low64(a) ^= b; return a; }
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator^=(double_int<bits, sign> &a, std::uint32_t b) noexcept { ref_low64(a) ^= b; return a; }
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator^=(double_int<bits, sign> &a, std::uint16_t b) noexcept { ref_low64(a) ^= b; return a; }
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator^=(double_int<bits, sign> &a, std::uint8_t b) noexcept { ref_low64(a) ^= b; return a; }

	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> operator^(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { double_int<bits, sign> res = a; res ^= b; return res; }

	SHORTERHAND_BINARY_FORMATTER(^)

	// -- shl/sal -- //

	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator<<=(double_int<bits, sign> &val, u64 count) noexcept
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

	template<u64 bits> inline constexpr double_int<bits, false> &operator>>=(double_int<bits, false> &val, u64 count) noexcept
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
	template<u64 bits> inline constexpr double_int<bits, true> &operator>>=(double_int<bits, true> &val, u64 count) noexcept
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

	template<u64 bits> inline constexpr double_int<bits, true> operator-(const double_int<bits, true> &a) noexcept { double_int<bits, true> res = ~a; ++res; return res; }

	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> operator~(const double_int<bits, sign> &a) noexcept { double_int<bits, sign> res; res.high = ~a.high; res.low = ~a.low; return res; }

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

	// utility function used internally - no divide by zero check
	template<u64 bits> inline constexpr std::pair<double_int<bits, false>, double_int<bits, false>> divmod_unchecked(const double_int<bits, false> &num, const double_int<bits, false> &den)
	{
		std::pair<double_int<bits, false>, double_int<bits, false>> res(0, 0);
		u64 bit = bits - 1;
		while (true)
		{
			res.second <<= 1;
			if (bit_test(num, bit)) bit_set(res.second, 0);

			if (res.second >= den)
			{
				res.second -= den;
				bit_set(res.first, bit);
			}

			if (bit-- == 0) break;
		}
		return res;
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

	template<u64 bits> inline constexpr int cmp(const double_int<bits, false> &a, const double_int<bits, false> &b) noexcept { return a.high < b.high ? -1 : a.high > b.high ? 1 : a.low < b.low ? -1 : a.low > b.low ? 1 : 0; }
	template<u64 bits> inline constexpr int cmp(const double_int<bits, true> &a, const double_int<bits, true> &b) noexcept { int ucmp = a.high < b.high ? -1 : a.high > b.high ? 1 : a.low < b.low ? -1 : a.low > b.low ? 1 : 0; return bit_test(a, bits - 1) ^ bit_test(b, bits - 1) ? -ucmp : ucmp; }
	
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

	SHORTHAND_CMP_FORMATTER(false, std::uint64_t)
	SHORTHAND_CMP_FORMATTER(false, std::uint32_t)
	SHORTHAND_CMP_FORMATTER(false, std::uint16_t)
	SHORTHAND_CMP_FORMATTER(false, std::uint8_t)

	SHORTHAND_CMP_FORMATTER(true, std::int64_t)
	SHORTHAND_CMP_FORMATTER(true, std::int32_t)
	SHORTHAND_CMP_FORMATTER(true, std::int16_t)
	SHORTHAND_CMP_FORMATTER(true, std::int8_t)

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
		double_int<bits, false> cpy = val;
		std::string str;
		int digit, dcount;
		u64 block;

		std::ostream::sentry sentry(ostr);
		if (!sentry) return ostr;

		// build the string
		if (ostr.flags() & std::ios::oct)
		{
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
			static constexpr double_int<bits, false> base = 10000000000000000000;
			std::pair<double_int<bits, false>, double_int<bits, false>> temp;
			while (true)
			{
				// get a block
				temp = divmod_unchecked(cpy, base);
				block = low64(temp.second);
				cpy = temp.first;
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
				if (cpy) { for (; dcount < 19; ++dcount) str.push_back('0'); }
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

	// adds b to a and optionally increments if carry is true. returns carry out
	inline constexpr bool ADC(u64 &a, u64 b, bool carry) noexcept
	{
		// do the addition and get carry out
		a += b;
		bool _c = a < b;

		// account for carry in
		if (carry && !++a) _c = true;

		// return carry out
		return _c;
	}
	template<u64 bits, bool sign> inline constexpr bool ADC(double_int<bits, sign> &a, const double_int<bits, sign> &b, bool carry) noexcept { return ADC(a.high, b.high, ADC(a.low, b.low, carry)); }

	// as ADC but adds the bitwise not of b
	inline constexpr bool ADCN(u64 &a, u64 b, bool carry) noexcept
	{
		// do the addition and get carry out
		a += ~b;
		bool _c = a < ~b;

		// account for carry in
		if (carry && !++a) _c = true;

		// return carry out
		return _c;
	}
	template<u64 bits, bool sign> inline constexpr bool ADCN(double_int<bits, sign> &a, const double_int<bits, sign> &b, bool carry) noexcept { return ADCN(a.high, b.high, ADCN(a.low, b.low, carry)); }

	// counts the number of set bits
	inline constexpr u64 num_set_bits(u64 val) noexcept
	{
		u64 num = 0;
		for (; val; ++num, val = val & (val - 1));
		return num;
	}
	template<u64 bits, bool sign> inline constexpr u64 num_set_bits(const double_int<bits, sign> &val) noexcept { return num_set_bits(val.low) + num_set_bits(val.high); }

	// gets the number of 64-bit blocks with significant bits (i.e. not including leading zero blocks)
	template<u64 bits, bool sign> inline constexpr u64 num_blocks(const double_int<bits, sign> &val) noexcept
	{
		#if USE_C_PTR_TO_FIRST && USE_COMPACTNESS
		// we need to make sure it's actually compact
		static_assert(sizeof(double_int<bits, sign>) == bits / CHAR_BIT, "type was not compact");

		// get pointer to u64 blocks (if compact, will be little-endian wrt 64-bit blocks)
		u64 *ptr = (u64*)&val;
		// go backwards through blocks (starting at high blocks) and stop on the first non-zero one
		for (int num = bits / 64; num > 0; --num)
			if (ptr[num - 1]) return num;
		// otherwise, value is bit-zero
		return 0;

		#else

		// start at zero, return when we've shifted out all the significant blocks
		double_int<bits, sign> cpy = val;
		int num;
		for (num = 0; cpy; ++num, cpy >>= 64);
		return num;

		#endif
	}

	// gets the low 64-bit word of a given value
	template<bool sign> u64 low64(const double_int<128, sign> &val) noexcept { return val.low; }
	template<u64 bits, bool sign> u64 low64(const double_int<bits, sign> &val) noexcept
	{
		#if USE_C_PTR_TO_FIRST
		// because the low half always comes first, the entire recursive definition for any number of bits is little-endian wrt 64-bit blocks
		return *(u64*)&val;
		#else
		return low64(val.low);
		#endif
	}

	// same as low64() but returns by reference
	template<bool sign> inline constexpr u64 &ref_low64(double_int<128, sign> &val) noexcept { return val.low; }
	template<u64 bits, bool sign> inline constexpr u64 &ref_low64(double_int<bits, sign> &val) noexcept
	{
		#if USE_C_PTR_TO_FIRST
		// because the low half always comes first, the entire recursive definition for any number of bits is little-endian wrt 64-bit blocks
		return *(u64*)&val;
		#else
		return ref_low64(val.low);
		#endif
	}

	template<bool sign> inline constexpr void make_not(double_int<128, sign> &val) noexcept { val.low = ~val.low; val.high = ~val.high; }
	template<u64 bits, bool sign> inline constexpr void make_not(double_int<bits, sign> &val) noexcept { make_not(val.low); make_not(val.high); }

	template<u64 bits, bool sign> inline constexpr void make_neg(double_int<bits, sign> &val) noexcept { make_not(val); ++val; }

	inline constexpr bool bit_test(std::uint64_t val, u64 bit) noexcept { return (val >> bit) & 1; }
	template<u64 bits, bool sign> inline constexpr bool bit_test(const double_int<bits, sign> &val, u64 bit) noexcept { bit &= bits - 1; return bit >= bits / 2 ? bit_test(val.high, bit - bits / 2) : bit_test(val.low, bit); }

	inline constexpr void bit_set(std::uint64_t &val, u64 bit) noexcept { val |= (std::uint64_t)1 << bit; }
	template<u64 bits, bool sign> inline constexpr void bit_set(double_int<bits, sign> &val, u64 bit) noexcept { bit &= bits - 1; if (bit >= bits / 2) bit_set(val.high, bit - bits / 2); else bit_set(val.low, bit); }

	template<u64 bits, bool sign> inline constexpr bool is_neg(const double_int<bits, sign> &val) noexcept { return bit_test(val, bits - 1); }
}

// -- std info definitions -- //

namespace std
{
	// -- double_int info -- //

	template<std::uint64_t bits, bool sign> struct std::is_integral<BiggerInts::double_int<bits, sign>> : std::true_type {};

	template<std::uint64_t bits, bool sign> struct std::is_signed<BiggerInts::double_int<bits, sign>> : std::bool_constant<sign> {};
	template<std::uint64_t bits, bool sign> struct std::is_unsigned<BiggerInts::double_int<bits, sign>> : std::bool_constant<!sign> {};

	template<std::uint64_t bits, bool sign> struct std::make_signed<BiggerInts::double_int<bits, sign>> { typedef BiggerInts::double_int<bits, true> type; };
	template<std::uint64_t bits, bool sign> struct std::make_unsigned<BiggerInts::double_int<bits, sign>> { typedef BiggerInts::double_int<bits, false> type; };
	
	template<std::uint64_t bits, bool sign> struct std::numeric_limits<BiggerInts::double_int<bits, sign>>
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

		static constexpr const BiggerInts::double_int<bits, sign> &min() { static const BiggerInts::double_int<bits, sign> val = sign ? (BiggerInts::double_int<bits, sign>)1 << (bits - 1) : 0; return val; }
		static constexpr const BiggerInts::double_int<bits, sign> &lowest() { return min(); }
		static constexpr const BiggerInts::double_int<bits, sign> &max() { static const BiggerInts::double_int<bits, sign> val = sign ? ~((BiggerInts::double_int<bits, sign>)1 << (bits - 1)) : ~(BiggerInts::double_int<bits, sign>)0; return val; }
		static constexpr const BiggerInts::double_int<bits, sign> &epsilon() { static const BiggerInts::double_int<bits, sign> val = 0; return val; }
		static constexpr const BiggerInts::double_int<bits, sign> &round_error() { return epsilon(); }
		static constexpr const BiggerInts::double_int<bits, sign> &infinity() { return epsilon(); }
		static constexpr const BiggerInts::double_int<bits, sign> &quiet_NaN() { return epsilon(); }
		static constexpr const BiggerInts::double_int<bits, sign> &signaling_NaN() { return epsilon(); }
		static constexpr const BiggerInts::double_int<bits, sign> &denorm_min() { return epsilon(); }
	};

	// -- masked_single_int info -- //

	template<typename T, std::uint64_t bits> struct std::is_integral<BiggerInts::masked_single_int<T, bits>> : std::is_integral<T> {};

	template<typename T, std::uint64_t bits> struct std::is_signed<BiggerInts::masked_single_int<T, bits>> : std::is_signed<T> {};
	template<typename T, std::uint64_t bits> struct std::is_unsigned<BiggerInts::masked_single_int<T, bits>> : std::is_unsigned<T> {};

	template<typename T, std::uint64_t bits> struct std::make_signed<BiggerInts::masked_single_int<T, bits>> { typedef std::make_signed_t<T> type; };
	template<typename T, std::uint64_t bits> struct std::make_unsigned<BiggerInts::masked_single_int<T, bits>> { typedef std::make_unsigned_t<T> type; };

	template<typename T, std::uint64_t bits> struct std::numeric_limits<BiggerInts::masked_single_int<T, bits>>
	{
		static constexpr bool is_specialized = true;
		static constexpr bool is_signed = std::is_signed_v<T>;
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
		static constexpr int digits = std::is_signed_v<T> ? bits - 1 : bits;
		static constexpr int digits10 = (int)(digits * 0.301029995663981195213738894724493026768189881462108541310); // log10(2)
		static constexpr int max_digits10 = 0;
		static constexpr int radix = 2;
		static constexpr int min_exponent = 0;
		static constexpr int min_exponent10 = 0;
		static constexpr int max_exponent = 0;
		static constexpr int max_exponent10 = 0;
		static constexpr bool traps = true;
		static constexpr bool tinyness_before = false;

		static constexpr const BiggerInts::masked_single_int<T, bits> &min() { static const BiggerInts::masked_single_int<T, bits> val = std::is_signed_v<T> ? (T)1 << (bits - 1) : 0; return val; }
		static constexpr const BiggerInts::masked_single_int<T, bits> &lowest() { return min(); }
		static constexpr const BiggerInts::masked_single_int<T, bits> &max() { static const BiggerInts::masked_single_int<T, bits> val = std::is_signed_v<T> ? ~((T)1 << (bits - 1)) : ~(T)0; return val; }
		static constexpr const BiggerInts::masked_single_int<T, bits> &epsilon() { static const BiggerInts::masked_single_int<T, bits> val = 0; return val; }
		static constexpr const BiggerInts::masked_single_int<T, bits> &round_error() { return epsilon(); }
		static constexpr const BiggerInts::masked_single_int<T, bits> &infinity() { return epsilon(); }
		static constexpr const BiggerInts::masked_single_int<T, bits> &quiet_NaN() { return epsilon(); }
		static constexpr const BiggerInts::masked_single_int<T, bits> &signaling_NaN() { return epsilon(); }
		static constexpr const BiggerInts::masked_single_int<T, bits> &denorm_min() { return epsilon(); }
	};
}

#endif
