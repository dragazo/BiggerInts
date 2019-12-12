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
#include <sstream>
#include <utility>
#include <exception>
#include <stdexcept>
#include <algorithm>
#include <string_view>
#include <cctype>

/*

Adds signed and unsigned types accessible by alias templates uint_t and int_t for any given number of bits.
These aliased types display behavior identical to the built-in integral types (moreover, uint_t<> and int_t<> for 8, 16, 32, or 64 bits returns a built-in type).
Any type returned has a full set of specialized type_traits and numeric_limits.
For a given number of bits, the signed/unsigned structures are identical (e.g. you can safely bind a reference from signed to unsigned with reinterpret_cast to avoid a (potentially expensive) copy).
Negative signed values are stored via 2's complement, just as with built-in integral types on anything made since the dark ages.

Report bugs to https://github.com/dragazo/BiggerInts/issues

*/

namespace BiggerInts
{
	namespace detail
	{
		template<typename T, std::uint64_t bits> struct masked_single_int;
		template<std::uint64_t bits, bool sign> struct double_int;

		// -- block access -- //

		// gets the low 64-bit word of a given value
		template<std::uint64_t bits, bool sign>
		constexpr std::uint64_t low64(const double_int<bits, sign> &val) noexcept { if constexpr (bits == 128) return val.low; else return low64(val.low); }
		// as low64() but returns by reference
		template<std::uint64_t bits, bool sign>
		constexpr std::uint64_t &ref_low64(double_int<bits, sign> &val) noexcept { if constexpr (bits == 128) return val.low; else return ref_low64(val.low); }

		// gets the high 64-bit word of a given value
		template<std::uint64_t bits, bool sign>
		constexpr std::uint64_t high64(const double_int<bits, sign> &val) noexcept { if constexpr (bits == 128) return val.high; else return high64(val.high); }

		// -- bit manip -- //

		// checks if a bit is set
		inline constexpr bool bit_test_unchecked(std::uint64_t val, std::uint64_t bit) noexcept { return (val >> bit) & 1; }
		template<std::uint64_t bits, bool sign>
		constexpr bool bit_test_unchecked(const double_int<bits, sign> &val, std::uint64_t bit) noexcept
		{
			return bit >= bits / 2 ? bit_test_unchecked(val.high, bit - bits / 2) : bit_test_unchecked(val.low, bit);
		}
		template<std::uint64_t bits, bool sign>
		constexpr bool bit_test(const double_int<bits, sign> &val, std::uint64_t bit) noexcept { return bit_test_unchecked(val, bit & (bits - 1)); }

		// sets the specified bit
		inline constexpr void bit_set_unchecked(std::uint64_t &val, std::uint64_t bit) noexcept { val |= (std::uint64_t)1 << bit; }
		template<std::uint64_t bits, bool sign>
		constexpr void bit_set_unchecked(double_int<bits, sign> &val, std::uint64_t bit) noexcept
		{
			if (bit >= bits / 2) bit_set_unchecked(val.high, bit - bits / 2); else bit_set_unchecked(val.low, bit);
		}
		template<std::uint64_t bits, bool sign>
		constexpr void bit_set(double_int<bits, sign> &val, std::uint64_t bit) noexcept { bit_set_unchecked(val, bit & (bits - 1)); }

		// determines if the value is negative under 2's complement (i.e. high bit is set)
		template<std::uint64_t bits, bool sign>
		constexpr bool is_neg(const double_int<bits, sign> &val) noexcept { return high64(val) & 0x8000000000000000; }

		// -- helpers -- //

		// returns true if val is a power of 2
		inline constexpr bool is_pow2(std::uint64_t val) noexcept { return val != 0 && (val & (val - 1)) == 0; }

		// given a size in bits, returns a power of 2 (also in bits) large enough to contain it and that is no smaller than 8
		inline constexpr std::uint64_t round_bits_up(std::uint64_t size)
		{
			if (size < 8) return 8;
			if (size > 0x8000000000000000ull) throw std::domain_error("bit size was too large");
			while (!is_pow2(size)) size = (size | (size - 1)) + 1;
			return size;
		}

		// -- #hax typedef omega synergy highjinks - don't even need to be defined, boi -- //

		template<std::uint64_t bits, bool sign, std::enable_if_t<(bits > 64 && is_pow2(bits)), int> = 0> double_int<bits, sign> returns_proper_type();
		template<std::uint64_t bits, bool sign, std::enable_if_t<(bits > 64 && !is_pow2(bits)), int> = 0> masked_single_int<double_int<round_bits_up(bits), sign>, bits> returns_proper_type();

		template<std::uint64_t bits, bool sign, std::enable_if_t<(bits == 64), int> = 0> std::conditional_t<sign, std::int64_t, std::uint64_t> returns_proper_type();
		template<std::uint64_t bits, bool sign, std::enable_if_t<(bits > 32 && bits < 64), int> = 0> masked_single_int<std::conditional_t<sign, std::int64_t, std::uint64_t>, bits> returns_proper_type();

		template<std::uint64_t bits, bool sign, std::enable_if_t<(bits == 32), int> = 0> std::conditional_t<sign, std::int32_t, std::uint32_t> returns_proper_type();
		template<std::uint64_t bits, bool sign, std::enable_if_t<(bits > 16 && bits < 32), int> = 0> masked_single_int<std::conditional_t<sign, std::int32_t, std::uint32_t>, bits> returns_proper_type();

		template<std::uint64_t bits, bool sign, std::enable_if_t<(bits == 16), int> = 0> std::conditional_t<sign, std::int16_t, std::uint16_t> returns_proper_type();
		template<std::uint64_t bits, bool sign, std::enable_if_t<(bits > 8 && bits < 16), int> = 0> masked_single_int<std::conditional_t<sign, std::int16_t, std::uint16_t>, bits> returns_proper_type();

		template<std::uint64_t bits, bool sign, std::enable_if_t<(bits == 8), int> = 0> std::conditional_t<sign, std::int8_t, std::uint8_t> returns_proper_type();
		template<std::uint64_t bits, bool sign, std::enable_if_t<(bits < 8), int> = 0> masked_single_int<std::conditional_t<sign, std::int8_t, std::uint8_t>, bits> returns_proper_type();

		// -- op utilities -- //

		// helper for divmod_unchecked - takes the bit to start with - don't use this directly
		template<std::uint64_t bits>
		constexpr std::pair<double_int<bits, false>, double_int<bits, false>> _divmod_unchecked(const double_int<bits, false> &num, const double_int<bits, false> &den, std::uint64_t bit)
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
		template<std::uint64_t bits>
		constexpr std::pair<double_int<bits, false>, double_int<bits, false>> divmod_unchecked(const double_int<bits, false> &num, const double_int<bits, false> &den)
		{
			return _divmod_unchecked(num, den, bits - 1);
		}

		// -- extra utilities -- //

		// helper function for bool conversion (SFINAE wasn't working inside the class definition)
		template<std::uint64_t bits, bool sign>
		constexpr bool to_bool(const double_int<bits, sign> &val) noexcept { return val.low || val.high; }

		// counts the number of set bits
		inline constexpr std::uint64_t num_set_bits(std::uint64_t val) noexcept
		{
			std::uint64_t num = 0;
			for (; val; ++num, val = val & (val - 1));
			return num;
		}

		template<std::uint64_t bits, bool sign>
		constexpr std::uint64_t num_set_bits(const double_int<bits, sign> &val) noexcept { return num_set_bits(val.low) + num_set_bits(val.high); }

		// gets the number of 64-bit blocks with significant bits (i.e. not including leading zero blocks)
		template<std::uint64_t bits, bool sign>
		constexpr std::uint64_t num_blocks(const double_int<bits, sign> &val) noexcept
		{
			// start at zero, return when we've shifted out all the significant blocks
			double_int<bits, sign> cpy = val;
			std::uint64_t num = 0;
			for (num = 0; cpy; ++num, cpy >>= 64);
			return num;
		}

		// -- in-place unary -- //

		// performs an in-place bitwise not
		template<std::uint64_t bits, bool sign>
		constexpr void make_not(double_int<bits, sign> &val) noexcept
		{
			if constexpr (bits == 128)
			{
				val.low = ~val.low;
				val.high = ~val.high;
			}
			else
			{
				make_not(val.low);
				make_not(val.high);
			}
		}

		// takes the 2's complement negative of the value in-place
		template<std::uint64_t bits, bool sign>
		constexpr void make_neg(double_int<bits, sign> &val) noexcept { make_not(val); ++val; }
	}

	template<std::uint64_t bits> using uint_t = decltype(detail::returns_proper_type<bits, false>());
	template<std::uint64_t bits> using int_t = decltype(detail::returns_proper_type<bits, true>());

	namespace detail
	{

		// contains an integral constant that represents the number of bits in a given (supported) (integral) type
		template<typename T> struct bit_count {};

		template<> struct bit_count<unsigned char> : std::integral_constant<std::uint64_t, sizeof(unsigned char) * CHAR_BIT> {};
		template<> struct bit_count<unsigned short> : std::integral_constant<std::uint64_t, sizeof(unsigned short) * CHAR_BIT> {};
		template<> struct bit_count<unsigned int> : std::integral_constant<std::uint64_t, sizeof(unsigned int) * CHAR_BIT> {};
		template<> struct bit_count<unsigned long> : std::integral_constant<std::uint64_t, sizeof(unsigned long) * CHAR_BIT> {};
		template<> struct bit_count<unsigned long long> : std::integral_constant<std::uint64_t, sizeof(unsigned long long) * CHAR_BIT> {};

		template<> struct bit_count<signed char> : std::integral_constant<std::uint64_t, sizeof(signed char) * CHAR_BIT> {};
		template<> struct bit_count<signed short> : std::integral_constant<std::uint64_t, sizeof(signed short) * CHAR_BIT> {};
		template<> struct bit_count<signed int> : std::integral_constant<std::uint64_t, sizeof(signed int) * CHAR_BIT> {};
		template<> struct bit_count<signed long> : std::integral_constant<std::uint64_t, sizeof(signed long) * CHAR_BIT> {};
		template<> struct bit_count<signed long long> : std::integral_constant<std::uint64_t, sizeof(signed long long) * CHAR_BIT> {};

		template<typename T, std::uint64_t bits> struct bit_count<masked_single_int<T, bits>> : std::integral_constant<std::uint64_t, bits> {};
		template<std::uint64_t bits, bool sign> struct bit_count<double_int<bits, sign>> : std::integral_constant<std::uint64_t, bits> {};

		// -- container impl -- //

		// holds an integral type T that, upon assignment, is masked to the specified number of bits. bits must be in the range (bit count T / 2, bit count T).
		template<typename T, std::uint64_t bits> struct masked_single_int
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
		template<std::uint64_t bits, bool sign> struct double_int
		{
		public: // -- data -- //

			static_assert(detail::is_pow2(bits) && bits >= 128, "double_int: bits was out of valid domain");

			typedef uint_t<bits / 2> half_t;

			half_t low;
			half_t high;

		public: // -- core -- //

			constexpr double_int() noexcept = default;

			constexpr double_int(const double_int&) noexcept = default;
			constexpr double_int &operator=(const double_int&) noexcept = default;

			constexpr explicit operator double_int<bits, !sign>() const noexcept { double_int<bits, !sign> cpy = *this; return cpy; }

		public: // -- promotion constructors -- //

			constexpr double_int(unsigned short val) noexcept : low(val), high(0u) {}
			constexpr double_int(unsigned int val) noexcept : low(val), high(0u) {}
			constexpr double_int(unsigned long val) noexcept : low(val), high(0u) {}
			constexpr double_int(unsigned long long val) noexcept : low(val), high(0u) {}

			constexpr double_int(signed short val) noexcept : low(val), high(val < 0 ? -1 : 0) {}
			constexpr double_int(signed int val) noexcept : low(val), high(val < 0 ? -1 : 0) {}
			constexpr double_int(signed long val) noexcept : low(val), high(val < 0 ? -1 : 0) {}
			constexpr double_int(signed long long val) noexcept : low(val), high(val < 0 ? -1 : 0) {}

			template<std::uint64_t _bits, std::enable_if_t<(bits > _bits), int> = 0>
				constexpr double_int(const double_int<_bits, false> &val) noexcept : low(val), high(0u) {}
				template<std::uint64_t _bits, std::enable_if_t<(bits > _bits), int> = 0>
				constexpr double_int(const double_int<_bits, true> &val) noexcept : low(val), high(is_neg(val) ? -1 : 0) {}

				template<bool other_sign>
				constexpr double_int(const double_int<bits, other_sign> &other) noexcept { low = other.low; high = other.high; }

		public: // -- promotion assignment -- //

			constexpr double_int &operator=(unsigned short val) noexcept { low = val; high = 0u; return *this; }
			constexpr double_int &operator=(unsigned int val) noexcept { low = val; high = 0u; return *this; }
			constexpr double_int &operator=(unsigned long val) noexcept { low = val; high = 0u; return *this; }
			constexpr double_int &operator=(unsigned long long val) noexcept { low = val; high = 0u; return *this; }

			constexpr double_int &operator=(signed short val) noexcept { low = val; high = val < 0 ? -1 : 0; return *this; }
			constexpr double_int &operator=(signed int val) noexcept { low = val; high = val < 0 ? -1 : 0; return *this; }
			constexpr double_int &operator=(signed long val) noexcept { low = val; high = val < 0 ? -1 : 0; return *this; }
			constexpr double_int &operator=(signed long long val) noexcept { low = val; high = val < 0 ? -1 : 0; return *this; }

			template<std::uint64_t _bits, std::enable_if_t<(bits > _bits), int> = 0>
				constexpr double_int &operator=(const double_int<_bits, false> &val) noexcept { low = val; high = 0u; return *this; }
				template<std::uint64_t _bits, std::enable_if_t<(bits > _bits), int> = 0>
				constexpr double_int &operator=(const double_int<_bits, true> &val) noexcept { low = val; high = is_neg(val) ? -1 : 0; return *this; }

				template<bool other_sign>
				constexpr double_int &operator=(const double_int<bits, other_sign> &other) noexcept { low = other.low; high = other.high; return *this; }

		public: // -- demotion conversion -- //

			constexpr explicit operator unsigned short() const noexcept { return (unsigned short)low; }
			constexpr explicit operator unsigned int() const noexcept { return (unsigned int)low; }
			constexpr explicit operator unsigned long() const noexcept { return (unsigned long)low; }
			constexpr explicit operator unsigned long long() const noexcept { return (unsigned long long)low; }

			constexpr explicit operator signed short() const noexcept { return (signed short)low; }
			constexpr explicit operator signed int() const noexcept { return (signed int)low; }
			constexpr explicit operator signed long() const noexcept { return (signed long)low; }
			constexpr explicit operator signed long long() const noexcept { return (signed long long)low; }

			template<std::uint64_t other_bits, bool other_sign, std::enable_if_t<(other_bits < bits), int> = 0>
				constexpr explicit operator double_int<other_bits, other_sign>() const noexcept
				{
					return (double_int<other_bits, other_sign>)low; // recursively select the low branch
				}

		public: // -- bool conversion -- //

			constexpr explicit operator bool() const noexcept { return detail::to_bool(*this); }
			constexpr friend bool operator!(const double_int &a) noexcept { return !(bool)a; }

		public: // -- static utilities -- //

			// attempts to parse the string into an integer value - returns true on successful parse.
			// base must be 10 (dec), 16 (hex), 8 (oct), or 0 to automatically determine base from C-style prefix in str.
			static bool try_parse(double_int &res, std::string_view str, int base = 10)
			{
				std::istringstream ss(std::string(str.begin(), str.end())); // unfortunately istringstream requires a copy because stdlib is dumb

				if (base == 10) {}
				else if (base == 16) ss.setf(std::ios::hex, std::ios::basefield);
				else if (base == 8) ss.setf(std::ios::oct, std::ios::basefield);
				else if (base == 0) ss.unsetf(std::ios::basefield);
				else throw std::invalid_argument("unrecognized base specified");

				ss >> res;
				if (!ss) return false; // parse needs to succeed
				while (std::isspace((unsigned char)ss.peek())) ss.get();
				if (!ss.eof()) return false; // we need to have parsed the entire string
				return true;
			}
			// as try_parse() but throws std::invalid_argument on failure
			static double_int parse(std::string_view str, int base = 10)
			{
				double_int res;
				if (!try_parse(res, str, base)) throw std::invalid_argument("failed to parse string");
				return res;
			}
		};



		// shorthand for binary ops +, -, *, etc. will refer to a form that uses both sides of double_int
		#define SHORTHAND_BINARY_FORMATTER(op, sign, type) \
		template<std::uint64_t bits> constexpr double_int<bits, sign> operator op(const double_int<bits, sign> &a, type b) noexcept { return a op (double_int<bits, sign>)b; } \
		template<std::uint64_t bits> constexpr double_int<bits, sign> operator op(type a, const double_int<bits, sign> &b) noexcept { return (double_int<bits, sign>)a op b; }
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

		template<std::uint64_t bits, bool sign>
		constexpr double_int<bits, sign> &operator++(double_int<bits, sign> &a) noexcept
		{
			if (!++a.low) ++a.high;
			return a;
		}
		template<std::uint64_t bits, bool sign>
		constexpr double_int<bits, sign> operator++(double_int<bits, sign> &a, int) noexcept { auto cpy = a; ++a; return cpy; }

		template<std::uint64_t bits, bool sign>
		constexpr double_int<bits, sign> &operator--(double_int<bits, sign> &a) noexcept
		{
			if (!a.low) --a.high;
			--a.low;
			return a;
		}
		template<std::uint64_t bits, bool sign>
		constexpr double_int<bits, sign> operator--(double_int<bits, sign> &a, int) noexcept { auto cpy = a; --a; return cpy; }

		// -- add -- //

		template<std::uint64_t bits, bool sign1, bool sign2> constexpr double_int<bits, sign1> &operator+=(double_int<bits, sign1> &a, const double_int<bits, sign2> &b) noexcept
		{
			a.low += b.low;
			if (a.low < b.low) ++a.high;
			a.high += b.high;
			return a;
		}

		template<std::uint64_t bits, bool sign, typename T>
		constexpr double_int<bits, sign> &operator+=(double_int<bits, sign> &a, const T &b) noexcept { a += (double_int<bits, sign>)b; return a; }

		template<std::uint64_t bits, bool sign>
		constexpr double_int<bits, sign> operator+(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { double_int<bits, sign> res = a; res += b; return res; }

		SHORTERHAND_BINARY_FORMATTER(+)

			// -- sub -- //

			template<std::uint64_t bits, bool sign1, bool sign2>
		constexpr double_int<bits, sign1> &operator-=(double_int<bits, sign1> &a, const double_int<bits, sign2> &b) noexcept
		{
			if (a.low < b.low) --a.high;
			a.low -= b.low;
			a.high -= b.high;
			return a;
		}

		template<std::uint64_t bits, bool sign, typename T>
		constexpr double_int<bits, sign> &operator-=(double_int<bits, sign> &a, const T &b) noexcept { a -= (double_int<bits, sign>)b; return a; }

		template<std::uint64_t bits, bool sign>
		constexpr double_int<bits, sign> operator-(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { double_int<bits, sign> res = a; res -= b; return res; }

		SHORTERHAND_BINARY_FORMATTER(-)

			// -- and -- //

			template<std::uint64_t bits, bool sign1, bool sign2>
		constexpr double_int<bits, sign1> &operator&=(double_int<bits, sign1> &a, const double_int<bits, sign2> &b) noexcept
		{
			a.low &= b.low;
			a.high &= b.high;
			return a;
		}

		template<std::uint64_t bits, bool sign, typename T>
		constexpr double_int<bits, sign> &operator&=(double_int<bits, sign> &a, const T &b) noexcept { a &= (double_int<bits, sign>)b; return a; }

		template<std::uint64_t bits, bool sign>
		constexpr double_int<bits, sign> operator&(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { double_int<bits, sign> res = a; res &= b; return res; }

		SHORTERHAND_BINARY_FORMATTER(&)

			// -- or -- //

			template<std::uint64_t bits, bool sign1, bool sign2>
		constexpr double_int<bits, sign1> &operator|=(double_int<bits, sign1> &a, const double_int<bits, sign2> &b) noexcept
		{
			a.low |= b.low;
			a.high |= b.high;
			return a;
		}

		template<std::uint64_t bits, bool sign, typename T>
		constexpr double_int<bits, sign> &operator|=(double_int<bits, sign> &a, const T &b) noexcept { a |= (double_int<bits, sign>)b; return a; }

		template<std::uint64_t bits, bool sign>
		constexpr double_int<bits, sign> operator|(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { double_int<bits, sign> res = a; res |= b; return res; }

		SHORTERHAND_BINARY_FORMATTER(| )

			// -- xor -- //

			template<std::uint64_t bits, bool sign1, bool sign2>
		constexpr double_int<bits, sign1> &operator^=(double_int<bits, sign1> &a, const double_int<bits, sign2> &b) noexcept
		{
			a.low ^= b.low;
			a.high ^= b.high;
			return a;
		}

		template<std::uint64_t bits, bool sign, typename T>
		constexpr double_int<bits, sign> &operator^=(double_int<bits, sign> &a, const T &b) noexcept { a ^= (double_int<bits, sign>)b; return a; }

		template<std::uint64_t bits, bool sign>
		constexpr double_int<bits, sign> operator^(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { auto res = a; res ^= b; return res; }

		SHORTERHAND_BINARY_FORMATTER(^)

			// -- shl/sal -- //

			template<std::uint64_t bits, bool sign>
		constexpr double_int<bits, sign> &operator<<=(double_int<bits, sign> &val, std::uint64_t count) noexcept
		{
			// count masked to bits limit - no-op on zero
			if (count &= bits - 1)
			{
				if (count >= bits / 2)
				{
					val.high = val.low;
					val.high <<= count - bits / 2;
					val.low = 0u;
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

		template<std::uint64_t bits, bool sign, typename T>
		constexpr double_int<bits, sign> operator<<(const double_int<bits, sign> &a, const T &count) noexcept { double_int<bits, sign> res = a; res <<= (std::uint64_t)count; return res; }

		// -- shr/sar -- //

		// shr
		template<std::uint64_t bits>
		constexpr double_int<bits, false> &operator>>=(double_int<bits, false> &val, std::uint64_t count) noexcept
		{
			// count masked to bits limit - no-op on zero
			if (count &= bits - 1)
			{
				if (count >= bits / 2)
				{
					val.low = val.high;
					val.low >>= count - bits / 2;
					val.high = 0u;
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
		template<std::uint64_t bits>
		constexpr double_int<bits, true> &operator>>=(double_int<bits, true> &val, std::uint64_t count) noexcept
		{
			// count masked to bits limit - no-op on zero
			if (count &= bits - 1)
			{
				if (count >= bits / 2)
				{
					val.low = val.high;
					val.low >>= count - bits / 2;
					if (detail::is_neg(val)) // fill with sign bit
					{
						val.high = -1;
						if (count != bits / 2) val.low |= val.high << (bits - count);
					}
					else val.high = 0;
				}
				else
				{
					val.low >>= count;
					val.low |= val.high << (bits / 2 - count);
					if (detail::is_neg(val)) // fill with sign bit
					{
						val.high >>= count;
						val.high |= (decltype(val.high))(-1) << (bits / 2 - count);
					}
					else val.high >>= count;
				}
			}
			return val;
		}

		template<std::uint64_t bits, bool sign, typename T>
		constexpr double_int<bits, sign> operator>>(const double_int<bits, sign> &a, const T &count) noexcept { double_int<bits, sign> res = a; res >>= (std::uint64_t)count; return res; }

		// -- unary -- //

		template<std::uint64_t bits, bool sign>
		constexpr double_int<bits, sign> operator+(const double_int<bits, sign> &a) noexcept { return a; }

		template<std::uint64_t bits>
		constexpr double_int<bits, true> operator-(const double_int<bits, true> &a) noexcept { double_int<bits, true> res = a; detail::make_neg(res); return res; }

		template<std::uint64_t bits, bool sign>
		constexpr double_int<bits, sign> operator~(const double_int<bits, sign> &a) noexcept { double_int<bits, sign> res = a; detail::make_not(res); return res; }

		// -- mul -- //

		template<std::uint64_t bits>
		constexpr double_int<bits, false> operator*(const double_int<bits, false> &a, const double_int<bits, false> &b) noexcept
		{
			double_int<bits, false> res = 0u, _a = a;
			for (std::uint64_t bit = 0; _a; ++bit, _a <<= 1)
				if (detail::bit_test(b, bit)) res += _a;
			return res;
		}
		template<std::uint64_t bits>
		constexpr double_int<bits, true> operator*(const double_int<bits, true> &a, const double_int<bits, true> &b) noexcept
		{
			// we'll do signed multiply in terms of unsigned
			double_int<bits, false> ua = a;
			double_int<bits, false> ub = b;

			bool neg = false;
			if (detail::bit_test(ua, bits - 1)) { detail::make_neg(ua); neg = !neg; }
			if (detail::bit_test(ub, bits - 1)) { detail::make_neg(ub); neg = !neg; }

			ua = ua * ub;

			if (neg) detail::make_neg(ua);

			return { ua };
		}

		SHORTERHAND_BINARY_FORMATTER(*)

			template<std::uint64_t bits, bool sign>
		constexpr double_int<bits, sign> &operator*=(double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { a = a * b; return a; }

		template<std::uint64_t bits, bool sign, typename T>
		constexpr double_int<bits, sign> &operator*=(double_int<bits, sign> &a, const T &b) noexcept { a *= (double_int<bits, sign>)b; return a; }

		inline constexpr std::pair<std::uint8_t, std::uint8_t> divmod(std::uint8_t num, std::uint8_t den) { return { num / den, num % den }; }
		inline constexpr std::pair<std::uint16_t, std::uint16_t> divmod(std::uint16_t num, std::uint16_t den) { return { num / den, num % den }; }
		inline constexpr std::pair<std::uint32_t, std::uint32_t> divmod(std::uint32_t num, std::uint32_t den) { return { num / den, num % den }; }
		inline constexpr std::pair<std::uint64_t, std::uint64_t> divmod(std::uint64_t num, std::uint64_t den) { return { num / den, num % den }; }

		template<std::uint64_t bits>
		constexpr std::pair<double_int<bits, false>, double_int<bits, false>> divmod(const double_int<bits, false> &num, const double_int<bits, false> &den)
		{
			if (!den) throw std::domain_error("divide by zero");
			return detail::divmod_unchecked(num, den);
		}

		inline constexpr std::pair<std::int8_t, std::int8_t> divmod(std::int8_t num, std::int8_t den) { return { num / den, num % den }; }
		inline constexpr std::pair<std::int16_t, std::int16_t> divmod(std::int16_t num, std::int16_t den) { return { num / den, num % den }; }
		inline constexpr std::pair<std::int32_t, std::int32_t> divmod(std::int32_t num, std::int32_t den) { return { num / den, num % den }; }
		inline constexpr std::pair<std::int64_t, std::int64_t> divmod(std::int64_t num, std::int64_t den) { return { num / den, num % den }; }

		template<std::uint64_t bits>
		constexpr std::pair<double_int<bits, true>, double_int<bits, true>> divmod(const double_int<bits, true> &num, const double_int<bits, true> &den)
		{
			// we'll do signed divmod in terms of unsigned
			double_int<bits, false> unum = num;
			double_int<bits, false> uden = den;

			bool neg = false;
			if (detail::bit_test(unum, bits - 1)) { detail::make_neg(unum); neg = !neg; }
			if (detail::bit_test(uden, bits - 1)) { detail::make_neg(uden); neg = !neg; }

			auto res = divmod(unum, uden);

			if (neg)
			{
				detail::make_neg(res.first);
				detail::make_neg(res.second);
			}

			return { res.first, res.second };
		}

		template<std::uint64_t bits, bool sign> inline constexpr double_int<bits, sign> &operator/=(double_int<bits, sign> &num, const double_int<bits, sign> &den) { num = divmod(num, den).first; return num; }
		template<std::uint64_t bits, bool sign, typename T> inline constexpr double_int<bits, sign> &operator/=(double_int<bits, sign> &a, const T &b) noexcept { a /= (double_int<bits, sign>)b; return a; }

		template<std::uint64_t bits, bool sign> inline constexpr double_int<bits, sign> &operator%=(double_int<bits, sign> &num, const double_int<bits, sign> &den) { num = divmod(num, den).second; return num; }
		template<std::uint64_t bits, bool sign, typename T> inline constexpr double_int<bits, sign> &operator%=(double_int<bits, sign> &a, const T &b) noexcept { a %= (double_int<bits, sign>)b; return a; }

		template<std::uint64_t bits, bool sign> inline constexpr double_int<bits, sign> operator/(const double_int<bits, sign> &num, const double_int<bits, sign> &den) { return divmod(num, den).first; }
		template<std::uint64_t bits, bool sign> inline constexpr double_int<bits, sign> operator%(const double_int<bits, sign> &num, const double_int<bits, sign> &den) { return divmod(num, den).second; }

		SHORTERHAND_BINARY_FORMATTER(/ )
			SHORTERHAND_BINARY_FORMATTER(%)

		// -- cmp -- //

		// cmp() will be equivalent to <=> but works for any version of C++

		// unsigned compare
		template<std::uint64_t bits>
		constexpr int cmp(const double_int<bits, false> &a, const double_int<bits, false> &b) noexcept
		{
			if constexpr (bits == 128)
			{
				return a.high < b.high ? -1 : a.high > b.high ? 1 : a.low < b.low ? -1 : a.low > b.low ? 1 : 0;
			}
			else
			{
				int temp = cmp(a.high, b.high);
				return temp != 0 ? temp : cmp(a.low, b.low);
			}
		}

		// signed compare
		template<std::uint64_t bits>
		constexpr int cmp(const double_int<bits, true> &a, const double_int<bits, true> &b) noexcept
		{
			// do the comparison in terms of unsigned
			int ucmp = cmp((double_int<bits, false>)a, (double_int<bits, false>)b);

			return detail::is_neg(a) ^ detail::is_neg(b) ? -ucmp : ucmp; // do quick mafs to transform it into a signed comparison
		}

		template<std::uint64_t bits1, std::uint64_t bits2, bool sign> constexpr bool operator==(const double_int<bits1, sign> &a, const double_int<bits2, sign> &b) noexcept { return cmp<std::max(bits1, bits2)>(a, b) == 0; }
		template<std::uint64_t bits1, std::uint64_t bits2, bool sign> constexpr bool operator!=(const double_int<bits1, sign> &a, const double_int<bits2, sign> &b) noexcept { return cmp<std::max(bits1, bits2)>(a, b) != 0; }
		template<std::uint64_t bits1, std::uint64_t bits2, bool sign> constexpr bool operator<(const double_int<bits1, sign> &a, const double_int<bits2, sign> &b) noexcept { return cmp<std::max(bits1, bits2)>(a, b) < 0; }
		template<std::uint64_t bits1, std::uint64_t bits2, bool sign> constexpr bool operator<=(const double_int<bits1, sign> &a, const double_int<bits2, sign> &b) noexcept { return cmp<std::max(bits1, bits2)>(a, b) <= 0; }
		template<std::uint64_t bits1, std::uint64_t bits2, bool sign> constexpr bool operator>(const double_int<bits1, sign> &a, const double_int<bits2, sign> &b) noexcept { return cmp<std::max(bits1, bits2)>(a, b) > 0; }
		template<std::uint64_t bits1, std::uint64_t bits2, bool sign> constexpr bool operator>=(const double_int<bits1, sign> &a, const double_int<bits2, sign> &b) noexcept { return cmp<std::max(bits1, bits2)>(a, b) >= 0; }

		// shorthand for comparing double_int with specified signage to a primitive type
		#define SHORTHAND_CMP_FORMATTER(sign, type) \
		template<std::uint64_t bits> constexpr bool operator==(const double_int<bits, sign> &a, type b) noexcept { return cmp(a, (double_int<bits, sign>)b) == 0; } \
		template<std::uint64_t bits> constexpr bool operator!=(const double_int<bits, sign> &a, type b) noexcept { return cmp(a, (double_int<bits, sign>)b) != 0; } \
		template<std::uint64_t bits> constexpr bool operator<(const double_int<bits, sign> &a, type b) noexcept { return cmp(a, (double_int<bits, sign>)b) < 0; } \
		template<std::uint64_t bits> constexpr bool operator<=(const double_int<bits, sign> &a, type b) noexcept { return cmp(a, (double_int<bits, sign>)b) <= 0; } \
		template<std::uint64_t bits> constexpr bool operator>(const double_int<bits, sign> &a, type b) noexcept { return cmp(a, (double_int<bits, sign>)b) > 0; } \
		template<std::uint64_t bits> constexpr bool operator>=(const double_int<bits, sign> &a, type b) noexcept { return cmp(a, (double_int<bits, sign>)b) >= 0; } \
		\
		template<std::uint64_t bits> constexpr bool operator==(type a, const double_int<bits, sign> &b) noexcept { return cmp((double_int<bits, sign>)a, b) == 0; } \
		template<std::uint64_t bits> constexpr bool operator!=(type a, const double_int<bits, sign> &b) noexcept { return cmp((double_int<bits, sign>)a, b) != 0; } \
		template<std::uint64_t bits> constexpr bool operator<(type a, const double_int<bits, sign> &b) noexcept { return cmp((double_int<bits, sign>)a, b) < 0; } \
		template<std::uint64_t bits> constexpr bool operator<=(type a, const double_int<bits, sign> &b) noexcept { return cmp((double_int<bits, sign>)a, b) <= 0; } \
		template<std::uint64_t bits> constexpr bool operator>(type a, const double_int<bits, sign> &b) noexcept { return cmp((double_int<bits, sign>)a, b) > 0; } \
		template<std::uint64_t bits> constexpr bool operator>=(type a, const double_int<bits, sign> &b) noexcept { return cmp((double_int<bits, sign>)a, b) >= 0; }

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

		template<std::uint64_t bits>
		std::ostream &operator<<(std::ostream &ostr, const double_int<bits, false> &val)
		{
			std::string str;
			int digit, dcount;
			std::uint64_t block;

			std::ostream::sentry sentry(ostr);
			if (!sentry) return ostr;

			// build the string
			if (ostr.flags() & std::ios::oct)
			{
				double_int<bits, false> cpy = val;

				while (true)
				{
					// get a block
					block = detail::low64(cpy) & 0777777777777777777777;
					cpy >>= 63;
					dcount = 0;

					// write the block - do-while to ensure 0 gets printed
					do
					{
						str.push_back('0' + (block & 7));
						block >>= 3;
						++dcount;
					} while (block);

					// if there's still stuff, pad with zeroes and continue
					if (cpy) { for (; dcount < 21; ++dcount) str.push_back('0'); }
					// otherwise we're done
					else break;
				}
				if (ostr.flags() & std::ios::showbase) ostr.put('0');
			}
			else if (ostr.flags() & std::ios::hex)
			{
				double_int<bits, false> cpy = val;
				const char hex_alpha = ostr.flags() & std::ios::uppercase ? 'A' : 'a';

				while (true)
				{
					// get a block
					block = detail::low64(cpy);
					cpy >>= 64;
					dcount = 0;

					// write the block - do-while to ensure 0 gets printed
					do
					{
						digit = block & 15;
						str.push_back(digit < 10 ? '0' + digit : hex_alpha + digit - 10);
						block >>= 4;
						++dcount;
					} while (block);

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
					temp = detail::divmod_unchecked(temp.first, base);
					block = detail::low64(temp.second);
					dcount = 0;

					// write the block - do-while to ensure 0 gets printed
					do
					{
						str.push_back('0' + block % 10);
						block /= 10;
						++dcount;
					} while (block);

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
		template<std::uint64_t bits>
		std::ostream &operator<<(std::ostream &ostr, const double_int<bits, true> &val)
		{
			// we'll do signed io in terms of unsigned
			double_int<bits, true> cpy = val;

			// if it's negative make it positing and print the - sign
			if (detail::bit_test(cpy, bits - 1))
			{
				detail::make_neg(cpy);
				ostr.put('-'); // put the - sign and dec width by 1
				if (ostr.width() > 0) ostr.width(ostr.width() - 1);
			}
			// otherwise, is positive. print the + sign if showpos is set
			else if (ostr.flags() & std::ios::showpos)
			{
				ostr.put('+'); // put the + sign and dec width by 1
				if (ostr.width() > 0) ostr.width(ostr.width() - 1);
			}

			// refer to cpy as unsigned
			ostr << (double_int<bits, false>)cpy;

			return ostr;
		}

		template<std::uint64_t bits>
		std::istream &parse_udouble_int(std::istream &istr, double_int<bits, false> &val, bool noskipws = false)
		{
			val = 0u; // start by zeroing value
			int digit, num_digits;
			std::uint64_t block;

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
					if (num_digits != 0)
					{
						// detect overflow
						if (detail::high64(val) >> (std::uint64_t)(64 - 3 * num_digits)) { istr.setstate(std::ios::failbit); return istr; }

						// incorporate it into value
						val <<= (std::uint64_t)(3 * num_digits);
						val |= block;
					}

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
					if (num_digits != 0)
					{
						// detect overflow
						if (detail::high64(val) >> (std::uint64_t)(64 - 4 * num_digits)) { istr.setstate(std::ios::failbit); return istr; }

						// incorporate it into value
						val <<= (std::uint64_t)(4 * num_digits);
						val |= block;
					}

					if (num_digits < 16) break;
				}
			}
			else if (istr.flags() & std::ios::dec)
			{
			parse_dec:

				double_int<bits * 2, false> bigval = 0u; // perform multiply/add in terms of a larger int type to easily tell if we overflow

				// first char must be an dec digit - don't extract
				if ((digit = istr.peek()) == EOF) { istr.setstate(std::ios::failbit | std::ios::eofbit); return istr; }
				if (digit < '0' || digit > '9') { istr.setstate(std::ios::failbit); return istr; }

				while (true)
				{
					block = 0;               // read a block of 19 dec digits
					std::uint64_t scale = 1; // scale factor to apply to val for block incorporation (see below)
					for (num_digits = 0; num_digits < 19; ++num_digits)
					{
						// get the digit and make sure it's in range - only extract it if it's good
						if ((digit = istr.peek()) == EOF) break;
						if (digit < '0' || digit > '9') break;
						istr.get();

						scale *= 10;
						block *= 10;
						block += digit - '0';
					}
					if (num_digits != 0)
					{
						// incorporate it into big value
						bigval *= scale;
						bigval += block;

						// detect overflow
						if (bigval.high) { istr.setstate(std::ios::failbit); return istr; }
					}

					if (num_digits < 19) break;
				}

				val = bigval.low; // copy low half of bigval to result
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

		template<std::uint64_t bits> inline std::istream &operator>>(std::istream &istr, double_int<bits, false> &val) { parse_udouble_int(istr, val); return istr; }
		template<std::uint64_t bits> std::istream &operator>>(std::istream &istr, double_int<bits, true> &val)
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
			double_int<bits, false> tmp;
			parse_udouble_int(istr, tmp, true);
			val = tmp; // store the tmp parse location back to val

			// account for sign in result
			if (neg) detail::make_neg(val);

			// check for signed overflow
			if (detail::is_neg(val) != neg) { istr.setstate(std::ios::failbit); return istr; }

			return istr;
		}
	}
}

// -- std info definitions -- //

namespace std
{
	// -- double_int info -- //

	template<std::uint64_t bits, bool sign> struct is_integral<BiggerInts::detail::double_int<bits, sign>> : std::true_type {};

	template<std::uint64_t bits, bool sign> struct is_signed<BiggerInts::detail::double_int<bits, sign>> : std::integral_constant<bool, sign> {};
	template<std::uint64_t bits, bool sign> struct is_unsigned<BiggerInts::detail::double_int<bits, sign>> : std::integral_constant<bool, !sign> {};

	template<std::uint64_t bits, bool sign> struct make_signed<BiggerInts::detail::double_int<bits, sign>> { typedef BiggerInts::detail::double_int<bits, true> type; };
	template<std::uint64_t bits, bool sign> struct make_unsigned<BiggerInts::detail::double_int<bits, sign>> { typedef BiggerInts::detail::double_int<bits, false> type; };
	
	template<std::uint64_t bits, bool sign> struct numeric_limits<BiggerInts::detail::double_int<bits, sign>>
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

		static constexpr BiggerInts::detail::double_int<bits, sign> _min = sign ? (BiggerInts::detail::double_int<bits, sign>)1 << (bits - 1) : (BiggerInts::detail::double_int<bits, sign>)0;
		static constexpr BiggerInts::detail::double_int<bits, sign> _max = sign ? ~((BiggerInts::detail::double_int<bits, sign>)1 << (bits - 1)) : ~(BiggerInts::detail::double_int<bits, sign>)0;
		static constexpr BiggerInts::detail::double_int<bits, sign> _zero = (BiggerInts::detail::double_int<bits, sign>)0;

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

	template<typename T, std::uint64_t bits> struct is_integral<BiggerInts::detail::masked_single_int<T, bits>> : std::is_integral<T> {};

	template<typename T, std::uint64_t bits> struct is_signed<BiggerInts::detail::masked_single_int<T, bits>> : std::is_signed<T> {};
	template<typename T, std::uint64_t bits> struct is_unsigned<BiggerInts::detail::masked_single_int<T, bits>> : std::is_unsigned<T> {};

	template<typename T, std::uint64_t bits> struct make_signed<BiggerInts::detail::masked_single_int<T, bits>> { typedef std::make_signed_t<T> type; };
	template<typename T, std::uint64_t bits> struct make_unsigned<BiggerInts::detail::masked_single_int<T, bits>> { typedef std::make_unsigned_t<T> type; };

	template<typename T, std::uint64_t bits> struct numeric_limits<BiggerInts::detail::masked_single_int<T, bits>>
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
