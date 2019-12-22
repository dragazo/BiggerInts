#ifndef DRAGAZO_BIGGER_INTS_H
#define DRAGAZO_BIGGER_INTS_H

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
#include <vector>

// msvc
#ifdef _MSC_VER
#include <intrin.h>
#endif

/*

adds larger integer types:

the template aliases uint_t and int_t represent fixed-size integers of any number of bits which is a power of 2 (e.g. 32-bit, 64-bit, 128-bit, etc.).
the signed form int_t is stored in 2's complement.
these are non-dynamic types and all their operations are constexpr.
this also makes e.g. 128-bit arithmetic very fast due to loop unrolling optimizations by the compiler.

the type alias bigint represents an arbitrary-sized signed integer in 2's complement.

report any bugs to https://github.com/dragazo/BiggerInts/issues

*/

namespace BiggerInts
{
	namespace detail
	{
		// -- compilation switches -- //

		// divmod() performs lots of bit shifting. many of these operations (on average) can be grouped together to run faster.
		// however this _may_ have some overhead that would actually make it slower for smaller integer types.
		// bit counts at or above this value will use the algorithm that minimizes the total number of shifts during divmod().
		static inline constexpr std::uint64_t divmod_shift_count_thresh = 256;
		// same as above but for multiplication.
		static inline constexpr std::uint64_t multiply_shift_count_thresh = 512;

		// --------------------------------
		
		template<typename T, std::uint64_t bits> struct masked_single_int;
		template<std::uint64_t bits, bool sign> struct double_int;
		class bigint;

		// -- bit manip -- //

		// checks if a bit is set
		template<std::uint64_t bits, bool sign>
		constexpr bool bit_test_unchecked(const double_int<bits, sign> &val, std::uint64_t bit) noexcept
		{
			return (val.blocks[bit / 64] >> (bit % 64)) & 1;
		}
		template<std::uint64_t bits, bool sign>
		constexpr bool bit_test(const double_int<bits, sign> &val, std::uint64_t bit) noexcept { return bit_test_unchecked(val, bit & (bits - 1)); }

		// sets the specified bit
		template<std::uint64_t bits, bool sign>
		constexpr void bit_set_unchecked(double_int<bits, sign> &val, std::uint64_t bit) noexcept { val.blocks[bit / 64] |= (std::uint64_t)1 << (bit % 64); }
		template<std::uint64_t bits, bool sign>
		constexpr void bit_set(double_int<bits, sign> &val, std::uint64_t bit) noexcept { bit_set_unchecked(val, bit & (bits - 1)); }

		// determines if the value is negative under 2's complement (i.e. high bit is set)
		template<std::uint64_t bits, bool sign>
		constexpr bool is_neg(const double_int<bits, sign> &val) noexcept { return val.blocks[bits / 64 - 1] & 0x8000000000000000; }

		// determines the index of the highest set bit
		inline constexpr std::size_t highest_set_bit(std::uint64_t val)
		{
			std::size_t i = 32, res = 0;
			for (; i > 0; i >>= 1)
			{
				std::uint64_t v = val >> i;
				if (v) { val = v; res += i; }
			}
			return res;
		}
		template<std::uint64_t bits, bool sign>
		constexpr std::size_t highest_set_bit(const double_int<bits, sign> &val) noexcept
		{
			std::size_t i = bits / 64;
			while (i-- > 0) if (val.blocks[i]) return i * 64 + highest_set_bit(val.blocks[i]);
			return 0;
		}

		// computes the full (unsigned) product of a * b (first = low, second = high)
		constexpr std::pair<std::uint64_t, std::uint64_t> _mul_u64(std::uint64_t a, std::uint64_t b) noexcept
		{
			std::uint64_t a_lo = (std::uint32_t)a;
			std::uint64_t b_lo = (std::uint32_t)b;
			std::uint64_t a_hi = a >> 32;
			std::uint64_t b_hi = b >> 32;

			std::uint64_t lo = a_lo * b_lo;
			std::uint64_t hi = a_hi * b_hi;

			std::uint64_t mid = a_hi * b_lo;
			if (std::uint64_t temp = a_lo * b_hi; (mid += temp) < temp) hi += 0x100000000ull;

			hi += mid >> 32;
			if (std::uint64_t temp = mid << 32; (lo += temp) < temp) ++hi;

			return { lo, hi };
		}
		// as _mul_u64() except is allowed to use faster, platform-specific intrinsics that might not be constexpr
		inline std::pair<std::uint64_t, std::uint64_t> _mul_u64_fast(std::uint64_t a, std::uint64_t b) noexcept
		{
			#ifdef _MSC_VER
			a = _umul128(a, b, &b);
			return { a, b };
			#else
			return _mul_u64(a, b);
			#endif
		}

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
		constexpr std::pair<double_int<bits, false>, double_int<bits, false>> divmod_unchecked(const double_int<bits, false> &num, const double_int<bits, false> &den)
		{
			std::pair<double_int<bits, false>, double_int<bits, false>> res(0u, 0u);
			
			std::uint64_t bit = highest_set_bit(num);

			// select algorithm at compile time
			if constexpr (bits >= divmod_shift_count_thresh)
			{
				const std::size_t den_highest_bit = highest_set_bit(den);
				std::size_t res_second_highest_bit = 0;
				std::size_t shift_count = 0;

				while (true)
				{
					++shift_count;
					if (bit_test_unchecked(num, bit))
					{
						res.second <<= shift_count;
						res_second_highest_bit += shift_count;
						shift_count = 0;
						++res.second.blocks[0];
					}
					else if (res_second_highest_bit + shift_count >= den_highest_bit)
					{
						res.second <<= shift_count;
						res_second_highest_bit += shift_count;
						shift_count = 0;
					}

					if (res.second >= den)
					{
						res.second -= den;
						res_second_highest_bit = highest_set_bit(res.second);
						bit_set_unchecked(res.first, bit);
					}

					if (bit-- == 0) break;
				}

				res.second <<= shift_count;
			}
			else
			{
				while (true)
				{
					res.second <<= 1;
					if (bit_test_unchecked(num, bit)) ++res.second.blocks[0];

					if (res.second >= den)
					{
						res.second -= den;
						bit_set_unchecked(res.first, bit);
					}

					if (bit-- == 0) break;
				}
			}

			return res;
		}

		// -- in-place unary -- //

		// performs an in-place bitwise not
		template<std::uint64_t bits, bool sign>
		constexpr void make_not(double_int<bits, sign> &val) noexcept
		{
			for (std::size_t i = 0; i < bits / 64; ++i) val.blocks[i] = ~val.blocks[i];
		}

		// takes the 2's complement negative of the value in-place
		template<std::uint64_t bits, bool sign>
		constexpr void make_neg(double_int<bits, sign> &val) noexcept { make_not(val); ++val; }

		// checks if T is a built-in integer type (signed or unsigned)
		template<typename T> struct is_builtin_int : std::false_type {};
		template<> struct is_builtin_int<unsigned short> : std::true_type {};
		template<> struct is_builtin_int<unsigned int> : std::true_type {};
		template<> struct is_builtin_int<unsigned long> : std::true_type {};
		template<> struct is_builtin_int<unsigned long long> : std::true_type {};
		template<> struct is_builtin_int<short> : std::true_type {};
		template<> struct is_builtin_int<int> : std::true_type {};
		template<> struct is_builtin_int<long> : std::true_type {};
		template<> struct is_builtin_int<long long> : std::true_type {};

		// checks if T is a biggerints type (signed or unsigned)
		template<typename T> struct is_biggerints_type : std::false_type {};
		template<std::uint64_t bits, bool sign> struct is_biggerints_type<double_int<bits, sign>> : std::true_type {};
		template<> struct is_biggerints_type<bigint> : std::true_type {};

		// checks if T is a signed double_int
		template<typename T> struct is_signed_double_int : std::false_type {};
		template<std::uint64_t bits> struct is_signed_double_int<double_int<bits, true>> : std::true_type {};
	}

	template<std::uint64_t bits> using uint_t = decltype(detail::returns_proper_type<bits, false>());
	template<std::uint64_t bits> using int_t = decltype(detail::returns_proper_type<bits, true>());
	using bigint = detail::bigint;

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

		// we'll just elect to define bigint as having zero bits - it really has infinite, but that's not representable
		template<> struct bit_count<bigint> : std::integral_constant<std::uint64_t, 0ull> {};

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

			std::uint64_t blocks[bits / 64];

		public: // -- core -- //

			constexpr double_int() noexcept : blocks{} {};

			constexpr double_int(const double_int&) noexcept = default;
			constexpr double_int &operator=(const double_int&) noexcept = default;

			constexpr double_int(const double_int<bits, !sign> &other) noexcept
			{
				for (std::size_t i = 0; i < bits / 64; ++i) blocks[i] = other.blocks[i];
			}
			constexpr double_int &operator=(const double_int<bits, !sign> &other) noexcept
			{
				for (std::size_t i = 0; i < bits / 64; ++i) blocks[i] = other.blocks[i];
				return *this;
			}

		public: // -- promotion constructors -- //

			constexpr double_int(unsigned short val) noexcept : double_int((unsigned long long)val) {}
			constexpr double_int(unsigned int val) noexcept : double_int((unsigned long long)val) {}
			constexpr double_int(unsigned long val) noexcept : double_int((unsigned long long)val) {}
			constexpr double_int(unsigned long long val) noexcept : blocks{}
			{
				blocks[0] = val;
				for (std::size_t i = 1; i < bits / 64; ++i) blocks[i] = 0;
			}

			constexpr double_int(signed short val) noexcept : double_int((signed long long)val) {}
			constexpr double_int(signed int val) noexcept : double_int((signed long long)val) {}
			constexpr double_int(signed long val) noexcept : double_int((signed long long)val) {}
			constexpr double_int(signed long long val) noexcept : blocks{}
			{
				blocks[0] = val;
				std::uint64_t fill = val < 0 ? -1 : 0;
				for (std::size_t i = 1; i < bits / 64; ++i) blocks[i] = fill;
			}

			template<std::uint64_t _bits, bool _sign, std::enable_if_t<(bits > _bits), int> = 0>
			constexpr double_int(const double_int<_bits, _sign> &other) noexcept : blocks{}
			{
				for (std::size_t i = 0; i < _bits / 64; ++i) blocks[i] = other.blocks[i];
				std::uint64_t fill = 0;
				if constexpr (_sign) fill = detail::is_neg(other) ? -1 : 0;
				for (std::size_t i = _bits / 64; i < bits / 64; ++i) blocks[i] = fill;
			}

		public: // -- promotion assignment -- //

			constexpr double_int &operator=(unsigned short val) noexcept { *this = (unsigned long long)val; return *this; }
			constexpr double_int &operator=(unsigned int val) noexcept { *this = (unsigned long long)val; return *this; }
			constexpr double_int &operator=(unsigned long val) noexcept { *this = (unsigned long long)val; return *this; }
			constexpr double_int &operator=(unsigned long long val) noexcept
			{
				blocks[0] = val;
				for (std::size_t i = 1; i < bits / 64; ++i) blocks[i] = 0;
				return *this;
			}

			constexpr double_int &operator=(signed short val) noexcept { *this = (signed long long)val; return *this; }
			constexpr double_int &operator=(signed int val) noexcept { *this = (signed long long)val; return *this; }
			constexpr double_int &operator=(signed long val) noexcept { *this = (signed long long)val; return *this; }
			constexpr double_int &operator=(signed long long val) noexcept
			{
				blocks[0] = val;
				std::uint64_t fill = val < 0 ? -1 : 0;
				for (std::size_t i = 1; i < bits / 64; ++i) blocks[i] = fill;
				return *this;
			}

			template<std::uint64_t _bits, bool _sign, std::enable_if_t<(bits > _bits), int> = 0>
					constexpr double_int &operator=(const double_int<_bits, _sign> &other) noexcept
					{
						for (std::size_t i = 0; i < _bits / 64; ++i) blocks[i] = other.blocks[i];
						std::uint64_t fill;
						if constexpr (_sign) fill = detail::is_neg(other) ? -1 : 0; else fill = 0;
						for (std::size_t i = _bits / 64; i < bits / 64; ++i) blocks[i] = fill;
						return *this;
					}

		public: // -- double_int demotions -- //

			template<std::uint64_t _bits, bool _sign, std::enable_if_t<(bits < _bits), int> = 0>
				constexpr explicit double_int(const double_int<_bits, _sign> &other) noexcept : blocks{}
				{
					for (std::size_t i = 0; i < bits / 64; ++i) blocks[i] = other.blocks[i];
				}

		public: // -- demotion conversion -- //

			constexpr explicit operator unsigned short() const noexcept { return (unsigned short)blocks[0]; }
			constexpr explicit operator unsigned int() const noexcept { return (unsigned int)blocks[0]; }
			constexpr explicit operator unsigned long() const noexcept { return (unsigned long)blocks[0]; }
			constexpr explicit operator unsigned long long() const noexcept { return (unsigned long long)blocks[0]; }

			constexpr explicit operator signed short() const noexcept { return (signed short)blocks[0]; }
			constexpr explicit operator signed int() const noexcept { return (signed int)blocks[0]; }
			constexpr explicit operator signed long() const noexcept { return (signed long)blocks[0]; }
			constexpr explicit operator signed long long() const noexcept { return (signed long long)blocks[0]; }

		public: // -- bool conversion -- //

			constexpr explicit operator bool() const noexcept
			{
				for (std::size_t i = 0; i < bits / 64; ++i) if (blocks[i]) return true;
				return false;
			}
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
				return ss.eof(); // we need to have parsed the entire string
			}
			// as try_parse() but throws std::invalid_argument on failure
			static void parse(double_int &res, std::string_view str, int base = 10)
			{
				if (!try_parse(res, str, base)) throw std::invalid_argument("failed to parse string");
			}
			// as parse() but allocates a new object to hold the result (less efficient if you're just going to assign it to a variable)
			static double_int parse(std::string_view str, int base = 10)
			{
				double_int res;
				parse(res, str, base);
				return res;
			}
		};

		// represents an arbitrarily large (signed) integer type in 2's complement
		class bigint
		{
		public: // -- data -- //

			// the dynamic array of blocks (little-endian). empty array is treated as zero.
			// at all times must be as small as possible while still preserving the correct signage of the 2's complement value.
			// users should probably not touch this directly, but is exposed because it's useful for writing extensions.
			std::vector<std::uint64_t> blocks;

			friend void _collapse(bigint&);

		public: // -- core -- //

			// constructs a new bigint with a value of zero
			bigint() noexcept = default;

			// copies the value of other
			bigint(const bigint&) = default;
			bigint &operator=(const bigint&) = default;

			// moves the value of other to this object. the moved-from object is zero after this operation.
			bigint(bigint&&) noexcept = default;
			bigint &operator=(bigint&&) noexcept = default;

		public: // -- promotion constructors -- //

			bigint(unsigned short val) : bigint((unsigned long long)val) {}
			bigint(unsigned int val) : bigint((unsigned long long)val) {}
			bigint(unsigned long val) : bigint((unsigned long long)val) {}
			bigint(unsigned long long val)
			{
				if (val)
				{
					blocks.push_back(val);
					if (val & 0x8000000000000000ull) blocks.push_back(0ull);
				}
			}

			bigint(short val) : bigint((long long)val) {}
			bigint(int val) : bigint((long long)val) {}
			bigint(long val) : bigint((long long)val) {}
			bigint(long long val)
			{
				if (val) blocks.push_back((unsigned long long)val);
			}

			template<std::uint64_t bits, bool sign>
			bigint(const double_int<bits, sign> &other) { operator=(other); }

		public: // -- promotion assignment -- //

			bigint &operator=(unsigned short val) { return operator=((unsigned long long)val); }
			bigint &operator=(unsigned int val) { return operator=((unsigned long long)val); }
			bigint &operator=(unsigned long val) { return operator=((unsigned long long)val); }
			bigint &operator=(unsigned long long val)
			{
				blocks.clear();
				if (val)
				{
					blocks.push_back(val);
					if (val & 0x8000000000000000ull) blocks.push_back(0ull);
				}
				return *this;
			}
			
			bigint &operator=(short val) { return operator=((long long)val); }
			bigint &operator=(int val) { return operator=((long long)val); }
			bigint &operator=(long val) { return operator=((long long)val); }
			bigint &operator=(long long val)
			{
				blocks.clear();
				if (val) blocks.push_back((unsigned long long)val);
				return *this;
			}

			template<std::uint64_t bits, bool sign>
			bigint &operator=(const double_int<bits, sign> &other)
			{
				blocks.assign(std::begin(other.blocks), std::end(other.blocks));
				if constexpr (!sign)
				{
					while (!blocks.empty() && blocks.back() == 0) blocks.pop_back();
					if (!blocks.empty() && (blocks.back() & 0x8000000000000000ull)) blocks.push_back(0ull);
				}
				else _collapse(*this);
				return *this;
			}

		public: // -- demotion conversion -- //

			explicit operator unsigned short() const noexcept { return (unsigned short)(blocks.empty() ? 0ull : blocks[0]); }
			explicit operator unsigned int() const noexcept { return (unsigned int)(blocks.empty() ? 0ull : blocks[0]); }
			explicit operator unsigned long() const noexcept { return (unsigned long)(blocks.empty() ? 0ull : blocks[0]); }
			explicit operator unsigned long long() const noexcept { return (unsigned long long)(blocks.empty() ? 0ull : blocks[0]); }

			explicit operator short() const noexcept { return (short)(blocks.empty() ? 0ull : blocks[0]); }
			explicit operator int() const noexcept { return (int)(blocks.empty() ? 0ull : blocks[0]); }
			explicit operator long() const noexcept { return (long)(blocks.empty() ? 0ull : blocks[0]); }
			explicit operator long long() const noexcept { return (long long)(blocks.empty() ? 0ull : blocks[0]); }

		public: // -- bool conversion -- //

			explicit operator bool() const noexcept { return !blocks.empty(); }
			friend bool operator!(const bigint &a) { return !(bool)a; }

		public: // -- static utilities -- //

			friend std::istream &operator>>(std::istream&, bigint&);
			friend bigint _pow(bigint, const bigint&);

			// attempts to parse the string into an integer value - returns true on successful parse.
			// base must be 10 (dec), 16 (hex), 8 (oct), or 0 to automatically determine base from C-style prefix in str.
			static bool try_parse(bigint &res, std::string_view str, int base = 10)
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
				return ss.eof(); // we need to have parsed the entire string
			}
			// as try_parse() but throws std::invalid_argument on failure
			static void parse(bigint &res, std::string_view str, int base = 10)
			{
				if (!try_parse(res, str, base)) throw std::invalid_argument("failed to parse string");
			}
			// as parse() but allocates a new object to hold the result (less efficient if you're just going to assign it to a variable)
			static bigint parse(std::string_view str, int base = 10)
			{
				bigint res;
				parse(res, str, base);
				return res;
			}

			// computes the result of raising a to the power of b
			static bigint pow(const bigint &a, const bigint &b) { return _pow(a, b); }
			static bigint pow(bigint &&a, const bigint &b) { return _pow(std::move(a), b); }
		};

		// -- bigint utility definitions -- //

		inline bool is_neg(const bigint &val) noexcept { return !val.blocks.empty() && (val.blocks.back() & 0x8000000000000000); }

		// collapses the value in the given bigint to make it satisfy the requirements of the blocks array
		inline void _collapse(bigint &a)
		{
			if (detail::is_neg(a))
			{
				for (std::size_t i = a.blocks.size(); i-- > 0 && a.blocks[i] == 0xffffffffffffffffull; ) a.blocks.pop_back();
				if (a.blocks.empty() || !(a.blocks.back() & 0x8000000000000000ull)) a.blocks.push_back(0xffffffffffffffffull);
			}
			else
			{
				for (std::size_t i = a.blocks.size(); i-- > 0 && a.blocks[i] == 0; ) a.blocks.pop_back();
				if (!a.blocks.empty() && (a.blocks.back() & 0x8000000000000000ull)) a.blocks.push_back(0);
			}
		}

		inline std::size_t highest_set_bit(const bigint &val) noexcept
		{
			if (val.blocks.empty()) return 0;
			if (val.blocks.back()) return (val.blocks.size() - 1) * 64 + highest_set_bit(val.blocks.back());
			else return (val.blocks.size() - 2) * 64 + highest_set_bit(val.blocks[val.blocks.size() - 2]);
		}
		inline void make_not(bigint &a)
		{
			if (a.blocks.empty()) a.blocks.push_back(0xffffffffffffffffull);
			else
			{
				for (std::uint64_t &v : a.blocks) v = ~v;
				_collapse(a);
			}
		}

		inline bool bit_test(const bigint &val, std::uint64_t bit) noexcept
		{
			std::uint64_t block = bit / 64;
			if (block >= val.blocks.size()) return detail::is_neg(val); // lookup in ghost blocks
			return (val.blocks[block] >> (bit % 64)) & 1; // lookup in declared range
		}
		inline bool bit_test_in_bounds_nonzero(const bigint &val, std::uint64_t bit) noexcept
		{
			return (val.blocks[bit / 64] >> (bit % 64)) & 1;
		}
		inline void bit_set(bigint &val, std::uint64_t bit) noexcept
		{
			std::uint64_t block = bit / 64;
			const bool neg = detail::is_neg(val);

			if (block >= val.blocks.size())
			{
				if (neg) return; // if negative in ghost block region, already set
				val.blocks.resize(block + 1, 0ul); // otherwise resize so that block is a valid index
			}
			bit %= 64;
			val.blocks[block] |= 1ull << bit;
			if (!neg && detail::is_neg(val)) val.blocks.push_back(0ull); // make sure positive things remain positive

			_collapse(val);
		}

		// -- operators -- //

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
			for (std::size_t i = 0; i < bits / 64; ++i) if (++a.blocks[i]) break;
			return a;
		}
		template<std::uint64_t bits, bool sign>
		[[nodiscard]] constexpr double_int<bits, sign> operator++(double_int<bits, sign> &a, int) noexcept { auto cpy = a; ++a; return cpy; }

		inline bigint &operator++(bigint &a)
		{
			for (std::size_t i = 1; i < a.blocks.size(); ++i) if (++a.blocks[i - 1]) { _collapse(a); return a; }

			if (a.blocks.empty()) a.blocks.push_back(1ull);
			else
			{
				std::uint64_t high = ++a.blocks.back();
				if (high == 0) a.blocks.clear();
				else if (high == 0x8000000000000000ull) a.blocks.push_back(0ull);
				else _collapse(a);
			}

			return a;
		}
		[[nodiscard]] inline bigint operator++(bigint &a, int) { bigint cpy = a; ++a; return cpy; }

		template<std::uint64_t bits, bool sign>
		constexpr double_int<bits, sign> &operator--(double_int<bits, sign> &a) noexcept
		{
			for (std::size_t i = 0; i < bits / 64; ++i)	if (a.blocks[i]--) break;
			return a;
		}
		template<std::uint64_t bits, bool sign>
		[[nodiscard]] constexpr double_int<bits, sign> operator--(double_int<bits, sign> &a, int) noexcept { auto cpy = a; --a; return cpy; }

		inline bigint &operator--(bigint &a)
		{
			for (std::size_t i = 1; i < a.blocks.size(); ++i) if (a.blocks[i - 1]--) { _collapse(a); return a; }

			if (a.blocks.empty()) a.blocks.push_back(0xffffffffffffffffull);
			else
			{
				std::uint64_t high = a.blocks.back()--;
				/* high == 0 case cannot happen because 0 is represented by empty array and is therefore handled above */
				if (high == 0x8000000000000000ull) a.blocks.push_back(0xffffffffffffffffull);
				else _collapse(a);
			}

			return a;
		}
		[[nodiscard]] inline bigint operator--(bigint &a, int) { bigint cpy = a; --a; return cpy; }

		inline void make_neg(bigint &a) { make_not(a); ++a; }

		// -- add -- //

		template<std::uint64_t bits, bool sign1, bool sign2> constexpr double_int<bits, sign1> &operator+=(double_int<bits, sign1> &a, const double_int<bits, sign2> &b) noexcept
		{
			std::uint64_t carry = 0;
			for (std::size_t i = 0; i < bits / 64; ++i)
			{
				std::uint64_t v = b.blocks[i] + carry;
				carry = (a.blocks[i] += v) < v || v < carry ? 1 : 0;
			}
			return a;
		}

		template<bool subtract>
		bigint &_add(bigint &a, const bigint &b)
		{
			std::size_t min = std::min(a.blocks.size(), b.blocks.size());
			const bool a_neg = detail::is_neg(a);
			const bool b_neg = subtract ^ detail::is_neg(b);

			// compute addition on the mutually covered range
			std::uint64_t carry = subtract ? 1 : 0;
			for (std::size_t i = 0; i < min; ++i)
			{
				std::uint64_t v = (subtract ? ~b.blocks[i] : b.blocks[i]) + carry;
				carry = (a.blocks[i] += v) < v || v < carry ? 1 : 0;
			}

			// perform addition on the extended range (ghost blocks)
			if (a.blocks.size() < b.blocks.size())
			{
				a.blocks.resize(b.blocks.size(), a_neg ? -1 : 0);
				for (std::size_t i = min; i < a.blocks.size(); ++i)
				{
					std::uint64_t v = (subtract ? ~b.blocks[i] : b.blocks[i]) + carry;
					carry = (a.blocks[i] += v) < v || v < carry ? 1 : 0;
				}
			}
			else if (a.blocks.size() > b.blocks.size())
			{
				std::size_t ghost = b_neg ? -1 : 0;
				for (std::size_t i = min; i < a.blocks.size(); ++i)
				{
					std::uint64_t v = ghost + carry;
					carry = (a.blocks[i] += v) < v || v < carry ? 1 : 0;
				}
			}

			// perform addition on the infinite prefix blocks
			if (!a_neg && !b_neg)
			{
				if (carry) a.blocks.push_back(1ull);
				else if (detail::is_neg(a)) a.blocks.push_back(0ull); // value of a has changed and might be negative now
			}
			else if (a_neg && b_neg)
			{
				if (!carry) a.blocks.push_back(0xfffffffffffffffeull);
				else if (!detail::is_neg(a)) a.blocks.push_back(0xffffffffffffffffull); // value of a has changed and might be negative now
			}
			else
			{
				if (carry)
				{
					if (detail::is_neg(a)) a.blocks.push_back(0ull);
				}
				else
				{
					if (!detail::is_neg(a)) a.blocks.push_back(0xffffffffffffffffull);
				}
			}

			_collapse(a);
			return a;
		}

		inline bigint &operator+=(bigint &a, const bigint &b) { return _add<false>(a, b); }

		inline bigint operator+(const bigint &a, const bigint &b) { bigint cpy = a; cpy += b; return cpy; }
		inline bigint operator+(bigint &&a, const bigint &b) { bigint cpy = std::move(a); cpy += b; return cpy; }
		inline bigint operator+(const bigint &a, bigint &&b) { bigint cpy = std::move(b); cpy += a; return cpy; }
		inline bigint operator+(bigint &&a, bigint &&b)
		{
			if (a.blocks.size() >= b.blocks.size()) { bigint cpy = std::move(a); cpy += b; return cpy; }
			else { bigint cpy = std::move(b); cpy += a; return cpy; }
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
			std::uint64_t carry = 1;
			for (std::size_t i = 0; i < bits / 64; ++i)
			{
				std::uint64_t v = ~b.blocks[i] + carry;
				carry = (a.blocks[i] += v) < v || v < carry ? 1 : 0;
			}
			return a;
		}

		inline bigint &operator-=(bigint &a, const bigint &b) { return _add<true>(a, b); }

		inline bigint operator-(const bigint &a, const bigint &b) { bigint cpy = a; cpy -= b; return cpy; }
		inline bigint operator-(bigint &&a, const bigint &b) { bigint cpy = std::move(a); cpy -= b; return cpy; }
		inline bigint operator-(const bigint &a, bigint &&b) { bigint cpy = std::move(b); cpy -= a; detail::make_neg(cpy); return cpy; } // not sure if this is faster than just making a dynamically-allocated copy
		inline bigint operator-(bigint &&a, bigint &&b) { bigint cpy = std::move(a); cpy -= b; return cpy; }

		template<std::uint64_t bits, bool sign, typename T>
		constexpr double_int<bits, sign> &operator-=(double_int<bits, sign> &a, const T &b) noexcept { a -= (double_int<bits, sign>)b; return a; }

		template<std::uint64_t bits, bool sign>
		constexpr double_int<bits, sign> operator-(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { double_int<bits, sign> res = a; res -= b; return res; }

		SHORTERHAND_BINARY_FORMATTER(-)

		// -- and -- //

		template<std::uint64_t bits, bool sign1, bool sign2>
		constexpr double_int<bits, sign1> &operator&=(double_int<bits, sign1> &a, const double_int<bits, sign2> &b) noexcept
		{
			for (std::size_t i = 0; i < bits / 64; ++i) a.blocks[i] &= b.blocks[i];
			return a;
		}

		inline bigint &operator&=(bigint &a, const bigint &b)
		{
			const bool a_neg = detail::is_neg(a);
			const bool b_neg = detail::is_neg(b);

			const std::size_t min = std::min(a.blocks.size(), b.blocks.size());
			
			for (std::size_t i = 0; i < min; ++i) a.blocks[i] &= b.blocks[i]; // perform on range declared by both

			if (a.blocks.size() < b.blocks.size())
			{
				           // if a was positive then rest is all 0's and thus result is already computed
				if (a_neg) // otherwise a was negative, so rest is all 1's - just copy remaining blocks from b
				{
					a.blocks.resize(b.blocks.size());
					for (std::size_t i = min; i < b.blocks.size(); ++i) a.blocks[i] = b.blocks[i];
				}
			}
			else if (b.blocks.size() < a.blocks.size())
			{
							// if b was negative then rest is all 1's - we keep everything in a that was already there
				if (!b_neg) // otherwise b was positive - rest is all 0's and thus result is already computed - just truncate down to computed region and make sure result is positive
				{
					a.blocks.resize(min);
					if (detail::is_neg(a)) a.blocks.push_back(0ull);
				}
			}

			_collapse(a);
			return a;
		}

		inline bigint operator&(const bigint &a, const bigint &b) { bigint cpy = a; cpy &= b; return cpy; }
		inline bigint operator&(bigint &&a, const bigint &b) { bigint cpy = std::move(a); cpy &= b; return cpy; }
		inline bigint operator&(const bigint &a, bigint &&b) { bigint cpy = std::move(b); cpy &= a; return cpy; }
		inline bigint operator&(bigint &&a, bigint &&b) { bigint cpy = std::move(a); cpy &= b; return cpy; }

		template<std::uint64_t bits, bool sign, typename T>
		constexpr double_int<bits, sign> &operator&=(double_int<bits, sign> &a, const T &b) noexcept { a &= (double_int<bits, sign>)b; return a; }

		template<std::uint64_t bits, bool sign>
		constexpr double_int<bits, sign> operator&(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { double_int<bits, sign> res = a; res &= b; return res; }

		SHORTERHAND_BINARY_FORMATTER(&)

		// -- or -- //

		template<std::uint64_t bits, bool sign1, bool sign2>
		constexpr double_int<bits, sign1> &operator|=(double_int<bits, sign1> &a, const double_int<bits, sign2> &b) noexcept
		{
			for (std::size_t i = 0; i < bits / 64; ++i) a.blocks[i] |= b.blocks[i];
			return a;
		}

		inline bigint &operator|=(bigint &a, const bigint &b)
		{
			const bool a_neg = detail::is_neg(a);
			const bool b_neg = detail::is_neg(b);

			const std::size_t min = std::min(a.blocks.size(), b.blocks.size());

			for (std::size_t i = 0; i < min; ++i) a.blocks[i] |= b.blocks[i]; // perform on range declared by both

			if (a.blocks.size() < b.blocks.size())
			{
							// if a was negative then rest is all 1's, and thus result is already computed
				if (!a_neg) // otherwise a positive, so rest is all 0's - just copy remaining blocks from b
				{
					a.blocks.resize(b.blocks.size());
					for (std::size_t i = min; i < b.blocks.size(); ++i) a.blocks[i] = b.blocks[i];
				}
			}
			else if (b.blocks.size() < a.blocks.size())
			{
						   // if b was positive then rest is all 0's - we keep everything in a that was already there
				if (b_neg) // otherwise b was negative, so rest is all 1's - just truncate down to computed region and make sure result is negative
				{
					a.blocks.resize(min);
					if (!detail::is_neg(a)) a.blocks.push_back(0xffffffffffffffffull);
				}
			}

			_collapse(a);
			return a;
		}

		inline bigint operator|(const bigint &a, const bigint &b) { bigint cpy = a; cpy |= b; return cpy; }
		inline bigint operator|(bigint &&a, const bigint &b) { bigint cpy = std::move(a); cpy |= b; return cpy; }
		inline bigint operator|(const bigint &a, bigint &&b) { bigint cpy = std::move(b); cpy |= a; return cpy; }
		inline bigint operator|(bigint &&a, bigint &&b) { bigint cpy = std::move(a); cpy |= b; return cpy; }

		template<std::uint64_t bits, bool sign, typename T>
		constexpr double_int<bits, sign> &operator|=(double_int<bits, sign> &a, const T &b) noexcept { a |= (double_int<bits, sign>)b; return a; }

		template<std::uint64_t bits, bool sign>
		constexpr double_int<bits, sign> operator|(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { double_int<bits, sign> res = a; res |= b; return res; }

		SHORTERHAND_BINARY_FORMATTER(| )

		// -- xor -- //

		template<std::uint64_t bits, bool sign1, bool sign2>
		constexpr double_int<bits, sign1> &operator^=(double_int<bits, sign1> &a, const double_int<bits, sign2> &b) noexcept
		{
			for (std::size_t i = 0; i < bits / 64; ++i) a.blocks[i] ^= b.blocks[i];
			return a;
		}

		inline bigint &operator^=(bigint &a, const bigint &b)
		{
			const bool a_neg = detail::is_neg(a);
			const bool b_neg = detail::is_neg(b);

			const std::size_t min = std::min(a.blocks.size(), b.blocks.size());

			for (std::size_t i = 0; i < min; ++i) a.blocks[i] ^= b.blocks[i]; // perform on range declared by both

			if (a.blocks.size() < b.blocks.size())
			{
				if (a_neg) // if a was negative then rest is all 1's - just copy the inverted bits of b
				{
					a.blocks.resize(b.blocks.size());
					for (std::size_t i = min; i < b.blocks.size(); ++i) a.blocks[i] = ~b.blocks[i];
				}
				else // otherwise a was positive, so rest is all 0's - just copy the bits of b
				{
					a.blocks.resize(b.blocks.size());
					for (std::size_t i = min; i < b.blocks.size(); ++i) a.blocks[i] = b.blocks[i];
				}
			}
			else if (b.blocks.size() < a.blocks.size())
			{
						   // if b was positive then rest is all 0's - so we already have the result
				if (b_neg) // otherwise b negative, so rest is all 1's - just invert the bits of a
				{
					for (std::size_t i = min; i < a.blocks.size(); ++i) a.blocks[i] = ~a.blocks[i];
				}
			}

			_collapse(a);
			return a;
		}

		inline bigint operator^(const bigint &a, const bigint &b) { bigint cpy = a; cpy ^= b; return cpy; }
		inline bigint operator^(bigint &&a, const bigint &b) { bigint cpy = std::move(a); cpy ^= b; return cpy; }
		inline bigint operator^(const bigint &a, bigint &&b) { bigint cpy = std::move(b); cpy ^= a; return cpy; }
		inline bigint operator^(bigint &&a, bigint &&b) { bigint cpy = std::move(a); cpy ^= b; return cpy; }

		template<std::uint64_t bits, bool sign, typename T>
		constexpr double_int<bits, sign> &operator^=(double_int<bits, sign> &a, const T &b) noexcept { a ^= (double_int<bits, sign>)b; return a; }

		template<std::uint64_t bits, bool sign>
		constexpr double_int<bits, sign> operator^(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { auto res = a; res ^= b; return res; }

		SHORTERHAND_BINARY_FORMATTER(^)

		// -- shl/sal -- //

		template<std::uint64_t bits, bool sign>
		constexpr double_int<bits, sign> &operator<<=(double_int<bits, sign> &val, std::uint64_t count) noexcept
		{
			count &= bits - 1;
			const std::size_t full = count / 64;
			if (full)
			{
				for (std::size_t i = bits / 64; i-- > full; ) val.blocks[i] = val.blocks[i - full];
				for (std::size_t i = full; i-- > 0; ) val.blocks[i] = 0;
				count &= 63;
			}
			if (count)
			{
				for (std::size_t i = bits / 64 - 1; i > full; --i) val.blocks[i] = (val.blocks[i] << count) | (val.blocks[i - 1] >> (64 - count));
				val.blocks[full] <<= count;
			}
			return val;
		}

		inline bigint &operator<<=(bigint &val, std::uint64_t count)
		{
			if (!val.blocks.empty())
			{
				const std::uint64_t fill = detail::is_neg(val) ? -1 : 0;
				const std::size_t full = count / 64;
				if (full)
				{
					val.blocks.resize(val.blocks.size() + full); // make room for the shifting blocks
					for (std::size_t i = val.blocks.size(); i-- > full; ) val.blocks[i] = val.blocks[i - full];
					for (std::size_t i = full; i-- > 0; ) val.blocks[i] = 0;
					count &= 63;
				}
				if (count)
				{
					std::uint64_t overflow_block = (fill << count) | (val.blocks.back() >> (64 - count)); // compute overflow block from high block (include ghost bits based on sign)
					for (std::size_t i = val.blocks.size() - 1; i > full; --i) val.blocks[i] = (val.blocks[i] << count) | (val.blocks[i - 1] >> (64 - count));
					val.blocks[full] <<= count;

					const std::uint64_t new_fill = detail::is_neg(val) ? -1 : 0; // get new fill pattern
					if (overflow_block != new_fill) val.blocks.push_back(overflow_block); // if the overflow block isn't generated by the new fill pattern, we need to append it
				}
			}
			return val;
		}

		inline bigint operator<<(const bigint &val, std::uint64_t count) { bigint cpy = val; cpy <<= count; return cpy; }
		inline bigint operator<<(bigint &&val, std::uint64_t count) { bigint cpy = std::move(val); cpy <<= count; return cpy; }

		template<std::uint64_t bits, bool sign, typename T>
		constexpr double_int<bits, sign> operator<<(const double_int<bits, sign> &a, const T &count) noexcept { double_int<bits, sign> res = a; res <<= (std::uint64_t)count; return res; }

		// -- shr/sar -- //

		// shr
		template<std::uint64_t bits>
		constexpr double_int<bits, false> &operator>>=(double_int<bits, false> &val, std::uint64_t count) noexcept
		{
			count &= bits - 1;
			const std::size_t full = count / 64;
			if (full)
			{
				for (std::size_t i = 0; i < bits / 64 - full; ++i) val.blocks[i] = val.blocks[i + full];
				for (std::size_t i = bits / 64 - full; i < bits / 64; ++i) val.blocks[i] = 0;
				count &= 63;
			}
			if (count)
			{
				for (std::size_t i = full; i < bits / 64 - 1; ++i) val.blocks[i] = (val.blocks[i] >> count) | (val.blocks[i + 1] << (64 - count));
				val.blocks[bits / 64 - 1] >>= count;
			}
			return val;
		}

		// sar
		template<std::uint64_t bits>
		constexpr double_int<bits, true> &operator>>=(double_int<bits, true> &val, std::uint64_t count) noexcept
		{
			count &= bits - 1;
			const std::uint64_t fill = (val.blocks[bits / 64 - 1] & 0x8000000000000000) ? -1 : 0;
			const std::size_t full = count / 64;
			if (full)
			{
				for (std::size_t i = 0; i < bits / 64 - full; ++i) val.blocks[i] = val.blocks[i + full];
				for (std::size_t i = bits / 64 - full; i < bits / 64; ++i) val.blocks[i] = fill;
				count &= 63;
			}
			if (count)
			{
				for (std::size_t i = full; i < bits / 64 - 1; ++i) val.blocks[i] = (val.blocks[i] >> count) | (val.blocks[i + 1] << (64 - count));
				val.blocks[bits / 64 - 1] = (val.blocks[bits / 64 - 1] >> count) | (fill << (64 - count));
			}
			return val;
		}

		inline bigint &operator>>=(bigint &val, std::uint64_t count)
		{
			if (!val.blocks.empty())
			{
				const std::uint64_t fill = detail::is_neg(val) ? -1 : 0;
				const std::size_t full = count / 64;
				if (full)
				{
					for (std::size_t i = 0; i < val.blocks.size() - full; ++i) val.blocks[i] = val.blocks[i + full];
					for (std::size_t i = val.blocks.size() - full; i < val.blocks.size(); ++i) val.blocks[i] = fill;
					count &= 63;
				}
				if (count)
				{
					for (std::size_t i = full; i < val.blocks.size() - 1; ++i) val.blocks[i] = (val.blocks[i] >> count) | (val.blocks[i + 1] << (64 - count));
					val.blocks[val.blocks.size() - 1] = (val.blocks[val.blocks.size() - 1] >> count) | (fill << (64 - count));
				}
				_collapse(val);
			}
			return val;
		}

		inline bigint operator>>(const bigint &val, std::uint64_t count) { bigint cpy = val; cpy >>= count; return cpy; }
		inline bigint operator>>(bigint &&val, std::uint64_t count) { bigint cpy = std::move(val); cpy >>= count; return cpy; }

		template<std::uint64_t bits, bool sign, typename T>
		constexpr double_int<bits, sign> operator>>(const double_int<bits, sign> &a, const T &count) noexcept { double_int<bits, sign> res = a; res >>= (std::uint64_t)count; return res; }

		// -- unary -- //

		template<std::uint64_t bits, bool sign>
		constexpr double_int<bits, sign> operator+(const double_int<bits, sign> &a) noexcept { return a; }

		inline bigint operator+(const bigint &a) { return a; }
		inline bigint operator+(bigint &&a) { return std::move(a); }

		template<std::uint64_t bits>
		constexpr double_int<bits, true> operator-(const double_int<bits, true> &a) noexcept { double_int<bits, true> res = a; detail::make_neg(res); return res; }

		inline bigint operator-(const bigint &a) { bigint res = a; detail::make_neg(res); return res; }
		inline bigint operator-(bigint &&a) { bigint res = std::move(a); detail::make_neg(res); return res; }

		template<std::uint64_t bits, bool sign>
		constexpr double_int<bits, sign> operator~(const double_int<bits, sign> &a) noexcept { double_int<bits, sign> res = a; detail::make_not(res); return res; }

		inline bigint operator~(const bigint &a) { bigint res = a; detail::make_not(res); return res; }
		inline bigint operator~(bigint &&a) { bigint res = std::move(a); detail::make_not(res); return res; }

		// -- mul -- //

		template<std::uint64_t bits, bool sign>
		constexpr double_int<bits, sign> _multiply(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept
		{
			double_int<bits, sign> res = 0ull; // the final result object - zero initialized for the add logic later

			for (std::size_t i = 0; i < bits / 64; ++i) // loop through each block in a
			{
				// here we construct temp_1, which is the value of a->blocks[i] * b, scaled by i words (think long multiplication from grade school)
				for (std::size_t j = i; j < bits / 64 - 1; ++j)
				{
					auto p = detail::_mul_u64_fast(a.blocks[i], b.blocks[j - i]); // compute the product with this block from b
					if ((res.blocks[j] += p.first) < p.first)
					{
						for (std::size_t k = j + 1; k < bits / 64 && !++res.blocks[k]; ++k);
					}
					if ((res.blocks[j + 1] += p.second) < p.second)
					{
						for (std::size_t k = j + 2; k < bits / 64 && !++res.blocks[k]; ++k);
					}
				}
				res.blocks[bits / 64 - 1] += a.blocks[i] * b.blocks[bits / 64 - 1 - i]; // last iteration of the loop is special since it overflows the (finite) array
			}
			return res;
		}
		inline bigint _multiply_positive(const bigint &a, const bigint &b)
		{
			const std::size_t a_size = a.blocks.size();
			const std::size_t b_size = b.blocks.size();

			bigint res; // the final result object - bastardized for efficiency right now
			res.blocks.resize(a_size + b_size, 0ull); // allocate all the space we'll need, zero initialized for the add logic later

			for (std::size_t i = 0; i < a_size; ++i) // loop through each block in a
			{
				// here we construct temp_1, which is the value of a->blocks[i] * b, scaled by i words (think long multiplication from grade school)
				for (std::size_t j = 0; j < b_size; ++j)
				{
					auto p = detail::_mul_u64_fast(a.blocks[i], b.blocks[j]); // compute the product with this block from b
					if ((res.blocks[j + i] += p.first) < p.first) // add the low half to the corresponding block in res - overflow propagates a carry up to higher blocks
					{
						for (std::size_t k = j + i + 1; !++res.blocks[k]; ++k);
					}
					if ((res.blocks[j + i + 1] += p.second) < p.second) // add the high half to the next block in res - overflow propagates a carry up to higher blocks
					{
						for (std::size_t k = j + i + 2; !++res.blocks[k]; ++k);
					}
				}
			}

			detail::_collapse(res); // we need to perform one collapse operation on the finished result to put it into a valid state
			return res;
		}
		template<typename U, typename V, std::enable_if_t<std::is_same_v<std::decay_t<U>, bigint> && std::is_same_v<std::decay_t<V>, bigint>, int> = 0>
		inline bigint _multiply_unknown(U &&a, V &&b)
		{
			const bool a_neg = detail::is_neg(a);
			const bool b_neg = detail::is_neg(b);

			if (a_neg && b_neg)
			{
				bigint a_cpy = std::forward<U>(a);
				bigint b_cpy = std::forward<V>(b);
				detail::make_neg(a_cpy);
				detail::make_neg(b_cpy);

				return detail::_multiply_positive(a_cpy, b_cpy);
			}
			else if (a_neg)
			{
				bigint a_cpy = std::forward<U>(a);
				detail::make_neg(a_cpy);

				bigint res = detail::_multiply_positive(a_cpy, b);
				detail::make_neg(res);
				return res;
			}
			else if (b_neg)
			{
				bigint b_cpy = std::forward<V>(b);
				detail::make_neg(b_cpy);

				bigint res = detail::_multiply_positive(a, b_cpy);
				detail::make_neg(res);
				return res;
			}
			else return detail::_multiply_positive(a, b);
		}

		inline bigint operator*(const bigint &a, const bigint &b) { return detail::_multiply_unknown(a, b); }
		inline bigint operator*(bigint &&a, const bigint &b) { return detail::_multiply_unknown(std::move(a), b); }
		inline bigint operator*(const bigint &a, bigint &&b) { return detail::_multiply_unknown(a, std::move(b)); }
		inline bigint operator*(bigint &&a, bigint &&b) { return detail::_multiply_unknown(std::move(a), std::move(b)); }

		template<std::uint64_t bits, bool sign>
		constexpr double_int<bits, sign> operator*(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { return detail::_multiply(a, b); }

		SHORTERHAND_BINARY_FORMATTER(*)

		inline bigint &operator*=(bigint &a, const bigint &b)
		{
			if (&a != &b) a = std::move(a) * b;
			else a = a * b;
			return a;
		}
		inline bigint &operator*=(bigint &a, bigint &&b) { a = std::move(a) * std::move(b); return a; }

		template<std::uint64_t bits, bool sign>
		constexpr double_int<bits, sign> &operator*=(double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { a = a * b; return a; }

		template<std::uint64_t bits, bool sign, typename T>
		constexpr double_int<bits, sign> &operator*=(double_int<bits, sign> &a, const T &b) noexcept { a *= (double_int<bits, sign>)b; return a; }

		inline int cmp(const bigint &a, const bigint &b) noexcept
		{
			const bool a_neg = detail::is_neg(a);
			const bool b_neg = detail::is_neg(b);
			if (a_neg ^ b_neg) return a_neg ? -1 : 1; // if they have different signs we know right away

			if (a.blocks.size() > b.blocks.size()) return a_neg ? -1 : 1; // if a has higher magnitude than b we can tell from the sign
			if (a.blocks.size() < b.blocks.size()) return a_neg ? 1 : -1; // similarly

			for (std::size_t i = a.blocks.size(); i-- > 0; ) if (a.blocks[i] != b.blocks[i]) return a.blocks[i] < b.blocks[i] ? -1 : 1; // otherwise base decision on first different block

			return 0; // otherwise they're equal
		}

		inline bool operator<(const bigint &a, const bigint &b) noexcept { return cmp(a, b) < 0; }
		inline bool operator<=(const bigint &a, const bigint &b) noexcept { return cmp(a, b) <= 0; }
		inline bool operator>(const bigint &a, const bigint &b) noexcept { return cmp(a, b) > 0; }
		inline bool operator>=(const bigint &a, const bigint &b) noexcept { return cmp(a, b) >= 0; }

		inline std::pair<bigint, bigint> divmod_unchecked_positive(const bigint &num, const bigint &den)
		{
			std::pair<bigint, bigint> res; // default constructed to (0, 0)

			std::uint64_t bit = highest_set_bit(num);

			const std::size_t den_highest_bit = highest_set_bit(den);
			std::size_t res_second_highest_bit = 0;
			std::size_t shift_count = 0;

			while (true)
			{
				++shift_count;
				if (bit_test(num, bit)) // no unchecked equivalent because arbitrary precision
				{
					res.second <<= shift_count;
					res_second_highest_bit += shift_count;
					shift_count = 0;
					++res.second;
				}
				else if (res_second_highest_bit + shift_count >= den_highest_bit)
				{
					res.second <<= shift_count;
					res_second_highest_bit += shift_count;
					shift_count = 0;
				}

				if (res.second >= den)
				{
					res.second -= den;
					res_second_highest_bit = highest_set_bit(res.second);
					bit_set(res.first, bit); // no unchecked equivalent because arbitrary precision
				}

				if (bit-- == 0) break;
			}

			res.second <<= shift_count;

			return res;
		}
		template<typename U, typename V, std::enable_if_t<std::is_same_v<std::decay_t<U>, bigint> && std::is_same_v<std::decay_t<V>, bigint>, int> = 0>
		inline std::pair<bigint, bigint> divmod_unchecked(U &&a, V &&b)
		{
			const bool a_neg = detail::is_neg(a);
			const bool b_neg = detail::is_neg(b);

			if (a_neg && b_neg)
			{
				bigint a_cpy = std::forward<U>(a);
				bigint b_cpy = std::forward<V>(b);
				detail::make_neg(a_cpy);
				detail::make_neg(b_cpy);

				return detail::divmod_unchecked_positive(a_cpy, b_cpy);
			}
			else if (a_neg)
			{
				bigint a_cpy = std::forward<U>(a);
				detail::make_neg(a_cpy);

				auto res = detail::divmod_unchecked_positive(a_cpy, b);
				detail::make_neg(res.first);
				detail::make_neg(res.second);
				return res;
			}
			else if (b_neg)
			{
				bigint b_cpy = std::forward<V>(b);
				detail::make_neg(b_cpy);

				auto res = detail::divmod_unchecked_positive(a, b_cpy);
				detail::make_neg(res.first);
				detail::make_neg(res.second);
				return res;
			}
			else return detail::divmod_unchecked_positive(a, b);
		}
		template<typename U, typename V, std::enable_if_t<std::is_same_v<std::decay_t<U>, bigint> && std::is_same_v<std::decay_t<V>, bigint>, int> = 0>
		inline std::pair<bigint, bigint> divmod(U &&a, V &&b)
		{
			if (!b) throw std::domain_error("divide by zero");
			return divmod_unchecked(std::forward<U>(a), std::forward<V>(b));
		}

		inline bigint operator/(const bigint &num, const bigint &den) { return detail::divmod(num, den).first; }
		inline bigint operator/(bigint &&num, const bigint &den) { return detail::divmod(std::move(num), den).first; }
		inline bigint operator/(const bigint &num, bigint &&den) { return detail::divmod(num, std::move(den)).first; }
		inline bigint operator/(bigint &&num, bigint &&den) { return detail::divmod(std::move(num), std::move(den)).first; }

		inline bigint operator%(const bigint &num, const bigint &den) { return detail::divmod(num, den).second; }
		inline bigint operator%(bigint &&num, const bigint &den) { return detail::divmod(std::move(num), den).second; }
		inline bigint operator%(const bigint &num, bigint &&den) { return detail::divmod(num, std::move(den)).second; }
		inline bigint operator%(bigint &&num, bigint &&den) { return detail::divmod(std::move(num), std::move(den)).second; }

		inline bigint &operator/=(bigint &num, const bigint &den)
		{
			if (&num != &den) num = std::move(num) / den;
			else num = num / den;
			return num;
		}
		inline bigint &operator/=(bigint &num, bigint &&den) { num = std::move(num) / std::move(den); return num; }

		inline bigint &operator%=(bigint &num, const bigint &den)
		{
			if (&num != &den) num = std::move(num) % den;
			else num = num % den;
			return num;
		}
		inline bigint &operator%=(bigint &num, bigint &&den) { num = std::move(num) % std::move(den); return num; }

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
			if (detail::is_neg(unum)) { detail::make_neg(unum); neg = !neg; }
			if (detail::is_neg(uden)) { detail::make_neg(uden); neg = !neg; }

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

		SHORTERHAND_BINARY_FORMATTER(/)
		SHORTERHAND_BINARY_FORMATTER(%)

		// -- cmp -- //

		// cmp() will be equivalent to <=> but works for any version of C++

		template<std::uint64_t bits_1, std::uint64_t bits_2, bool sign>
		constexpr int cmp(const double_int<bits_1, sign> &a, const double_int<bits_2, sign> &b) noexcept
		{
			if constexpr (bits_2 > bits_1) return -cmp(b, a); // wolog let a be at least as large as b (in physical size)
			else
			{
				if constexpr (sign)
				{
					const bool a_neg = detail::is_neg(a);
					const bool b_neg = detail::is_neg(b);
					if (a_neg ^ b_neg) return a_neg ? -1 : 1; // if they have different signs we know right away
				}

				if constexpr (bits_1 > bits_2)
				{
					std::uint64_t fill = 0;
					if constexpr (sign) fill = detail::is_neg(b) ? 0xffffffffffffffffull : 0ull;
					for (std::size_t i = bits_1 / 64; i-- > bits_2 / 64; ) if (a.blocks[i] != fill) // if a has higher magnitude than b we can tell (immediately for unsigned, or from sign bit for signed)
					{
						if constexpr (!sign) return 1;
						else return detail::is_neg(a) ? -1 : 1;
					}
				}

				for (std::size_t i = bits_2 / 64; i-- > 0; ) if (a.blocks[i] != b.blocks[i]) return a.blocks[i] < b.blocks[i] ? -1 : 1; // otherwise base becision on first different block

				return 0; // otherwise they're equal
			}
		}
		template<std::uint64_t bits, bool sign>
		constexpr int cmp(const double_int<bits, sign> &a, std::conditional_t<sign, long long, unsigned long long> val) noexcept
		{
			std::uint64_t fill = 0;
			if constexpr (sign) fill = val < 0 ? 0xffffffffffffffffull : 0ull;
			for (std::size_t i = bits / 64 - 1; i > 0; --i) if (a.blocks[i] != fill)
			{
				if constexpr (!sign) return 1;
				else return detail::is_neg(a) ? -1 : 1;
			}
			if (a.blocks[0] != (std::uint64_t)val)
			{
				if constexpr (!sign) return a.blocks[0] < (std::uint64_t)val ? -1 : 1;
				else return a.blocks[0] < (std::uint64_t)val ? 1 : -1;
			}
			return 0;
		}
		template<std::uint64_t bits, bool sign>
		constexpr int cmp(const double_int<bits, sign> &a, std::conditional_t<sign, long, unsigned long> val) noexcept { return cmp(a, (std::conditional_t<sign, long long, unsigned long long>)val); }
		template<std::uint64_t bits, bool sign>
		constexpr int cmp(const double_int<bits, sign> &a, std::conditional_t<sign, int, unsigned int> val) noexcept { return cmp(a, (std::conditional_t<sign, long long, unsigned long long>)val); }
		template<std::uint64_t bits, bool sign>
		constexpr int cmp(const double_int<bits, sign> &a, std::conditional_t<sign, short, unsigned short> val) noexcept { return cmp(a, (std::conditional_t<sign, long long, unsigned long long>)val); }

		template<std::uint64_t bits, bool sign>
		constexpr int cmp(std::conditional_t<sign, long long, unsigned long long> val, const double_int<bits, sign> &a) noexcept { return -cmp(a, (std::conditional_t<sign, long long, unsigned long long>)val); }
		template<std::uint64_t bits, bool sign>
		constexpr int cmp(std::conditional_t<sign, long, unsigned long> val, const double_int<bits, sign> &a) noexcept { return -cmp(a, (std::conditional_t<sign, long long, unsigned long long>)val); }
		template<std::uint64_t bits, bool sign>
		constexpr int cmp(std::conditional_t<sign, int, unsigned int> val, const double_int<bits, sign> &a) noexcept { return -cmp(a, (std::conditional_t<sign, long long, unsigned long long>)val); }
		template<std::uint64_t bits, bool sign>
		constexpr int cmp(std::conditional_t<sign, short, unsigned short> val, const double_int<bits, sign> &a) noexcept { return -cmp(a, (std::conditional_t<sign, long long, unsigned long long>)val); }

		template<std::uint64_t bits, bool sign> constexpr int cmp(const double_int<bits, sign> &a, std::conditional_t<!sign, long long, unsigned long long> val) noexcept = delete;
		template<std::uint64_t bits, bool sign> constexpr int cmp(const double_int<bits, sign> &a, std::conditional_t<!sign, long, unsigned long> val) noexcept = delete;
		template<std::uint64_t bits, bool sign> constexpr int cmp(const double_int<bits, sign> &a, std::conditional_t<!sign, int, unsigned int> val) noexcept = delete;
		template<std::uint64_t bits, bool sign> constexpr int cmp(const double_int<bits, sign> &a, std::conditional_t<!sign, short, unsigned short> val) noexcept = delete;

		template<std::uint64_t bits, bool sign> constexpr int cmp(std::conditional_t<!sign, long long, unsigned long long> val, const double_int<bits, sign> &a) noexcept = delete;
		template<std::uint64_t bits, bool sign> constexpr int cmp(std::conditional_t<!sign, long, unsigned long> val, const double_int<bits, sign> &a) noexcept = delete;
		template<std::uint64_t bits, bool sign> constexpr int cmp(std::conditional_t<!sign, int, unsigned int> val, const double_int<bits, sign> &a) noexcept = delete;
		template<std::uint64_t bits, bool sign> constexpr int cmp(std::conditional_t<!sign, short, unsigned short> val, const double_int<bits, sign> &a) noexcept = delete;
		
		template<std::uint64_t bits1, std::uint64_t bits2, bool sign> constexpr bool operator==(const double_int<bits1, sign> &a, const double_int<bits2, sign> &b) noexcept { return cmp(a, b) == 0; }
		template<std::uint64_t bits1, std::uint64_t bits2, bool sign> constexpr bool operator!=(const double_int<bits1, sign> &a, const double_int<bits2, sign> &b) noexcept { return cmp(a, b) != 0; }
		template<std::uint64_t bits1, std::uint64_t bits2, bool sign> constexpr bool operator<(const double_int<bits1, sign> &a, const double_int<bits2, sign> &b) noexcept { return cmp(a, b) < 0; }
		template<std::uint64_t bits1, std::uint64_t bits2, bool sign> constexpr bool operator<=(const double_int<bits1, sign> &a, const double_int<bits2, sign> &b) noexcept { return cmp(a, b) <= 0; }
		template<std::uint64_t bits1, std::uint64_t bits2, bool sign> constexpr bool operator>(const double_int<bits1, sign> &a, const double_int<bits2, sign> &b) noexcept { return cmp(a, b) > 0; }
		template<std::uint64_t bits1, std::uint64_t bits2, bool sign> constexpr bool operator>=(const double_int<bits1, sign> &a, const double_int<bits2, sign> &b) noexcept { return cmp(a, b) >= 0; }

		#define DRAGAZO_BIGGERINTS_DI_BI_CMP_SHORTHAND(T) \
		template<std::uint64_t bits, bool sign> constexpr bool operator==(const double_int<bits, sign> &a, std::conditional_t<sign, T, unsigned T> b) { return cmp(a, b) == 0; } \
		template<std::uint64_t bits, bool sign> constexpr bool operator!=(const double_int<bits, sign> &a, std::conditional_t<sign, T, unsigned T> b) { return cmp(a, b) != 0; } \
		template<std::uint64_t bits, bool sign> constexpr bool operator<(const double_int<bits, sign> &a, std::conditional_t<sign, T, unsigned T> b) { return cmp(a, b) < 0; } \
		template<std::uint64_t bits, bool sign> constexpr bool operator<=(const double_int<bits, sign> &a, std::conditional_t<sign, T, unsigned T> b) { return cmp(a, b) <= 0; } \
		template<std::uint64_t bits, bool sign> constexpr bool operator>(const double_int<bits, sign> &a, std::conditional_t<sign, T, unsigned T> b) { return cmp(a, b) > 0; } \
		template<std::uint64_t bits, bool sign> constexpr bool operator>=(const double_int<bits, sign> &a, std::conditional_t<sign, T, unsigned T> b) { return cmp(a, b) >= 0; } \
		\
		template<std::uint64_t bits, bool sign> constexpr bool operator==(std::conditional_t<sign, T, unsigned T> a, const double_int<bits, sign> &b) { return cmp(a, b) == 0; } \
		template<std::uint64_t bits, bool sign> constexpr bool operator!=(std::conditional_t<sign, T, unsigned T> a, const double_int<bits, sign> &b) { return cmp(a, b) != 0; } \
		template<std::uint64_t bits, bool sign> constexpr bool operator<(std::conditional_t<sign, T, unsigned T> a, const double_int<bits, sign> &b) { return cmp(a, b) < 0; } \
		template<std::uint64_t bits, bool sign> constexpr bool operator<=(std::conditional_t<sign, T, unsigned T> a, const double_int<bits, sign> &b) { return cmp(a, b) <= 0; } \
		template<std::uint64_t bits, bool sign> constexpr bool operator>(std::conditional_t<sign, T, unsigned T> a, const double_int<bits, sign> &b) { return cmp(a, b) > 0; } \
		template<std::uint64_t bits, bool sign> constexpr bool operator>=(std::conditional_t<sign, T, unsigned T> a, const double_int<bits, sign> &b) { return cmp(a, b) >= 0; } \
		\
		template<std::uint64_t bits, bool sign> constexpr bool operator==(const double_int<bits, sign> &a, std::conditional_t<!sign, T, unsigned T> b) = delete; \
		template<std::uint64_t bits, bool sign> constexpr bool operator!=(const double_int<bits, sign> &a, std::conditional_t<!sign, T, unsigned T> b) = delete; \
		template<std::uint64_t bits, bool sign> constexpr bool operator<(const double_int<bits, sign> &a, std::conditional_t<!sign, T, unsigned T> b) = delete; \
		template<std::uint64_t bits, bool sign> constexpr bool operator<=(const double_int<bits, sign> &a, std::conditional_t<!sign, T, unsigned T> b) = delete; \
		template<std::uint64_t bits, bool sign> constexpr bool operator>(const double_int<bits, sign> &a, std::conditional_t<!sign, T, unsigned T> b) = delete; \
		template<std::uint64_t bits, bool sign> constexpr bool operator>=(const double_int<bits, sign> &a, std::conditional_t<!sign, T, unsigned T> b) = delete; \
		\
		template<std::uint64_t bits, bool sign> constexpr bool operator==(std::conditional_t<!sign, T, unsigned T> a, const double_int<bits, sign> &b) = delete; \
		template<std::uint64_t bits, bool sign> constexpr bool operator!=(std::conditional_t<!sign, T, unsigned T> a, const double_int<bits, sign> &b) = delete; \
		template<std::uint64_t bits, bool sign> constexpr bool operator<(std::conditional_t<!sign, T, unsigned T> a, const double_int<bits, sign> &b) = delete; \
		template<std::uint64_t bits, bool sign> constexpr bool operator<=(std::conditional_t<!sign, T, unsigned T> a, const double_int<bits, sign> &b) = delete; \
		template<std::uint64_t bits, bool sign> constexpr bool operator>(std::conditional_t<!sign, T, unsigned T> a, const double_int<bits, sign> &b) = delete; \
		template<std::uint64_t bits, bool sign> constexpr bool operator>=(std::conditional_t<!sign, T, unsigned T> a, const double_int<bits, sign> &b) = delete;

		DRAGAZO_BIGGERINTS_DI_BI_CMP_SHORTHAND(long long)
		DRAGAZO_BIGGERINTS_DI_BI_CMP_SHORTHAND(long)
		DRAGAZO_BIGGERINTS_DI_BI_CMP_SHORTHAND(int)
		DRAGAZO_BIGGERINTS_DI_BI_CMP_SHORTHAND(short)

		inline bool operator==(const bigint &a, const bigint &b) noexcept { return a.blocks == b.blocks; }
		inline bool operator!=(const bigint &a, const bigint &b) noexcept { return a.blocks != b.blocks; }

		inline bool operator==(const bigint &a, long long val) noexcept
		{
			return (val == 0 && a.blocks.empty()) || (a.blocks.size() == 1 && a.blocks[0] == (unsigned long long)val);
		}
		inline bool operator==(const bigint &a, long val) noexcept { return a == (long long)val; }
		inline bool operator==(const bigint &a, int val) noexcept { return a == (long long)val; }
		inline bool operator==(const bigint &a, short val) noexcept { return a == (long long)val; }

		inline bool operator==(long long val, const bigint &a) noexcept { return a == val; }
		inline bool operator==(long val, const bigint &a) noexcept { return a == val; }
		inline bool operator==(int val, const bigint &a) noexcept { return a == val; }
		inline bool operator==(short val, const bigint &a) noexcept { return a == val; }

		inline bool operator==(const bigint &a, unsigned long long val) noexcept
		{
			return (val == 0 && a.blocks.empty()) || (a.blocks.size() == 1 && a.blocks[0] == val) || (a.blocks.size() == 2 && a.blocks[0] == val && a.blocks[1] == 0);
		}
		inline bool operator==(const bigint &a, unsigned long val) noexcept { return a == (unsigned long long)val; }
		inline bool operator==(const bigint &a, unsigned int val) noexcept { return a == (unsigned long long)val; }
		inline bool operator==(const bigint &a, unsigned short val) noexcept { return a == (unsigned long long)val; }

		inline bool operator==(unsigned long long val, const bigint &a) noexcept { return a == val; }
		inline bool operator==(unsigned long val, const bigint &a) noexcept { return a == val; }
		inline bool operator==(unsigned int val, const bigint &a) noexcept { return a == val; }
		inline bool operator==(unsigned short val, const bigint &a) noexcept { return a == val; }

		inline bool operator!=(const bigint &a, long long val) noexcept { return !(a == val); }
		inline bool operator!=(const bigint &a, long val) noexcept { return !(a == val); }
		inline bool operator!=(const bigint &a, int val) noexcept { return !(a == val); }
		inline bool operator!=(const bigint &a, short val) noexcept { return !(a == val); }

		inline bool operator!=(long long val, const bigint &a) noexcept { return !(a == val); }
		inline bool operator!=(long val, const bigint &a) noexcept { return !(a == val); }
		inline bool operator!=(int val, const bigint &a) noexcept { return !(a == val); }
		inline bool operator!=(short val, const bigint &a) noexcept { return !(a == val); }

		inline bool operator!=(const bigint &a, unsigned long long val) noexcept { return !(a == val); }
		inline bool operator!=(const bigint &a, unsigned long val) noexcept { return !(a == val); }
		inline bool operator!=(const bigint &a, unsigned int val) noexcept { return !(a == val); }
		inline bool operator!=(const bigint &a, unsigned short val) noexcept { return !(a == val); }

		inline bool operator!=(unsigned long long val, const bigint &a) noexcept { return !(a == val); }
		inline bool operator!=(unsigned long val, const bigint &a) noexcept { return !(a == val); }
		inline bool operator!=(unsigned int val, const bigint &a) noexcept { return !(a == val); }
		inline bool operator!=(unsigned short val, const bigint &a) noexcept { return !(a == val); }

		// -- io -- //

		// given a hex character, converts it to an integer [0, 15] - returns true if it was a valid hex digit. ch is only meaningful on success.
		inline constexpr bool ext_hex(int &ch) noexcept
		{
			if (ch >= '0' && ch <= '9') { ch -= '0'; return true; }
			ch |= 32;
			if (ch >= 'a' && ch <= 'f') { ch = ch - 'a' + 10; return true; }
			return false;
		}

		template<typename T> static inline const T stringify_decimal_base = 10000000000000000000ull; // base to use for decimal stringification

		template<typename T, std::enable_if_t<detail::is_biggerints_type<std::decay_t<T>>::value, int> = 0>
		std::ostream &print_positive(std::ostream &ostr, T &&val)
		{
			std::string str;
			int digit, dcount;
			std::uint64_t block;

			std::ostream::sentry sentry(ostr);
			if (!sentry) return ostr;

			// build the string
			if (ostr.flags() & std::ios::oct)
			{
				auto cpy = std::forward<T>(val);

				while (true)
				{
					// get a block
					if constexpr (std::is_same_v<std::decay_t<T>, bigint>) block = cpy.blocks.empty() ? 0ull : cpy.blocks[0];
					else block = cpy.blocks[0];
					block &= 0777777777777777777777ull;
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
				auto cpy = std::forward<T>(val);
				const char hex_alpha = ostr.flags() & std::ios::uppercase ? 'A' : 'a';

				while (true)
				{
					// get a block
					if constexpr (std::is_same_v<std::decay_t<T>, bigint>) block = cpy.blocks.empty() ? 0ull : cpy.blocks[0];
					else block = cpy.blocks[0];
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
				std::pair<std::decay_t<T>, std::decay_t<T>> temp;
				temp.first = std::forward<T>(val); // temp.first takes the role of cpy in the other radicies

				while (true)
				{
					// get a block
					temp = detail::divmod_unchecked(std::move(temp.first), stringify_decimal_base<std::decay_t<T>>);
					if constexpr (std::is_same_v<std::decay_t<T>, bigint>) block = temp.second.blocks.empty() ? 0ull : temp.second.blocks[0];
					else block = temp.second.blocks[0];
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
		std::ostream &operator<<(std::ostream &ostr, const double_int<bits, false> &val) { print_positive(ostr, val); return ostr; }
		template<std::uint64_t bits>
		std::ostream &operator<<(std::ostream &ostr, const double_int<bits, true> &val)
		{
			double_int<bits, false> cpy = val; // we'll do signed io in terms of unsigned

			if (detail::is_neg(cpy)) // if it's negative make it positive and print the - sign
			{
				detail::make_neg(cpy);
				ostr.put('-');
				if (ostr.width() > 0) ostr.width(ostr.width() - 1); // dec width if was specified
			}
			else if (ostr.flags() & std::ios::showpos) // otherwise, is positive. print the + sign if showpos is set
			{
				ostr.put('+');
				if (ostr.width() > 0) ostr.width(ostr.width() - 1); // dec width if was specified
			}

			ostr << cpy; // then just refer to cpy as positive
			return ostr;
		}

		template<typename T, std::enable_if_t<std::is_same_v<std::decay_t<T>, bigint>, int> = 0>
		std::ostream &operator<<(std::ostream &ostr, T &&val)
		{
			if (detail::is_neg(val))
			{
				bigint cpy = std::forward<T>(val); // if it's negative, we need to make it positive to print it with the helper function
				detail::make_neg(cpy);
				ostr.put('-'); // put the minus sign to denote that it was negative
				if (ostr.width() > 0) ostr.width(ostr.width() - 1); // dec width if was specified
				print_positive(ostr, std::move(cpy)); // then print the positive value
			}
			else
			{
				if (ostr.flags() & std::ios::showpos) // if showpos is set, we print the implicit '+' sign
				{
					ostr.put('+');
					if (ostr.width() > 0) ostr.width(ostr.width() - 1); // dec width if was specified
				}
				print_positive(ostr, std::forward<T>(val)); // then print the positive value
			}

			return ostr;
		}

		template<typename T, std::enable_if_t<detail::is_biggerints_type<T>::value, int> = 0>
		std::istream &parse_positive(std::istream &istr, T &val, bool noskipws = false)
		{
			constexpr std::uint64_t bits = detail::bit_count<std::decay_t<T>>::value;

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
						if constexpr (!std::is_same_v<std::decay_t<T>, bigint>)
						{
							// detect overflow
							if (val.blocks[bits / 64 - 1] >> (std::uint64_t)(64 - 3 * num_digits)) { istr.setstate(std::ios::failbit); return istr; }
						}

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
						if constexpr (!std::is_same_v<std::decay_t<T>, bigint>)
						{
							// detect overflow
							if (val.blocks[bits / 64 - 1] >> (std::uint64_t)(64 - 4 * num_digits)) { istr.setstate(std::ios::failbit); return istr; }
						}

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

				// perform multiply/add in terms of a larger int type to easily tell if we overflow (only for finite integral types)
				std::conditional_t<std::is_same_v<std::decay_t<T>, bigint>, bigint, double_int<detail::bit_count<std::decay_t<T>>::value * 2, false>> bigval = std::move(val);

				// first char must be a dec digit - don't extract
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

						if constexpr (!std::is_same_v<std::decay_t<T>, bigint>)
						{
							// detect overflow
							if (bigval.blocks[bits / 64]) { istr.setstate(std::ios::failbit); return istr; }
						}
					}

					if (num_digits < 19) break;
				}

				if constexpr (!std::is_same_v<std::decay_t<T>, bigint>) val = (T)bigval; // copy low half of bigval to result
				else val = std::move(bigval);
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

		template<std::uint64_t bits> inline std::istream &operator>>(std::istream &istr, double_int<bits, false> &val) { parse_positive(istr, val); return istr; }
		
		template<typename T, std::enable_if_t<detail::is_signed_double_int<T>::value || std::is_same_v<T, bigint>, int> = 0>
		std::istream &parse_signed(std::istream &istr, T &val)
		{
			int ch;
			bool neg = false; // we'll do signed io in terms of unsigned

			std::istream::sentry sentry(istr);
			if (!sentry) return istr;

			// look at the first char - fail if we can't get one
			if ((ch = istr.peek()) == EOF) { istr.setstate(std::ios::failbit | std::ios::eofbit); return istr; }

			// account for optional sign - extract it if present
			if (ch == '-') { istr.get(); neg = true; }
			else if (ch == '+') istr.get();

			if constexpr (!std::is_same_v<std::decay_t<T>, bigint>)
			{
				// parse the value (and don't skip ws since we already did that)
				double_int<detail::bit_count<std::decay_t<T>>::value, false> tmp;
				if (!parse_positive(istr, tmp, true)) return istr;
				val = tmp; // store the tmp parse location back to val
			}
			else
			{
				if (!parse_positive(istr, val, true)) return istr;
			}

			// account for sign in result
			if (neg) detail::make_neg(val);

			if constexpr (!std::is_same_v<std::decay_t<T>, bigint>)
			{
				// check for signed overflow
				if (detail::is_neg(val) != neg) { istr.setstate(std::ios::failbit); return istr; }
			}

			return istr;
		}

		template<std::uint64_t bits>
		std::istream &operator>>(std::istream &istr, double_int<bits, true> &val) { parse_signed(istr, val); return istr; }
		inline std::istream &operator>>(std::istream &istr, bigint &val) { parse_signed(istr, val); return istr; }

		// -- misc functions -- //

		inline bigint _pow(bigint a, const bigint &b) // pass by value is intentional
		{
			if (detail::is_neg(b)) return {}; // if exponent is negative just return 0

			bigint res = 1;
			std::uint64_t high_bit = highest_set_bit(b);

			for (std::uint64_t bit = 0; bit <= high_bit+1; ++bit)
			{
				if (detail::bit_test(b, bit)) res *= a;
				a *= a;
			}

			return res;
		}
	}
}

// -- std info definitions -- //

namespace std
{
	// -- bigint info -- //

	template<> struct is_integral<BiggerInts::detail::bigint> : std::true_type {};

	template<> struct is_signed<BiggerInts::detail::bigint> : std::true_type {};
	template<> struct is_unsigned<BiggerInts::detail::bigint> : std::true_type {};

	template<> struct make_signed<BiggerInts::detail::bigint> { typedef BiggerInts::detail::bigint type; };
	template<> struct make_unsigned<BiggerInts::detail::bigint> { typedef BiggerInts::detail::bigint type; };

	template<> struct numeric_limits<BiggerInts::detail::bigint>
	{
		static constexpr bool is_specialized = true;
		static constexpr bool is_signed = true;
		static constexpr bool is_integer = true;
		static constexpr bool is_exact = true;
		static constexpr bool has_infinity = false;
		static constexpr bool has_quiet_NaN = false;
		static constexpr bool has_signaling_NaN = false;
		static constexpr bool has_denorm = false;
		static constexpr bool has_denorm_loss = false;
		static constexpr std::float_round_style round_style = std::round_toward_zero;
		static constexpr bool is_iec559 = false;
		static constexpr bool is_bounded = false;
		static constexpr bool is_modulo = false;
		static constexpr int digits = 0;
		static constexpr int digits10 = (int)(digits * 0.301029995663981195213738894724493026768189881462108541310); // log10(2)
		static constexpr int max_digits10 = 0;
		static constexpr int radix = 2;
		static constexpr int min_exponent = 0;
		static constexpr int min_exponent10 = 0;
		static constexpr int max_exponent = 0;
		static constexpr int max_exponent10 = 0;
		static constexpr bool traps = true;
		static constexpr bool tinyness_before = false;

		static BiggerInts::detail::bigint min() { return {}; }
		static BiggerInts::detail::bigint lowest() { return {}; }
		static BiggerInts::detail::bigint max() { return {}; }
		static BiggerInts::detail::bigint epsilon() { return {}; }
		static BiggerInts::detail::bigint round_error() { return {}; }
		static BiggerInts::detail::bigint infinity() { return {}; }
		static BiggerInts::detail::bigint quiet_NaN() { return {}; }
		static BiggerInts::detail::bigint signaling_NaN() { return {}; }
		static BiggerInts::detail::bigint denorm_min() { return {}; }
	};

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
