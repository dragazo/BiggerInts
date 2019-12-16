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
#include <vector>

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
					if (bit_test(num, bit))
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
						bit_set(res.first, bit);
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
					if (bit_test(num, bit)) ++res.second.blocks[0];

					if (res.second >= den)
					{
						res.second -= den;
						bit_set(res.first, bit);
					}

					if (bit-- == 0) break;
				}
			}

			return res;
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

			constexpr double_int() noexcept = default;

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
						return *this;
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
			static double_int parse(std::string_view str, int base = 10)
			{
				double_int res;
				if (!try_parse(res, str, base)) throw std::invalid_argument("failed to parse string");
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
			bigint() noexcept {}

			// copies the value of other
			bigint(const bigint&) = default;
			bigint &operator=(const bigint&) = default;

			// moves the value of other to this object. the moved-from object is zero after this operation.
			bigint(bigint&&) = default;
			bigint &operator=(bigint&&) = default;

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
			}

			bigint &operator=(short val) { return operator=((long long)val); }
			bigint &operator=(int val) { return operator=((long long)val); }
			bigint &operator=(long val) { return operator=((long long)val); }
			bigint &operator=(long long val)
			{
				blocks.clear();
				if (val) blocks.push_back((unsigned long long)val);
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

		inline bigint &operator+=(bigint &a, const bigint &b)
		{
			std::size_t min = std::min(a.blocks.size(), b.blocks.size());
			const bool a_neg = detail::is_neg(a);
			const bool b_neg = detail::is_neg(b);

			// compute addition on the mutually covered range
			std::uint64_t carry = 0;
			for (std::size_t i = 0; i < min; ++i)
			{
				std::uint64_t v = b.blocks[i] + carry;
				carry = (a.blocks[i] += v) < v || v < carry ? 1 : 0;
			}

			// perform addition on the extended range (ghost blocks)
			if (a.blocks.size() < b.blocks.size())
			{
				a.blocks.resize(b.blocks.size(), a_neg ? -1 : 0);
				for (std::size_t i = min; i < a.blocks.size(); ++i)
				{
					std::uint64_t v = b.blocks[i] + carry;
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

			// perform addition on te infinite prefix blocks
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
			if (count >= 64)
			{
				std::size_t full = count / 64;
				for (std::size_t i = bits / 64; i-- > full; ) val.blocks[i] = val.blocks[i - full];
				for (std::size_t i = full; i > 0;) val.blocks[--i] = 0;
				count &= 63;
			}
			if (count)
			{
				for (std::size_t i = bits / 64 - 1; i > 0; --i) val.blocks[i] = (val.blocks[i] << count) | (val.blocks[i - 1] >> (64 - count));
				val.blocks[0] <<= count;
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
			count &= bits - 1;
			if (count >= 64)
			{
				std::size_t full = count / 64;
				for (std::size_t i = 0; i < bits / 64 - full; ++i) val.blocks[i] = val.blocks[i + full];
				for (std::size_t i = bits / 64 - full; i < bits / 64; ++i) val.blocks[i] = 0;
				count &= 63;
			}
			if (count)
			{
				for (std::size_t i = 0; i < bits / 64 - 1; ++i) val.blocks[i] = (val.blocks[i] >> count) | (val.blocks[i + 1] << (64 - count));
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
			if (count >= 64)
			{
				std::size_t full = count / 64;
				for (std::size_t i = 0; i < bits / 64 - full; ++i) val.blocks[i] = val.blocks[i + full];
				for (std::size_t i = bits / 64 - full; i < bits / 64; ++i) val.blocks[i] = fill;
				count &= 63;
			}
			if (count)
			{
				for (std::size_t i = 0; i < bits / 64 - 1; ++i) val.blocks[i] = (val.blocks[i] >> count) | (val.blocks[i + 1] << (64 - count));
				val.blocks[bits / 64 - 1] = (val.blocks[bits / 64 - 1] >> count) | (fill << (64 - count));
			}
			return val;
		}

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

		template<std::uint64_t bits>
		constexpr double_int<bits, false> _multiply(const double_int<bits, false> *a, const double_int<bits, false> *b) noexcept
		{
			std::uint64_t b_high_bit = highest_set_bit(*b);

			if (std::uint64_t a_high_bit = highest_set_bit(*a); a_high_bit < b_high_bit)
			{
				std::swap(a, b);
				b_high_bit = a_high_bit;
			}

			if constexpr (bits >= multiply_shift_count_thresh)
			{
				double_int<bits, false> res = 0u, _a = *a;
				std::size_t shift_count = 0;
				for (std::uint64_t bit = 0; bit <= b_high_bit; ++bit, ++shift_count)
					if (detail::bit_test(*b, bit))
					{
						_a <<= shift_count;
						shift_count = 0;
						res += _a;
					}
				return res;
			}
			else
			{
				double_int<bits, false> res = 0u, _a = *a;
				for (std::uint64_t bit = 0; bit <= b_high_bit; ++bit, _a <<= 1)
					if (detail::bit_test(*b, bit)) res += _a;
				return res;
			}
		}

		template<std::uint64_t bits>
		constexpr double_int<bits, false> operator*(const double_int<bits, false> &a, const double_int<bits, false> &b) noexcept
		{
			return detail::_multiply(&a, &b);
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

		template<std::uint64_t bits_1, std::uint64_t bits_2, bool sign>
		constexpr int cmp(const double_int<bits_1, sign> &a, const double_int<bits_2, sign> &b) noexcept
		{
			if constexpr (bits_2 > bits_1) return -cmp(b, a); // wolog let a be at least as large as b (in physical size)
			else
			{
				if constexpr (bits_1 > bits_2)
				{
					std::uint64_t fill = 0;
					if constexpr (sign) fill = detail::is_neg(b) ? 0xffffffffffffffffull : 0ull;
					for (std::size_t i = bits_1 / 64; i-- > bits_2 / 64; ) if (a.blocks[i] != fill)
					{
						if constexpr (!sign) return 1;
						else return detail::is_neg(a) ? -1 : 1;
					}
				}
				for (std::size_t i = bits_2 / 64; i-- > 0;) if (a.blocks[i] != b.blocks[i])
				{
					if constexpr (!sign) return a.blocks[i] < b.blocks[i] ? -1 : 1;
					else return a.blocks[i] < b.blocks[i] ? 1 : -1;
				}
				return 0;
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
			return val == 0 && a.blocks.empty() || a.blocks.size() == 1 && a.blocks[0] == val;
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
			return val == 0 && a.blocks.empty() || a.blocks[0] == val && (a.blocks.size() == 1 || a.blocks.size() == 2 && a.blocks[1] == 0);
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
					block = cpy.blocks[0] & 0777777777777777777777;
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
					block = cpy.blocks[0];
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
					block = temp.second.blocks[0];
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
						if (val.blocks[bits / 64 - 1] >> (std::uint64_t)(64 - 3 * num_digits)) { istr.setstate(std::ios::failbit); return istr; }

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
						if (val.blocks[bits / 64 - 1] >> (std::uint64_t)(64 - 4 * num_digits)) { istr.setstate(std::ios::failbit); return istr; }

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
						if (bigval.blocks[bits / 64]) { istr.setstate(std::ios::failbit); return istr; }
					}

					if (num_digits < 19) break;
				}

				val = (decltype(val))bigval; // copy low half of bigval to result
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
			if (!parse_udouble_int(istr, tmp, true)) return istr;
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
