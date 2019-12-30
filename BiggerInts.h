#ifndef DRAGAZO_BIGGER_INTS_H
#define DRAGAZO_BIGGER_INTS_H

#include <cstdint>
#include <climits>
#include <type_traits>
#include <utility>
#include <iostream>
#include <string>
#include <string_view>
#include <stdexcept>
#include <vector>

// detect arch bits
#if _WIN64 || __x86_64__ || __ppc64__
#define DRAGAZO_BIGGER_INTS_ARCHBIT 64
#else
#define DRAGAZO_BIGGER_INTS_ARCHBIT 32
#endif

// 64-bit msvc needs intrin.h for the u128 intrinsics
#if DRAGAZO_BIGGER_INTS_ARCHBIT == 64 && _MSC_VER
#include <intrin.h>
#endif

/*
adds larger integer types:

the template aliases uint_t<N> and int_t<N> represent fixed-size N-bit integer types of any number of bits which is a power of 2 (e.g. 32, 64, 128, etc.).
the signed form int_t<N> is stored in 2's complement.
these are (non-dynamic) value types and all their core operations are constexpr.

the type alias bigint is an arbitrary-sized (signed) integer type stored in 2's complement.

all of these types expose a member called blocks, which is the internal representation of the value (little-endian u64 array).
modifying this directly for uint_t<N> or int_t<N> is relatively safe, but bigint has several data contracts that must be maintained, so you're probably better off not touching it directly.

report any bugs to https://github.com/dragazo/BiggerInts/issues
*/

namespace BiggerInts
{
	namespace detail
	{
		template<std::uint64_t bits, bool sign> struct fixed_int;
		class bigint;
		
		// -- fixed-sized int type selection -- //

		template<std::uint64_t bits, bool sign> struct fixed_int_selector { typedef detail::fixed_int<bits, sign> type; };
		template<bool sign> struct fixed_int_selector<64, sign> { typedef std::conditional_t<sign, std::int64_t, std::uint64_t> type; };
		template<bool sign> struct fixed_int_selector<32, sign> { typedef std::conditional_t<sign, std::int32_t, std::uint32_t> type; };
		template<bool sign> struct fixed_int_selector<16, sign> { typedef std::conditional_t<sign, std::int16_t, std::uint16_t> type; };
		template<bool sign> struct fixed_int_selector<8, sign> { typedef std::conditional_t<sign, std::int8_t, std::uint8_t> type; };

		// -- fixed-sized wrapper types -- //

		struct const_fixed_int_wrapper
		{
			const std::uint64_t *blocks;   // pointer to the blocks array
			std::size_t          blocks_n; // the number of blocks
		};
		struct fixed_int_wrapper
		{
			std::uint64_t *blocks;   // pointer to the blocks array
			std::size_t    blocks_n; // the number of blocks
			constexpr inline operator const_fixed_int_wrapper() noexcept { return { blocks, blocks_n }; }
		};

		// -- core fixed-sizes util -- //

		constexpr bool is_pow2(std::uint64_t val) noexcept { return val != 0 && (val & (val - 1)) == 0; }

		constexpr bool bit_test_unchecked(const_fixed_int_wrapper val, std::uint64_t bit) noexcept { return (val.blocks[bit / 64] >> (bit % 64)) & 1; }
		constexpr void bit_set_unchecked(fixed_int_wrapper val, std::uint64_t bit) noexcept { val.blocks[bit / 64] |= (std::uint64_t)1 << (bit % 64); }

		constexpr bool is_neg(const_fixed_int_wrapper val) noexcept { return val.blocks[val.blocks_n - 1] & 0x8000000000000000; }

		constexpr bool nonzero(const_fixed_int_wrapper val) noexcept
		{
			for (std::size_t i = 0; i < val.blocks_n; ++i) if (val.blocks[i]) return true;
			return false;
		}

		constexpr std::size_t highest_set_bit(std::uint64_t val)
		{
			std::size_t i = 32, res = 0;
			for (; i > 0; i >>= 1)
			{
				std::uint64_t v = val >> i;
				if (v) { val = v; res += i; }
			}
			return res;
		}
		constexpr std::size_t highest_set_bit(const_fixed_int_wrapper val) noexcept
		{
			for (std::size_t i = val.blocks_n; i-- > 0; ) if (val.blocks[i]) return i * 64 + highest_set_bit(val.blocks[i]);
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
		// as _mul_u64() except is allowed to use faster, platform-specific intrinsics
		inline std::pair<std::uint64_t, std::uint64_t> _mul_u64_fast(std::uint64_t a, std::uint64_t b) noexcept
		{
			#if DRAGAZO_BIGGER_INTS_ARCHBIT == 64
				#if _MSC_VER
				a = _umul128(a, b, &b);
				return { a, b };
				#elif defined(__GNUC__)
				auto r = (__uint128_t)a * b;
				return { (std::uint64_t)r, r >> 64 };
				#else
				return _mul_u64(a, b);
				#endif
			#else
			return _mul_u64(a, b);
			#endif
		}

		// -- fixed-size op impl -- //
		
		constexpr void increment(detail::fixed_int_wrapper a) noexcept
		{
			for (std::size_t i = 0; i < a.blocks_n; ++i) if (++a.blocks[i]) break;
		}
		constexpr void decrement(detail::fixed_int_wrapper a) noexcept
		{
			for (std::size_t i = 0; i < a.blocks_n; ++i) if (a.blocks[i]--) break;
		}

		constexpr void make_not(detail::fixed_int_wrapper a) noexcept
		{
			for (std::size_t i = 0; i < a.blocks_n; ++i) a.blocks[i] = ~a.blocks[i];
		}
		constexpr void make_neg(detail::fixed_int_wrapper a) noexcept { detail::make_not(a); detail::increment(a); }

		constexpr void add_same_size(detail::fixed_int_wrapper a, detail::const_fixed_int_wrapper b) noexcept
		{
			std::uint64_t carry = 0;
			for (std::size_t i = 0; i < a.blocks_n; ++i)
			{
				std::uint64_t v = b.blocks[i] + carry;
				carry = (a.blocks[i] += v) < v || v < carry ? 1 : 0;
			}
		}
		constexpr std::uint64_t add_u64(detail::fixed_int_wrapper a, std::uint64_t b) noexcept // returns overflow from high block (0 or 1)
		{
			if ((a.blocks[0] += b) < b)
			{
				std::size_t i = 1;
				for (; i < a.blocks_n; ++i) if (++a.blocks[i]) break; // on overflow carry a 1 into higher blocks
				if (i == a.blocks_n) return 1;                        // if we exhausted all loop iterations, we overflow from the high block
			}
			return 0;
		}

		constexpr void sub_same_size(detail::fixed_int_wrapper a, detail::const_fixed_int_wrapper b) noexcept
		{
			std::uint64_t carry = 1;
			for (std::size_t i = 0; i < a.blocks_n; ++i)
			{
				std::uint64_t v = ~b.blocks[i] + carry;
				carry = (a.blocks[i] += v) < v || v < carry ? 1 : 0;
			}
		}
		constexpr void and_same_size(detail::fixed_int_wrapper a, detail::const_fixed_int_wrapper b) noexcept
		{
			for (std::size_t i = 0; i < a.blocks_n; ++i) a.blocks[i] &= b.blocks[i];
		}
		constexpr void or_same_size(detail::fixed_int_wrapper a, detail::const_fixed_int_wrapper b) noexcept
		{
			for (std::size_t i = 0; i < a.blocks_n; ++i) a.blocks[i] |= b.blocks[i];
		}
		constexpr void xor_same_size(detail::fixed_int_wrapper a, detail::const_fixed_int_wrapper b) noexcept
		{
			for (std::size_t i = 0; i < a.blocks_n; ++i) a.blocks[i] ^= b.blocks[i];
		}

		constexpr void shl(detail::fixed_int_wrapper a, std::uint64_t count) noexcept
		{
			count &= a.blocks_n * 64 - 1;
			const std::size_t full = count / 64;
			if (full)
			{
				for (std::size_t i = a.blocks_n; i-- > full; ) a.blocks[i] = a.blocks[i - full];
				for (std::size_t i = full; i-- > 0; ) a.blocks[i] = 0;
				count &= 63;
			}
			if (count)
			{
				for (std::size_t i = a.blocks_n - 1; i > full; --i) a.blocks[i] = (a.blocks[i] << count) | (a.blocks[i - 1] >> (64 - count));
				a.blocks[full] <<= count;
			}
		}
		constexpr void shr(detail::fixed_int_wrapper a, std::uint64_t count) noexcept
		{
			count &= a.blocks_n * 64 - 1;
			const std::size_t full = count / 64;
			if (full)
			{
				for (std::size_t i = 0; i < a.blocks_n - full; ++i) a.blocks[i] = a.blocks[i + full];
				for (std::size_t i = a.blocks_n - full; i < a.blocks_n; ++i) a.blocks[i] = 0;
				count &= 63;
			}
			if (count)
			{
				for (std::size_t i = full; i < a.blocks_n - 1; ++i) a.blocks[i] = (a.blocks[i] >> count) | (a.blocks[i + 1] << (64 - count));
				a.blocks[a.blocks_n - 1] >>= count;
			}
		}
		constexpr void sar(detail::fixed_int_wrapper a, std::uint64_t count) noexcept
		{
			count &= a.blocks_n * 64 - 1;
			const std::uint64_t fill = (a.blocks[a.blocks_n - 1] & 0x8000000000000000) ? -1 : 0;
			const std::size_t full = count / 64;
			if (full)
			{
				for (std::size_t i = 0; i < a.blocks_n - full; ++i) a.blocks[i] = a.blocks[i + full];
				for (std::size_t i = a.blocks_n - full; i < a.blocks_n; ++i) a.blocks[i] = fill;
				count &= 63;
			}
			if (count)
			{
				for (std::size_t i = full; i < a.blocks_n - 1; ++i) a.blocks[i] = (a.blocks[i] >> count) | (a.blocks[i + 1] << (64 - count));
				a.blocks[a.blocks_n - 1] = (a.blocks[a.blocks_n - 1] >> count) | (fill << (64 - count));
			}
		}

		// res must already be initialized to zero before calling this
		inline void multiply_same_size_already_zero(fixed_int_wrapper res, const_fixed_int_wrapper a, const_fixed_int_wrapper b) noexcept
		{
			for (std::size_t i = 0; i < res.blocks_n; ++i) // loop through each block in a
			{
				// here we construct temp_1, which is the value of a->blocks[i] * b, scaled by i words (think long multiplication from grade school)
				for (std::size_t j = i; j < res.blocks_n - 1; ++j)
				{
					auto p = detail::_mul_u64_fast(a.blocks[i], b.blocks[j - i]); // compute the product with this block from b
					if ((res.blocks[j] += p.first) < p.first)
					{
						for (std::size_t k = j + 1; k < res.blocks_n && !++res.blocks[k]; ++k);
					}
					if ((res.blocks[j + 1] += p.second) < p.second)
					{
						for (std::size_t k = j + 2; k < res.blocks_n && !++res.blocks[k]; ++k);
					}
				}
				res.blocks[res.blocks_n - 1] += a.blocks[i] * b.blocks[res.blocks_n - 1 - i]; // last iteration of the loop is special since it overflows the (finite) array
			}
		}
		inline std::uint64_t multiply_u64_already_zero(fixed_int_wrapper res, const_fixed_int_wrapper a, std::uint64_t b) noexcept // returns overflow from high block (unsigned multiply logic)
		{
			std::uint64_t overflow = 0; // overflow block out of the (fixed-sized) value
			for (std::size_t i = 0; i < res.blocks_n - 1; ++i)
			{
				auto p = detail::_mul_u64_fast(a.blocks[i], b); // compute the product with this block and b
				if ((res.blocks[i] += p.first) < p.first)
				{
					std::size_t j = i + 1;
					for (; j < res.blocks_n; ++j) if (++res.blocks[j]) break;
					if (j == res.blocks_n) ++overflow;
				}
				if ((res.blocks[i + 1] += p.second) < p.second)
				{
					std::size_t j = i + 2;
					for (; j < res.blocks_n; ++j) if (++res.blocks[j]) break;
					if (j == res.blocks_n) ++overflow;
				}
			}
			auto p = detail::_mul_u64_fast(a.blocks[res.blocks_n - 1], b); // compute the product with last block and b
			if ((res.blocks[res.blocks_n - 1] += p.first) < p.first) ++overflow;
			return overflow + p.second;
		}
		
		constexpr bool cmp_less_non_negative(const_fixed_int_wrapper a, const_fixed_int_wrapper b) noexcept
		{
			for (std::size_t i = a.blocks_n; i-- > 0; ) if (a.blocks[i] != b.blocks[i]) return a.blocks[i] < b.blocks[i];
			return false;
		}

		// quo, rem must already be initialized to 0, 0
		constexpr void divmod_unchecked_positive_already_zero(fixed_int_wrapper quo, fixed_int_wrapper rem, const_fixed_int_wrapper num, const_fixed_int_wrapper den) noexcept
		{
			std::uint64_t bit = highest_set_bit(num);

			const std::size_t den_highest_bit = highest_set_bit(den);
			std::size_t res_second_highest_bit = 0;
			std::size_t shift_count = 0;

			while (true)
			{
				++shift_count;
				if (res_second_highest_bit + shift_count >= den_highest_bit)
				{
					shl(rem, shift_count);
					res_second_highest_bit += shift_count;
					shift_count = 0;

					if (bit_test_unchecked(num, bit)) ++rem.blocks[0];

					if (!cmp_less_non_negative(rem, den))
					{
						sub_same_size(rem, den);
						res_second_highest_bit = highest_set_bit(rem);
						bit_set_unchecked(quo, bit);
					}
				}
				else if (bit_test_unchecked(num, bit))
				{
					shl(rem, shift_count);
					res_second_highest_bit += shift_count;
					shift_count = 0;
					++rem.blocks[0];
				}

				if (bit-- == 0) break;
			}

			shl(rem, shift_count);
		}
		constexpr void divmod_unchecked_already_zero(fixed_int_wrapper quo, fixed_int_wrapper rem, fixed_int_wrapper num, fixed_int_wrapper den)
		{
			bool neg = false;
			if (detail::is_neg(num)) { detail::make_neg(num); neg = !neg; }
			if (detail::is_neg(den)) { detail::make_neg(den); neg = !neg; }

			divmod_unchecked_positive_already_zero(quo, rem, num, den);

			if (neg)
			{
				detail::make_neg(quo);
				detail::make_neg(rem);
			}
		}

		// -- custom intrinsics -- //

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
		template<std::uint64_t bits, bool sign> struct is_biggerints_type<fixed_int<bits, sign>> : std::true_type {};
		template<> struct is_biggerints_type<bigint> : std::true_type {};

		// checks if T is a signed fixed_int
		template<typename T> struct is_signed_fixed_int : std::false_type {};
		template<std::uint64_t bits> struct is_signed_fixed_int<fixed_int<bits, true>> : std::true_type {};
	}

	template<std::uint64_t bits> using uint_t = typename detail::fixed_int_selector<bits, false>::type;
	template<std::uint64_t bits> using int_t = typename detail::fixed_int_selector<bits, true>::type;
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

		template<std::uint64_t bits, bool sign> struct bit_count<fixed_int<bits, sign>> : std::integral_constant<std::uint64_t, bits> {};
		template<> struct bit_count<bigint> : std::integral_constant<std::uint64_t, 0ull> {}; // we'll just elect to define bigint as having zero bits - it really has infinite, but that's not representable

		// -- container impl -- //

		bool try_parse_unsigned(fixed_int_wrapper res, std::string_view str, int base = 10);
		bool try_parse_signed(fixed_int_wrapper res, std::string_view str, int base = 10);

		// given a desired number of bits (power of 2), simulates it with 2 ints of half that size - sign flag marks signed/unsigned
		template<std::uint64_t bits, bool sign> struct fixed_int
		{
		public: // -- data -- //

			static_assert(detail::is_pow2(bits) && bits >= 128, "fixed_int: bits was out of valid domain");

			std::uint64_t blocks[bits / 64];

		public: // -- core -- //

			constexpr fixed_int() noexcept : blocks{} {};

			constexpr fixed_int(const fixed_int&) noexcept = default;
			constexpr fixed_int &operator=(const fixed_int&) noexcept = default;

			constexpr fixed_int(const fixed_int<bits, !sign> &other) noexcept
			{
				for (std::size_t i = 0; i < bits / 64; ++i) blocks[i] = other.blocks[i];
			}
			constexpr fixed_int &operator=(const fixed_int<bits, !sign> &other) noexcept
			{
				for (std::size_t i = 0; i < bits / 64; ++i) blocks[i] = other.blocks[i];
				return *this;
			}

		public: // -- promotion constructors -- //

			constexpr fixed_int(unsigned short val) noexcept : fixed_int((unsigned long long)val) {}
			constexpr fixed_int(unsigned int val) noexcept : fixed_int((unsigned long long)val) {}
			constexpr fixed_int(unsigned long val) noexcept : fixed_int((unsigned long long)val) {}
			constexpr fixed_int(unsigned long long val) noexcept : blocks{}
			{
				blocks[0] = val;
				for (std::size_t i = 1; i < bits / 64; ++i) blocks[i] = 0;
			}

			constexpr fixed_int(signed short val) noexcept : fixed_int((signed long long)val) {}
			constexpr fixed_int(signed int val) noexcept : fixed_int((signed long long)val) {}
			constexpr fixed_int(signed long val) noexcept : fixed_int((signed long long)val) {}
			constexpr fixed_int(signed long long val) noexcept : blocks{}
			{
				blocks[0] = val;
				std::uint64_t fill = val < 0 ? -1 : 0;
				for (std::size_t i = 1; i < bits / 64; ++i) blocks[i] = fill;
			}

			template<std::uint64_t _bits, bool _sign, std::enable_if_t<(bits > _bits), int> = 0>
					constexpr fixed_int(const fixed_int<_bits, _sign> &other) noexcept : blocks{}
					{
						for (std::size_t i = 0; i < _bits / 64; ++i) blocks[i] = other.blocks[i];
						std::uint64_t fill = 0;
						if constexpr (_sign) fill = detail::is_neg(other) ? -1 : 0;
						for (std::size_t i = _bits / 64; i < bits / 64; ++i) blocks[i] = fill;
					}

		public: // -- promotion assignment -- //

			constexpr fixed_int &operator=(unsigned short val) noexcept { *this = (unsigned long long)val; return *this; }
			constexpr fixed_int &operator=(unsigned int val) noexcept { *this = (unsigned long long)val; return *this; }
			constexpr fixed_int &operator=(unsigned long val) noexcept { *this = (unsigned long long)val; return *this; }
			constexpr fixed_int &operator=(unsigned long long val) noexcept
			{
				blocks[0] = val;
				for (std::size_t i = 1; i < bits / 64; ++i) blocks[i] = 0;
				return *this;
			}

			constexpr fixed_int &operator=(signed short val) noexcept { *this = (signed long long)val; return *this; }
			constexpr fixed_int &operator=(signed int val) noexcept { *this = (signed long long)val; return *this; }
			constexpr fixed_int &operator=(signed long val) noexcept { *this = (signed long long)val; return *this; }
			constexpr fixed_int &operator=(signed long long val) noexcept
			{
				blocks[0] = val;
				std::uint64_t fill = val < 0 ? -1 : 0;
				for (std::size_t i = 1; i < bits / 64; ++i) blocks[i] = fill;
				return *this;
			}

			template<std::uint64_t _bits, bool _sign, std::enable_if_t<(bits > _bits), int> = 0>
					constexpr fixed_int &operator=(const fixed_int<_bits, _sign> &other) noexcept
					{
						for (std::size_t i = 0; i < _bits / 64; ++i) blocks[i] = other.blocks[i];
						std::uint64_t fill = 0;
						if constexpr (_sign) fill = detail::is_neg(other) ? -1 : 0;
						for (std::size_t i = _bits / 64; i < bits / 64; ++i) blocks[i] = fill;
						return *this;
					}

		public: // -- fixed_int demotions -- //

			template<std::uint64_t _bits, bool _sign, std::enable_if_t<(bits < _bits), int> = 0>
				constexpr explicit fixed_int(const fixed_int<_bits, _sign> &other) noexcept : blocks{}
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
			constexpr friend bool operator!(const fixed_int &a) noexcept { return !(bool)a; }

		public: // -- assigny things -- //

			constexpr fixed_int &operator++() noexcept { increment(wrap(*this)); return *this; }
			[[nodiscard]] constexpr fixed_int operator++(int) noexcept { auto cpy = *this; detail::increment(wrap(*this)); return cpy; }

			constexpr fixed_int &operator--() noexcept { decrement(wrap(*this)); return *this; }
			[[nodiscard]] constexpr fixed_int operator--(int) noexcept { auto cpy = *this; detail::decrement(wrap(*this)); return cpy; }

			constexpr fixed_int &operator+=(const fixed_int &other) noexcept { detail::add_same_size(wrap(*this), wrap(other)); return *this; }
			constexpr fixed_int &operator-=(const fixed_int &other) noexcept { detail::sub_same_size(wrap(*this), wrap(other)); return *this; }

			constexpr fixed_int &operator&=(const fixed_int &other) noexcept { detail::and_same_size(wrap(*this), wrap(other)); return *this; }
			constexpr fixed_int &operator|=(const fixed_int &other) noexcept { detail::or_same_size(wrap(*this), wrap(other)); return *this; }
			constexpr fixed_int &operator^=(const fixed_int &other) noexcept { detail::xor_same_size(wrap(*this), wrap(other)); return *this; }

			constexpr fixed_int &operator<<=(std::uint64_t count) noexcept { detail::shl(wrap(*this), count); return *this; }
			constexpr fixed_int &operator>>=(std::uint64_t count) noexcept
			{
				if constexpr (sign) sar(wrap(*this), count); else shr(wrap(*this), count);
				return *this;
			}

			constexpr fixed_int &operator*=(const fixed_int &other) noexcept { *this = *this * other; return *this; }

			constexpr fixed_int &operator/=(const fixed_int &other) noexcept { *this = divmod(*this, other).first; return *this; }
			constexpr fixed_int &operator%=(const fixed_int &other) noexcept { *this = divmod(*this, other).second; return *this; }

		public: // -- static utilities -- //

			// attempts to parse the string into an integer value - returns true on successful parse.
			// base must be 10 (dec), 16 (hex), 8 (oct), or 0 to automatically determine base from C-style prefix in str.
			static bool try_parse(fixed_int &res, std::string_view str, int base = 10) { if constexpr (sign) return detail::try_parse_signed(wrap(res), str, base); else return detail::try_parse_unsigned(wrap(res), str, base); }
			// as try_parse() but throws std::invalid_argument on failure
			static void parse(fixed_int &res, std::string_view str, int base = 10) { if (!try_parse(res, str, base)) throw std::invalid_argument("failed to parse string"); }
			// as parse() but allocates a new object to hold the result (less efficient if you're just going to assign it to a variable)
			static fixed_int parse(std::string_view str, int base = 10) { fixed_int res; parse(res, str, base); return res; }
		};

		void _collapse(bigint&); // collapses the value in the given bigint to make it satisfy the requirements of the blocks array

		// represents an arbitrarily large (signed) integer type in 2's complement
		class bigint
		{
		public: // -- data -- //

			// the dynamic array of blocks (little-endian). empty array is treated as zero.
			// at all times must be as small as possible while still preserving the correct signage of the 2's complement value.
			// users should probably not touch this directly, but is exposed because it's useful for writing extensions.
			std::vector<std::uint64_t> blocks;

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
			bigint(const fixed_int<bits, sign> &other) { operator=(other); }

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
			bigint &operator=(const fixed_int<bits, sign> &other)
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

		public: // -- assigny things -- //

			bigint &operator++();
			[[nodiscard]] bigint operator++(int);

			bigint &operator--();
			[[nodiscard]] bigint operator--(int);

			bigint &operator+=(const bigint &b);
			bigint &operator-=(const bigint &b);

			bigint &operator&=(const bigint &b);
			bigint &operator|=(const bigint &b);
			bigint &operator^=(const bigint &b);

			bigint &operator<<=(std::uint64_t count);
			bigint &operator>>=(std::uint64_t count);

			bigint &operator*=(const bigint &b);
			bigint &operator*=(bigint &&b);

		public: // -- static utilities -- //

			// attempts to parse the string into an integer value - returns true on successful parse.
			// base must be 10 (dec), 16 (hex), 8 (oct), or 0 to automatically determine base from C-style prefix in str.
			static bool try_parse(bigint &res, std::string_view str, int base = 10);
			// as try_parse() but throws std::invalid_argument on failure
			static void parse(bigint &res, std::string_view str, int base = 10);
			// as parse() but allocates a new object to hold the result (less efficient if you're just going to assign it to a variable)
			static bigint parse(std::string_view str, int base = 10);

			// computes the result of raising a to the power of b
			static bigint pow(bigint a, const bigint &b);
			// computes v!
			static bigint factorial(std::uint64_t v);
		};

		// -- wrapping utils -- //

		template<std::uint64_t bits, bool sign> constexpr fixed_int_wrapper wrap(fixed_int<bits, sign>&&) noexcept = delete;
		template<std::uint64_t bits, bool sign> constexpr fixed_int_wrapper wrap(fixed_int<bits, sign> &v) noexcept { return { v.blocks, bits / 64 }; }
		template<std::uint64_t bits, bool sign> constexpr const_fixed_int_wrapper wrap(const fixed_int<bits, sign> &v) noexcept { return { v.blocks, bits / 64 }; }

		// -- bigint utility definitions -- //

		bool is_neg(const bigint &val) noexcept;
		void make_neg(bigint &a);
		void make_not(bigint &a);

		inline bool nonzero(const bigint &val) noexcept { return !val.blocks.empty(); }

		std::size_t highest_set_bit(const bigint &val) noexcept;

		bool bit_test(const bigint &val, std::uint64_t bit) noexcept;
		bool bit_test_in_bounds_nonzero(const bigint &val, std::uint64_t bit) noexcept;
		void bit_set(bigint &val, std::uint64_t bit) noexcept;

		int cmp(const bigint &a, const bigint &b) noexcept;
		bool cmp_less_non_negative(const bigint &a, const bigint &b) noexcept;

		// -- add -- //

		template<std::uint64_t bits, bool sign> constexpr fixed_int<bits, sign> operator+(const fixed_int<bits, sign> &a, const fixed_int<bits, sign> &b) { auto cpy = a; cpy += b; return cpy; }
		template<std::uint64_t bits, bool sign, typename T, std::enable_if_t<detail::is_builtin_int<T>::value, int> = 0> constexpr fixed_int<bits, sign> operator+(const fixed_int<bits, sign> &a, T b) { return a + (fixed_int<bits, sign>)b; }
		template<std::uint64_t bits, bool sign, typename T, std::enable_if_t<detail::is_builtin_int<T>::value, int> = 0> constexpr fixed_int<bits, sign> operator+(T a, const fixed_int<bits, sign> &b) { return (fixed_int<bits, sign>)a + b; }

		bigint operator+(const bigint &a, const bigint &b);
		bigint operator+(bigint &&a, const bigint &b);
		bigint operator+(const bigint &a, bigint &&b);
		bigint operator+(bigint &&a, bigint &&b);

		// -- sub -- //

		template<std::uint64_t bits, bool sign> constexpr fixed_int<bits, sign> operator-(const fixed_int<bits, sign> &a, const fixed_int<bits, sign> &b) { auto cpy = a; cpy -= b; return cpy; }
		template<std::uint64_t bits, bool sign, typename T, std::enable_if_t<detail::is_builtin_int<T>::value, int> = 0> constexpr fixed_int<bits, sign> operator-(const fixed_int<bits, sign> &a, T b) { return a - (fixed_int<bits, sign>)b; }
		template<std::uint64_t bits, bool sign, typename T, std::enable_if_t<detail::is_builtin_int<T>::value, int> = 0> constexpr fixed_int<bits, sign> operator-(T a, const fixed_int<bits, sign> &b) { return (fixed_int<bits, sign>)a - b; }

		bigint operator-(const bigint &a, const bigint &b);
		bigint operator-(bigint &&a, const bigint &b);
		bigint operator-(const bigint &a, bigint &&b);
		bigint operator-(bigint &&a, bigint &&b);

		// -- and -- //

		template<std::uint64_t bits, bool sign> constexpr fixed_int<bits, sign> operator&(const fixed_int<bits, sign> &a, const fixed_int<bits, sign> &b) { auto cpy = a; cpy &= b; return cpy; }
		template<std::uint64_t bits, bool sign, typename T, std::enable_if_t<detail::is_builtin_int<T>::value, int> = 0> constexpr fixed_int<bits, sign> operator&(const fixed_int<bits, sign> &a, T b) { return a & (fixed_int<bits, sign>)b; }
		template<std::uint64_t bits, bool sign, typename T, std::enable_if_t<detail::is_builtin_int<T>::value, int> = 0> constexpr fixed_int<bits, sign> operator&(T a, const fixed_int<bits, sign> &b) { return (fixed_int<bits, sign>)a & b; }

		bigint operator&(const bigint &a, const bigint &b);
		bigint operator&(bigint &&a, const bigint &b);
		bigint operator&(const bigint &a, bigint &&b);
		bigint operator&(bigint &&a, bigint &&b);

		// -- or -- //

		template<std::uint64_t bits, bool sign> constexpr fixed_int<bits, sign> operator|(const fixed_int<bits, sign> &a, const fixed_int<bits, sign> &b) { auto cpy = a; cpy |= b; return cpy; }
		template<std::uint64_t bits, bool sign, typename T, std::enable_if_t<detail::is_builtin_int<T>::value, int> = 0> constexpr fixed_int<bits, sign> operator|(const fixed_int<bits, sign> &a, T b) { return a | (fixed_int<bits, sign>)b; }
		template<std::uint64_t bits, bool sign, typename T, std::enable_if_t<detail::is_builtin_int<T>::value, int> = 0> constexpr fixed_int<bits, sign> operator|(T a, const fixed_int<bits, sign> &b) { return (fixed_int<bits, sign>)a | b; }

		bigint operator|(const bigint &a, const bigint &b);
		bigint operator|(bigint &&a, const bigint &b);
		bigint operator|(const bigint &a, bigint &&b);
		bigint operator|(bigint &&a, bigint &&b);

		// -- xor -- //

		template<std::uint64_t bits, bool sign> constexpr fixed_int<bits, sign> operator^(const fixed_int<bits, sign> &a, const fixed_int<bits, sign> &b) { auto cpy = a; cpy ^= b; return cpy; }
		template<std::uint64_t bits, bool sign, typename T, std::enable_if_t<detail::is_builtin_int<T>::value, int> = 0> constexpr fixed_int<bits, sign> operator^(const fixed_int<bits, sign> &a, T b) { return a ^ (fixed_int<bits, sign>)b; }
		template<std::uint64_t bits, bool sign, typename T, std::enable_if_t<detail::is_builtin_int<T>::value, int> = 0> constexpr fixed_int<bits, sign> operator^(T a, const fixed_int<bits, sign> &b) { return (fixed_int<bits, sign>)a ^ b; }

		bigint operator^(const bigint &a, const bigint &b);
		bigint operator^(bigint &&a, const bigint &b);
		bigint operator^(const bigint &a, bigint &&b);
		bigint operator^(bigint &&a, bigint &&b);

		// -- shl/sal -- //

		template<std::uint64_t bits, bool sign>
		constexpr fixed_int<bits, sign> operator<<(const fixed_int<bits, sign> &a, std::uint64_t count) noexcept { auto cpy = a; cpy <<= count; return cpy; }

		bigint operator<<(const bigint &val, std::uint64_t count);
		bigint operator<<(bigint &&val, std::uint64_t count);

		// -- shr/sar -- //

		template<std::uint64_t bits, bool sign>
		constexpr fixed_int<bits, sign> operator>>(const fixed_int<bits, sign> &a, std::uint64_t count) noexcept { auto cpy = a; cpy >>= count; return cpy; }

		bigint operator>>(const bigint &val, std::uint64_t count);
		bigint operator>>(bigint &&val, std::uint64_t count);

		// -- unary -- //

		template<std::uint64_t bits, bool sign>
		constexpr fixed_int<bits, sign> operator+(const fixed_int<bits, sign> &a) noexcept { return a; }

		bigint operator+(const bigint &a);
		bigint operator+(bigint &&a);

		template<std::uint64_t bits>
		constexpr fixed_int<bits, true> operator-(const fixed_int<bits, true> &a) noexcept { fixed_int<bits, true> res = a; detail::make_neg(wrap(res)); return res; }

		bigint operator-(const bigint &a);
		bigint operator-(bigint &&a);

		template<std::uint64_t bits, bool sign>
		constexpr fixed_int<bits, sign> operator~(const fixed_int<bits, sign> &a) noexcept { fixed_int<bits, sign> res = a; detail::make_not(wrap(res)); return res; }

		bigint operator~(const bigint &a);
		bigint operator~(bigint &&a);

		// -- mul -- //

		template<std::uint64_t bits, bool sign>
		constexpr fixed_int<bits, sign> operator*(const fixed_int<bits, sign> &a, const fixed_int<bits, sign> &b) noexcept { fixed_int<bits, sign> res = 0; detail::multiply_same_size_already_zero(wrap(res), wrap(a), wrap(b)); return res; }

		bigint operator*(const bigint &a, const bigint &b);
		bigint operator*(bigint &&a, const bigint &b);
		bigint operator*(const bigint &a, bigint &&b);
		bigint operator*(bigint &&a, bigint &&b);

		// -- comparison -- //

		inline bool operator==(const bigint &a, const bigint &b) noexcept { return a.blocks == b.blocks; }
		inline bool operator!=(const bigint &a, const bigint &b) noexcept { return a.blocks != b.blocks; }
		inline bool operator<(const bigint &a, const bigint &b) noexcept { return cmp(a, b) < 0; }
		inline bool operator<=(const bigint &a, const bigint &b) noexcept { return cmp(a, b) <= 0; }
		inline bool operator>(const bigint &a, const bigint &b) noexcept { return cmp(a, b) > 0; }
		inline bool operator>=(const bigint &a, const bigint &b) noexcept { return cmp(a, b) >= 0; }

		// -- divmod -- //

		void divmod_unchecked_positive_already_zero(bigint &quo, bigint &rem, const bigint &num, const bigint &den);
		std::pair<bigint, bigint> divmod_unchecked_positive(const bigint &num, const bigint &den);

		std::pair<bigint, bigint> divmod(const bigint &a, const bigint &b);
		std::pair<bigint, bigint> divmod(bigint &&a, const bigint &b);
		std::pair<bigint, bigint> divmod(const bigint &a, bigint &&b);
		std::pair<bigint, bigint> divmod(bigint &&a, bigint &&b);

		bigint operator/(const bigint &num, const bigint &den);
		bigint operator/(bigint &&num, const bigint &den);
		bigint operator/(const bigint &num, bigint &&den);
		bigint operator/(bigint &&num, bigint &&den);

		bigint operator%(const bigint &num, const bigint &den);
		bigint operator%(bigint &&num, const bigint &den);
		bigint operator%(const bigint &num, bigint &&den);
		bigint operator%(bigint &&num, bigint &&den);

		bigint &operator/=(bigint &num, const bigint &den);
		bigint &operator/=(bigint &num, bigint &&den);

		bigint &operator%=(bigint &num, const bigint &den);
		bigint &operator%=(bigint &num, bigint &&den);

		template<std::uint64_t bits>
		constexpr std::pair<fixed_int<bits, false>, fixed_int<bits, false>> divmod(const fixed_int<bits, false> &num, const fixed_int<bits, false> &den)
		{
			if (!detail::nonzero(wrap(den))) throw std::domain_error("divide by zero");
			std::pair<fixed_int<bits, false>, fixed_int<bits, false>> res(0ull, 0ull); // must initialize these to zero
			detail::divmod_unchecked_positive_already_zero(wrap(res.first), wrap(res.second), wrap(num), wrap(den));
			return res;
		}
		template<std::uint64_t bits>
		constexpr std::pair<fixed_int<bits, true>, fixed_int<bits, true>> divmod(const fixed_int<bits, true> &num, const fixed_int<bits, true> &den)
		{
			if (!detail::nonzero(wrap(den))) throw std::domain_error("divide by zero");
			std::pair<fixed_int<bits, true>, fixed_int<bits, true>> res(0ull, 0ull); // must initialize these to zero
			fixed_int<bits, true> num_cpy = num;
			fixed_int<bits, true> den_cpy = den;
			detail::divmod_unchecked_already_zero(wrap(res.first), wrap(res.second), wrap(num_cpy), wrap(den_cpy));
			return res;
		}

		template<std::uint64_t bits, bool sign>
		constexpr fixed_int<bits, sign> operator/(const fixed_int<bits, sign> &num, const fixed_int<bits, sign> &den) { return divmod(num, den).first; }
		template<std::uint64_t bits, bool sign>
		constexpr fixed_int<bits, sign> operator%(const fixed_int<bits, sign> &num, const fixed_int<bits, sign> &den) { return divmod(num, den).second; }

		// -- cmp -- //

		// cmp() will be equivalent to <=> but works for any version of C++

		constexpr bool cmp_less_non_negative_same_size(const_fixed_int_wrapper a, const_fixed_int_wrapper b) noexcept // as cmp() but requires they both be non-negative
		{
			for (std::size_t i = a.blocks_n; i-- > 0; ) if (a.blocks[i] != b.blocks[i]) return a.blocks[i] < b.blocks[i];
			return false;
		}
		template<std::uint64_t bits, bool sign>
		constexpr bool cmp_less_non_negative(const fixed_int<bits, sign> &a, const fixed_int<bits, sign> &b) noexcept { return cmp_less_non_negative_same_size(wrap(a), wrap(b)); }

		template<bool sign>
		constexpr int cmp(const_fixed_int_wrapper a, const_fixed_int_wrapper b) noexcept
		{
			if (b.blocks_n > a.blocks_n) return -cmp<sign>(b, a); // wolog let a be at least as large as b (in physical size)
			else
			{
				if constexpr (sign)
				{
					const bool a_neg = detail::is_neg(a);
					const bool b_neg = detail::is_neg(b);
					if (a_neg ^ b_neg) return a_neg ? -1 : 1; // if they have different signs we know right away
				}

				if (a.blocks_n > b.blocks_n)
				{
					std::uint64_t fill = 0;
					if constexpr (sign) fill = detail::is_neg(b) ? 0xffffffffffffffffull : 0ull;
					for (std::size_t i = a.blocks_n; i-- > b.blocks_n; ) if (a.blocks[i] != fill) // if a has higher magnitude than b we can tell (immediately for unsigned, or from sign bit for signed)
					{
						if constexpr (!sign) return 1;
						else return detail::is_neg(a) ? -1 : 1;
					}
				}

				for (std::size_t i = b.blocks_n; i-- > 0; ) if (a.blocks[i] != b.blocks[i]) return a.blocks[i] < b.blocks[i] ? -1 : 1; // otherwise base becision on first different block

				return 0; // otherwise they're equal
			}
		}
		template<std::uint64_t bits_1, std::uint64_t bits_2, bool sign>
		constexpr int cmp(const fixed_int<bits_1, sign> &a, const fixed_int<bits_2, sign> &b) noexcept { return cmp<true>(wrap(a), wrap(b)); }

		template<bool sign>
		constexpr int cmp(const_fixed_int_wrapper a, std::conditional_t<sign, long long, unsigned long long> val) noexcept
		{
			std::uint64_t fill = 0;
			if constexpr (sign) fill = val < 0 ? 0xffffffffffffffffull : 0ull;
			for (std::size_t i = a.blocks_n - 1; i > 0; --i) if (a.blocks[i] != fill)
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
		constexpr int cmp(const fixed_int<bits, sign> &a, std::conditional_t<sign, long long, unsigned long long> val) noexcept { return cmp<sign>(wrap(a), val); }
		
		template<std::uint64_t bits1, std::uint64_t bits2, bool sign> constexpr bool operator==(const fixed_int<bits1, sign> &a, const fixed_int<bits2, sign> &b) noexcept { return cmp(a, b) == 0; }
		template<std::uint64_t bits1, std::uint64_t bits2, bool sign> constexpr bool operator!=(const fixed_int<bits1, sign> &a, const fixed_int<bits2, sign> &b) noexcept { return cmp(a, b) != 0; }
		template<std::uint64_t bits1, std::uint64_t bits2, bool sign> constexpr bool operator<(const fixed_int<bits1, sign> &a, const fixed_int<bits2, sign> &b) noexcept { return cmp(a, b) < 0; }
		template<std::uint64_t bits1, std::uint64_t bits2, bool sign> constexpr bool operator<=(const fixed_int<bits1, sign> &a, const fixed_int<bits2, sign> &b) noexcept { return cmp(a, b) <= 0; }
		template<std::uint64_t bits1, std::uint64_t bits2, bool sign> constexpr bool operator>(const fixed_int<bits1, sign> &a, const fixed_int<bits2, sign> &b) noexcept { return cmp(a, b) > 0; }
		template<std::uint64_t bits1, std::uint64_t bits2, bool sign> constexpr bool operator>=(const fixed_int<bits1, sign> &a, const fixed_int<bits2, sign> &b) noexcept { return cmp(a, b) >= 0; }

		template<std::uint64_t bits, bool sign, typename T, std::enable_if_t<detail::is_builtin_int<T>::value && std::is_signed_v<T> == sign, int> = 0> constexpr bool operator==(const fixed_int<bits, sign> &a, T b) noexcept { return cmp(a, b) == 0; }
		template<std::uint64_t bits, bool sign, typename T, std::enable_if_t<detail::is_builtin_int<T>::value && std::is_signed_v<T> == sign, int> = 0> constexpr bool operator!=(const fixed_int<bits, sign> &a, T b) noexcept { return cmp(a, b) != 0; }
		template<std::uint64_t bits, bool sign, typename T, std::enable_if_t<detail::is_builtin_int<T>::value && std::is_signed_v<T> == sign, int> = 0> constexpr bool operator<(const fixed_int<bits, sign> &a, T b) noexcept { return cmp(a, b) < 0; }
		template<std::uint64_t bits, bool sign, typename T, std::enable_if_t<detail::is_builtin_int<T>::value && std::is_signed_v<T> == sign, int> = 0> constexpr bool operator<=(const fixed_int<bits, sign> &a, T b) noexcept { return cmp(a, b) <= 0; }
		template<std::uint64_t bits, bool sign, typename T, std::enable_if_t<detail::is_builtin_int<T>::value && std::is_signed_v<T> == sign, int> = 0> constexpr bool operator>(const fixed_int<bits, sign> &a, T b) noexcept { return cmp(a, b) > 0; }
		template<std::uint64_t bits, bool sign, typename T, std::enable_if_t<detail::is_builtin_int<T>::value && std::is_signed_v<T> == sign, int> = 0> constexpr bool operator>=(const fixed_int<bits, sign> &a, T b) noexcept { return cmp(a, b) >= 0; }

		template<std::uint64_t bits, bool sign, typename T, std::enable_if_t<detail::is_builtin_int<T>::value && std::is_signed_v<T> == sign, int> = 0> constexpr bool operator==(T a, const fixed_int<bits, sign> &b) noexcept { return cmp(b, a) == 0; }
		template<std::uint64_t bits, bool sign, typename T, std::enable_if_t<detail::is_builtin_int<T>::value && std::is_signed_v<T> == sign, int> = 0> constexpr bool operator!=(T a, const fixed_int<bits, sign> &b) noexcept { return cmp(b, a) != 0; }
		template<std::uint64_t bits, bool sign, typename T, std::enable_if_t<detail::is_builtin_int<T>::value && std::is_signed_v<T> == sign, int> = 0> constexpr bool operator<(T a, const fixed_int<bits, sign> &b) noexcept { return cmp(b, a) > 0; }
		template<std::uint64_t bits, bool sign, typename T, std::enable_if_t<detail::is_builtin_int<T>::value && std::is_signed_v<T> == sign, int> = 0> constexpr bool operator<=(T a, const fixed_int<bits, sign> &b) noexcept { return cmp(b, a) >= 0; }
		template<std::uint64_t bits, bool sign, typename T, std::enable_if_t<detail::is_builtin_int<T>::value && std::is_signed_v<T> == sign, int> = 0> constexpr bool operator>(T a, const fixed_int<bits, sign> &b) noexcept { return cmp(b, a) < 0; }
		template<std::uint64_t bits, bool sign, typename T, std::enable_if_t<detail::is_builtin_int<T>::value && std::is_signed_v<T> == sign, int> = 0> constexpr bool operator>=(T a, const fixed_int<bits, sign> &b) noexcept { return cmp(b, a) <= 0; }

		inline bool operator==(const bigint &a, long long val) noexcept
		{
			return (val == 0 && a.blocks.empty()) || (a.blocks.size() == 1 && a.blocks[0] == (unsigned long long)val);
		}
		inline bool operator==(const bigint &a, long val) noexcept { return a == (long long)val; }
		inline bool operator==(const bigint &a, int val) noexcept { return a == (long long)val; }
		inline bool operator==(const bigint &a, short val) noexcept { return a == (long long)val; }

		inline bool operator==(long long val, const bigint &a) noexcept { return a == val; }
		inline bool operator==(long val, const bigint &a) noexcept { return a == (long long)val; }
		inline bool operator==(int val, const bigint &a) noexcept { return a == (long long)val; }
		inline bool operator==(short val, const bigint &a) noexcept { return a == (long long)val; }

		inline bool operator==(const bigint &a, unsigned long long val) noexcept
		{
			return (val == 0 && a.blocks.empty()) || (a.blocks.size() == 1 && a.blocks[0] == val) || (a.blocks.size() == 2 && a.blocks[0] == val && a.blocks[1] == 0);
		}
		inline bool operator==(const bigint &a, unsigned long val) noexcept { return a == (unsigned long long)val; }
		inline bool operator==(const bigint &a, unsigned int val) noexcept { return a == (unsigned long long)val; }
		inline bool operator==(const bigint &a, unsigned short val) noexcept { return a == (unsigned long long)val; }

		inline bool operator==(unsigned long long val, const bigint &a) noexcept { return a == val; }
		inline bool operator==(unsigned long val, const bigint &a) noexcept { return a == (unsigned long long)val; }
		inline bool operator==(unsigned int val, const bigint &a) noexcept { return a == (unsigned long long)val; }
		inline bool operator==(unsigned short val, const bigint &a) noexcept { return a == (unsigned long long)val; }

		inline bool operator!=(const bigint &a, long long val) noexcept { return !(a == val); }
		inline bool operator!=(const bigint &a, long val) noexcept { return !(a == (long long)val); }
		inline bool operator!=(const bigint &a, int val) noexcept { return !(a == (long long)val); }
		inline bool operator!=(const bigint &a, short val) noexcept { return !(a == (long long)val); }

		inline bool operator!=(long long val, const bigint &a) noexcept { return !(a == val); }
		inline bool operator!=(long val, const bigint &a) noexcept { return !(a == (long long)val); }
		inline bool operator!=(int val, const bigint &a) noexcept { return !(a == (long long)val); }
		inline bool operator!=(short val, const bigint &a) noexcept { return !(a == (long long)val); }

		inline bool operator!=(const bigint &a, unsigned long long val) noexcept { return !(a == val); }
		inline bool operator!=(const bigint &a, unsigned long val) noexcept { return !(a == (unsigned long long)val); }
		inline bool operator!=(const bigint &a, unsigned int val) noexcept { return !(a == (unsigned long long)val); }
		inline bool operator!=(const bigint &a, unsigned short val) noexcept { return !(a == (unsigned long long)val); }

		inline bool operator!=(unsigned long long val, const bigint &a) noexcept { return !(a == val); }
		inline bool operator!=(unsigned long val, const bigint &a) noexcept { return !(a == (unsigned long long)val); }
		inline bool operator!=(unsigned int val, const bigint &a) noexcept { return !(a == (unsigned long long)val); }
		inline bool operator!=(unsigned short val, const bigint &a) noexcept { return !(a == (unsigned long long)val); }

		// -- io -- //

		void print_positive_core(std::ostream &ostr, const_fixed_int_wrapper val);
		void print_signed_core(std::ostream &ostr, const_fixed_int_wrapper val);
		
		std::ostream &operator<<(std::ostream &ostr, const bigint &val);
		std::ostream &operator<<(std::ostream &ostr, bigint &&val);

		template<std::uint64_t bits>
		std::ostream &operator<<(std::ostream &ostr, const fixed_int<bits, false> &val) { print_positive_core(ostr, wrap(val)); return ostr; }
		template<std::uint64_t bits>
		std::ostream &operator<<(std::ostream &ostr, const fixed_int<bits, true> &val) { print_signed_core(ostr, wrap(val)); return ostr; }

		void parse_positive_core(std::istream &istr, fixed_int_wrapper val);
		void parse_signed_core(std::istream &istr, fixed_int_wrapper val);

		template<std::uint64_t bits>
		std::istream &operator>>(std::istream &istr, fixed_int<bits, false> &val) { parse_positive_core(istr, wrap(val)); return istr; }
		template<std::uint64_t bits>
		std::istream &operator>>(std::istream &istr, fixed_int<bits, true> &val) { parse_signed_core(istr, wrap(val)); return istr; }

		std::istream &operator>>(std::istream &istr, bigint &val);
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

	// -- fixed_int info -- //

	template<std::uint64_t bits, bool sign> struct is_integral<BiggerInts::detail::fixed_int<bits, sign>> : std::true_type {};

	template<std::uint64_t bits, bool sign> struct is_signed<BiggerInts::detail::fixed_int<bits, sign>> : std::integral_constant<bool, sign> {};
	template<std::uint64_t bits, bool sign> struct is_unsigned<BiggerInts::detail::fixed_int<bits, sign>> : std::integral_constant<bool, !sign> {};

	template<std::uint64_t bits, bool sign> struct make_signed<BiggerInts::detail::fixed_int<bits, sign>> { typedef BiggerInts::detail::fixed_int<bits, true> type; };
	template<std::uint64_t bits, bool sign> struct make_unsigned<BiggerInts::detail::fixed_int<bits, sign>> { typedef BiggerInts::detail::fixed_int<bits, false> type; };
	
	template<std::uint64_t bits, bool sign> struct numeric_limits<BiggerInts::detail::fixed_int<bits, sign>>
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

		static constexpr BiggerInts::detail::fixed_int<bits, sign> _min = sign ? (BiggerInts::detail::fixed_int<bits, sign>)1 << (bits - 1) : (BiggerInts::detail::fixed_int<bits, sign>)0;
		static constexpr BiggerInts::detail::fixed_int<bits, sign> _max = sign ? ~((BiggerInts::detail::fixed_int<bits, sign>)1 << (bits - 1)) : ~(BiggerInts::detail::fixed_int<bits, sign>)0;
		static constexpr BiggerInts::detail::fixed_int<bits, sign> _zero = (BiggerInts::detail::fixed_int<bits, sign>)0;

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
}

#endif
