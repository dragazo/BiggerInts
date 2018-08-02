#ifndef BIGGER_INTS_H
#define BIGGER_INTS_H

#include <cstdint>
#include <type_traits>

namespace BiggerInts
{
	typedef std::uint64_t u64;

	// returns true if val is a power of 2
	inline constexpr bool is_pow2(u64 val) { return val != 0 && (val & (val - 1)) == 0; }
	
	// given a size in bits, returns a power of 2 (also in bits) large enough to contain it and that is no smaller than 8
	inline constexpr u64 round_bits_up(u64 size)
	{
		if (size < 8) return 8;
		while (!is_pow2(size)) size = (size | (size - 1)) + 1;
		return size;
	}

	// given a desired size in bits, provides a typedef that will suffice.
	// 8, 16, 32, and 64 yield std::uint??_t directly. other values are wrappers defined here.
	template<u64 bits> struct uint;
	// identical to uint<bits>::type
	template<u64 bits> using uint_t = typename uint<bits>::type;
	
	// gets the number of bits in the specified type (only works for types accessible via uint_t)
	// bit_count<uint_t<bits>>::value == bits is always true
	template<typename T> struct bit_count;

	// ----------------------------------------- //

	// -- DON'T USE ANYTHING BELOW THIS POINT -- //

	// ----------------------------------------- //

	// holds a type T where assignments to it are masked
	template<typename T, u64 bits> struct masked_single_int;

	// wraps larger integer sizes in terms of 2 smaller ones
	template<u64 bits> struct double_int;

	// -- #hax typedef omega synergy highjinks - don't even need to be defined, boi -- //

	template<u64 bits, std::enable_if_t<(bits > 64 && is_pow2(bits))>* = 0> double_int<bits> returns_proper_type();
	template<u64 bits, std::enable_if_t<(bits > 64 && !is_pow2(bits))>* = 0> masked_single_int<double_int<round_bits_up(bits)>, bits> returns_proper_type();

	template<u64 bits, std::enable_if_t<(bits == 64)>* = 0> std::uint64_t returns_proper_type();
	template<u64 bits, std::enable_if_t<(bits > 32 && bits < 64)>* = 0> masked_single_int<std::uint64_t, bits> returns_proper_type();

	template<u64 bits, std::enable_if_t<(bits == 32)>* = 0> std::uint32_t returns_proper_type();
	template<u64 bits, std::enable_if_t<(bits > 16 && bits < 32)>* = 0> masked_single_int<std::uint32_t, bits> returns_proper_type();

	template<u64 bits, std::enable_if_t<(bits == 16)>* = 0> std::uint16_t returns_proper_type();
	template<u64 bits, std::enable_if_t<(bits > 8 && bits < 16)>* = 0> masked_single_int<std::uint16_t, bits> returns_proper_type();

	template<u64 bits, std::enable_if_t<(bits == 8)>* = 0> std::uint8_t returns_proper_type();
	template<u64 bits, std::enable_if_t<(bits < 8)>* = 0> masked_single_int<std::uint8_t, bits> returns_proper_type();

	template<u64 bits> struct uint { typedef decltype(returns_proper_type<bits>()) type; };

	// -- sigh, back to normal templating :( -- //
	
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
	template<u64 bits> struct bit_count<double_int<bits>> : std::integral_constant<u64, bits> {};

	// -- container iml -- //

	template<typename T, u64 bits> struct masked_single_int
	{
	private:

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

	template<u64 bits> struct double_int
	{
	public: // -- ctor/asgn -- //

		inline constexpr double_int() = default;

		inline constexpr double_int(const double_int&) = default;
		inline constexpr double_int &operator=(const double_int&) = default;

	private: // -- data -- //
		
		static_assert(is_pow2(bits) && bits >= 128, "double_int: bits was out of valid domain");

		typedef uint_t<bits / 2> half;

		half low, high;

	public: // -- conversion -- //

		template<typename U, std::enable_if_t<(bit_count<U>::value <= bits / 2)>* = 0>
		inline constexpr double_int(const U &val) noexcept : low(val), high(0) {}

		template<typename U, std::enable_if_t<(bit_count<U>::value <= bits / 2)>* = 0>
		inline constexpr double_int &operator=(const U &val) noexcept { low = val; high = 0; return *this; }

		template<typename U, std::enable_if_t<(bit_count<U>::value > bits / 2)>* = 0>
		inline constexpr double_int(const U &val) noexcept : low(val), high(val >> (bits / 2)) {}

		template<typename U, std::enable_if_t<(bit_count<U>::value > bits / 2)>* = 0>
		inline constexpr double_int &operator=(const U &val) noexcept { low = val; high = val >> (bits / 2); }

	public: // -- operators -- //

		inline constexpr double_int &operator&=(const double_int &other) noexcept { low &= other.low; high &= other.high; return *this; }
		inline constexpr double_int &operator|=(const double_int &other) noexcept { low |= other.low; high |= other.high; return *this; }
		inline constexpr double_int &operator^=(const double_int &other) noexcept { low ^= other.low; high ^= other.high; return *this; }

		inline constexpr double_int &operator<<=(u64 count) noexcept
		{
			// count masked to (bits / 2) limit - no-op on zero
			if (count &= (bits / 2) - 1)
			{
				high <<= count;
				high |= low >> (bits / 2 - count);
				low <<= count;
			}
			return *this;
		}
		inline constexpr double_int &operator>>=(u64 count) noexcept
		{
			// count masked to (bits / 2) limit - no-op on zero
			if (count &= (bits / 2) - 1)
			{
				low >>= count;
				low |= high << (bits / 2 - count);
				high >>= count;
			}
			return *this;
		}

		inline constexpr friend double_int operator&(const double_int &a, const double_int &b) noexcept { double_int res = a; res &= b; return res; }
		inline constexpr friend double_int operator|(const double_int &a, const double_int &b) noexcept { double_int res = a; res |= b; return res; }
		inline constexpr friend double_int operator^(const double_int &a, const double_int &b) noexcept { double_int res = a; res ^= b; return res; }

		inline constexpr friend double_int operator<<(const double_int &a, u64 count) noexcept { double_int res = a; res <<= b; return res; }
		inline constexpr friend double_int operator>>(const double_int &a, u64 count) noexcept { double_int res = a; res >>= b; return res; }

		// equivalent to <=> but works for any version of C++
		inline constexpr friend int cmp(const double_int &a, const double_int &b) noexcept { return a.high < b.high ? -1 : a.high > b.high ? 1 : a.low < b.low ? -1 : a.low > b.low ? 1 : 0; }
		inline constexpr friend bool operator==(const double_int &a, const double_int &b) noexcept { return cmp(a, b) == 0; }
		inline constexpr friend bool operator!=(const double_int &a, const double_int &b) noexcept { return cmp(a, b) != 0; }
		inline constexpr friend bool operator<(const double_int &a, const double_int &b) noexcept { return cmp(a, b) < 0; }
		inline constexpr friend bool operator<=(const double_int &a, const double_int &b) noexcept { return cmp(a, b) <= 0; }
		inline constexpr friend bool operator>(const double_int &a, const double_int &b) noexcept { return cmp(a, b) > 0; }
		inline constexpr friend bool operator>=(const double_int &a, const double_int &b) noexcept { return cmp(a, b) >= 0; }
	};
}

#endif
