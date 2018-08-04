#ifndef BIGGER_INTS_H
#define BIGGER_INTS_H

#include <cstdint>
#include <type_traits>
#include <iostream>
#include <string>
#include <utility>
#include <exception>

namespace BiggerInts
{
	typedef std::uint64_t u64;

	// given a desired size in bits, provides a typedef that will suffice.
	// 8, 16, 32, and 64 yield std::uint??_t directly. other values are wrappers defined here.
	template<u64 bits> struct uint;
	// identical to uint<bits>::type
	template<u64 bits> using uint_t = typename uint<bits>::type;
	
	// gets the number of bits in the specified type (only works for types accessible via uint_t)
	// bit_count<uint_t<bits>>::value == bits is always true
	template<typename T> struct bit_count;

	// given a desired size in bits, creates it from two ints of half that size
	template<u64 bits> inline constexpr uint_t<bits> build_uint(const uint_t<bits / 2> &high, const uint_t<bits / 2> &low) noexcept;

	// ----------------------------------------- //

	// -- DON'T USE ANYTHING BELOW THIS POINT -- //

	// ----------------------------------------- //

	// -- helpers and types -- //

	// returns true if val is a power of 2
	inline constexpr bool is_pow2(u64 val) { return val != 0 && (val & (val - 1)) == 0; }

	// given a size in bits, returns a power of 2 (also in bits) large enough to contain it and that is no smaller than 8
	inline constexpr u64 round_bits_up(u64 size)
	{
		if (size < 8) return 8;
		while (!is_pow2(size)) size = (size | (size - 1)) + 1;
		return size;
	}

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

	// -- builder impl -- //

	template<u64 bits> inline constexpr uint_t<bits> build_uint(const uint_t<bits / 2> &high, const uint_t<bits / 2> &low) noexcept { return {high, low}; }

	template<> inline std::uint64_t build_uint<64>(const std::uint32_t &high, const std::uint32_t &low) noexcept { return ((std::uint64_t)high << 32) | low; }
	template<> inline std::uint32_t build_uint<32>(const std::uint16_t &high, const std::uint16_t &low) noexcept { return ((std::uint32_t)high << 16) | low; }
	template<> inline std::uint16_t build_uint<16>(const std::uint8_t &high, const std::uint8_t &low) noexcept { return ((std::uint16_t)high << 8) | low; }

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

		template<typename U, std::enable_if_t<(bit_count<U>::value <= bits / 2)>* = 0>
		inline constexpr operator U() const noexcept { return low; }

		inline constexpr double_int(half _high, half _low) noexcept : high(_high), low(_low) {}

	public: // -- operators -- //

		inline constexpr double_int &operator++() { if (!++low) ++high; return *this; }
		inline constexpr double_int &operator--() { if (!low) --high; --low; return *this; }

		inline constexpr double_int operator++(int) { double_int val = *this; ++*this; return *this; }
		inline constexpr double_int operator--(int) { double_int val = *this; --*this; return *this; }

		inline constexpr double_int &operator+=(const double_int &other) noexcept
		{
			low += other.low;
			if (low < other.low) ++high;
			high += other.high;
			return *this;
		}
		inline constexpr double_int &operator-=(const double_int &other) noexcept
		{
			if (low < other.low) --high;
			low -= other.low;
			high -= other.high;
			return *this;
		}

		inline constexpr double_int &operator&=(const double_int &other) noexcept { low &= other.low; high &= other.high; return *this; }
		inline constexpr double_int &operator|=(const double_int &other) noexcept { low |= other.low; high |= other.high; return *this; }
		inline constexpr double_int &operator^=(const double_int &other) noexcept { low ^= other.low; high ^= other.high; return *this; }

		inline constexpr double_int &operator<<=(u64 count) noexcept
		{
			// count masked to bits limit - no-op on zero
			if (count &= bits - 1)
			{
				if (count >= bits / 2)
				{
					high = low;
					high <<= count - bits / 2;
					low = 0;
				}
				else
				{
					high <<= count;
					high |= low >> (bits / 2 - count);
					low <<= count;
				}
			}
			return *this;
		}
		inline constexpr double_int &operator>>=(u64 count) noexcept
		{
			// count masked to bits limit - no-op on zero
			if (count &= bits - 1)
			{
				if (count >= bits / 2)
				{
					low = high;
					low >>= count - bits / 2;
					high = 0;
				}
				else
				{
					low >>= count;
					low |= high << (bits / 2 - count);
					high >>= count;
				}
			}
			return *this;
		}

		inline constexpr explicit operator bool() const noexcept { return low || high; }
		inline constexpr friend bool operator!(const double_int &a) noexcept { return !(bool)a; }

		inline constexpr friend double_int operator~(const double_int &a) noexcept { double_int res; res.low = ~a.low; res.high = ~a.high; return res; }

		inline constexpr friend double_int operator+(const double_int &a, const double_int &b) noexcept { double_int res = a; res += b; return res; }
		inline constexpr friend double_int operator-(const double_int &a, const double_int &b) noexcept { double_int res = a; res -= b; return res; }

		inline constexpr friend double_int operator&(const double_int &a, const double_int &b) noexcept { double_int res = a; res &= b; return res; }
		inline constexpr friend double_int operator|(const double_int &a, const double_int &b) noexcept { double_int res = a; res |= b; return res; }
		inline constexpr friend double_int operator^(const double_int &a, const double_int &b) noexcept { double_int res = a; res ^= b; return res; }

		inline constexpr friend double_int operator<<(const double_int &a, u64 count) noexcept { double_int res = a; res <<= count; return res; }
		inline constexpr friend double_int operator>>(const double_int &a, u64 count) noexcept { double_int res = a; res >>= count; return res; }

		inline constexpr friend double_int operator*(double_int a, const double_int &b) noexcept
		{
			double_int res = 0;
			for (u64 bit = 0; a; ++bit, a <<= 1)
				if (bit_test(b, bit)) res += a;
			return res;
		}
		inline constexpr double_int &operator*=(const double_int &b) noexcept { *this = *this * b; return *this; }

		inline constexpr double_int &operator/=(const double_int &den) noexcept { *this = divmod(*this, den).first; return *this; }
		inline constexpr double_int &operator%=(const double_int &den) noexcept { *this = divmod(*this, den).second; return *this; }

		inline constexpr friend double_int operator/(const double_int &num, const double_int &den) noexcept { return divmod(num, den).first; }
		inline constexpr friend double_int operator%(const double_int &num, const double_int &den) noexcept { return divmod(num, den).second; }

		// equivalent to <=> but works for any version of C++
		inline constexpr friend int cmp(const double_int &a, const double_int &b) noexcept { return a.high < b.high ? -1 : a.high > b.high ? 1 : a.low < b.low ? -1 : a.low > b.low ? 1 : 0; }
		inline constexpr friend bool operator==(const double_int &a, const double_int &b) noexcept { return cmp(a, b) == 0; }
		inline constexpr friend bool operator!=(const double_int &a, const double_int &b) noexcept { return cmp(a, b) != 0; }
		inline constexpr friend bool operator<(const double_int &a, const double_int &b) noexcept { return cmp(a, b) < 0; }
		inline constexpr friend bool operator<=(const double_int &a, const double_int &b) noexcept { return cmp(a, b) <= 0; }
		inline constexpr friend bool operator>(const double_int &a, const double_int &b) noexcept { return cmp(a, b) > 0; }
		inline constexpr friend bool operator>=(const double_int &a, const double_int &b) noexcept { return cmp(a, b) >= 0; }

		inline friend std::ostream &operator<<(std::ostream &ostr, const double_int &val)
		{
			double_int copy = val;
			std::string str;
			int digit;
			
			// build the string
			if (ostr.flags() & std::ios::oct)
			{
				const double_int mask = 7;
				do
				{
					digit = copy & mask;
					copy >>= 3;
					str.push_back('0' + digit);
				}
				while (copy);
				if (ostr.flags() & std::ios::showbase) ostr.put('0');
			}
			else // default to hex mode
			{
				const double_int mask = 15;
				const char hex_alpha = ostr.flags() & std::ios::uppercase ? 'A' : 'a';
				do
				{
					digit = copy & mask;
					copy >>= 4;
					str.push_back(digit < 10 ? '0' + digit : hex_alpha + digit - 10);
				}
				while (copy);
				if (ostr.flags() & std::ios::showbase)
				{
					ostr.put('0'); // uses put() to make sure we don't clobber ostr.width()
					ostr.put('x');
				}
			}
			
			// write the string
			if (ostr.flags() & std::ios::left)
			{
				for (int i = str.size() - 1; i >= 0; --i) ostr.put(str[i]);
				for (int i = (int)ostr.width() - (int)str.size(); i > 0; --i) ostr.put(ostr.fill());
			}
			else // default to right
			{
				for (int i = (int)ostr.width() - (int)str.size(); i > 0; --i) ostr.put(ostr.fill());
				for (int i = str.size() - 1; i >= 0; --i) ostr.put(str[i]);
			}

			// clobber width - makes sure the next write isn't weird
			ostr.width(0);
			return ostr;
		}

	public: // -- utilities -- //

		inline constexpr friend std::pair<double_int, double_int> divmod(const double_int &num, const double_int &den) noexcept
		{
			if (!den) throw std::domain_error("divide by zero");

			auto res = std::make_pair<double_int, double_int>(0, 0);
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

		inline constexpr friend bool bit_test(const double_int &val, u64 bit) noexcept { bit &= bits - 1; return bit >= bits / 2 ? bit_test(val.high, bit - bits / 2) : bit_test(val.low, bit); }
		inline constexpr friend void bit_set(double_int &val, u64 bit) noexcept { bit &= bits - 1; if (bit >= bits / 2) bit_set(val.high, bit - bits / 2); else bit_set(val.low, bit); }

	};

	// returns a pair of <quotient, remainder>
	inline constexpr std::pair<std::uint8_t, std::uint8_t> divmod(std::uint8_t num, std::uint8_t den) { return {num / den, num % den}; }
	inline constexpr std::pair<std::uint16_t, std::uint16_t> divmod(std::uint16_t num, std::uint16_t den) { return {num / den, num % den}; }
	inline constexpr std::pair<std::uint32_t, std::uint32_t> divmod(std::uint32_t num, std::uint32_t den) { return {num / den, num % den}; }
	inline constexpr std::pair<std::uint64_t, std::uint64_t> divmod(std::uint64_t num, std::uint64_t den) { return {num / den, num % den}; }

	inline constexpr bool bit_test(u64 val, u64 bit) noexcept { return (val >> bit) & 1; }

	inline constexpr void bit_set(std::uint8_t &val, u64 bit) noexcept { val |= (std::uint8_t)1 << bit; }
	inline constexpr void bit_set(std::uint16_t &val, u64 bit) noexcept { val |= (std::uint16_t)1 << bit; }
	inline constexpr void bit_set(std::uint32_t &val, u64 bit) noexcept { val |= (std::uint32_t)1 << bit; }
	inline constexpr void bit_set(std::uint64_t &val, u64 bit) noexcept { val |= (std::uint64_t)1 << bit; }

}

#endif
