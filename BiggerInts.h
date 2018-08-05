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
	// -- helpers and types -- //

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

	// holds a type T where assignments to it are masked
	template<typename T, u64 bits> struct masked_single_int;

	// wraps larger integer sizes in terms of 2 smaller ones
	template<u64 bits, bool sign> struct double_int;

	// -- #hax typedef omega synergy highjinks - don't even need to be defined, boi -- //
	
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

	// gets an unsigned integer with the specified number of bits
	template<u64 bits> using uint_t = decltype(returns_proper_type<bits, false>());
	// gets a signed integer with the specified number of bits
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
	
	template<u64 bits, bool sign> struct double_int
	{
	public: // -- ctor/asgn -- //

		inline constexpr double_int() = default;

		inline constexpr double_int(const double_int&) = default;
		inline constexpr double_int &operator=(const double_int&) = default;

	public: // -- data -- //
		
		static_assert(is_pow2(bits) && bits >= 128, "double_int: bits was out of valid domain");

		typedef uint_t<bits / 2> half;

		half low, high;

	public: // -- conversion -- //

		template<typename U, std::enable_if_t<(bit_count<U>::value <= bits / 2), int> = 0>
		inline constexpr double_int(const U &val) noexcept : low(val), high(0) {}

		template<typename U, std::enable_if_t<(bit_count<U>::value <= bits / 2), int> = 0>
		inline constexpr double_int &operator=(const U &val) noexcept { low = val; high = 0; return *this; }

		template<typename U, std::enable_if_t<(bit_count<U>::value <= bits / 2), int> = 0>
		inline constexpr explicit operator U() const noexcept { return (U)low; }

		inline constexpr double_int(half _high, half _low) noexcept : high(_high), low(_low) {}

	public: // -- bool conversion -- //

		inline constexpr explicit operator bool() const noexcept { return low || high; }
		inline constexpr friend bool operator!(const double_int &a) noexcept { return !(bool)a; }
	};

	// -- operators -- //

	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator++(double_int<bits, sign> &a) noexcept { if (!++a.low) ++a.high; return a; }
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator--(double_int<bits, sign> &a) noexcept { if (!a.low) --a.high; --a.low; return a; }

	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> operator++(double_int<bits, sign> &a, int) noexcept { double_int<bits, sign> val = a; ++a; return val; }
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> operator--(double_int<bits, sign> &a, int) noexcept { double_int<bits, sign> val = a; --a; return val; }

	template<u64 bits, bool sign1, bool sign2> inline constexpr double_int<bits, sign1> &operator+=(double_int<bits, sign1> &a, const double_int<bits, sign2> &b) noexcept
	{
		a.low += b.low;
		if (a.low < b.low) ++a.high;
		a.high += b.high;
		return a;
	}
	template<u64 bits, bool sign1, bool sign2> inline constexpr double_int<bits, sign1> &operator-=(double_int<bits, sign1> &a, const double_int<bits, sign2> &b) noexcept
	{
		if (a.low < b.low) --a.high;
		a.low -= b.low;
		a.high -= b.high;
		return a;
	}

	template<u64 bits, bool sign1, bool sign2> inline constexpr double_int<bits, sign1> &operator&=(double_int<bits, sign1> &a, const double_int<bits, sign2> &b) noexcept { a.low &= b.low; a.high &= b.high; return a; }
	template<u64 bits, bool sign1, bool sign2> inline constexpr double_int<bits, sign1> &operator|=(double_int<bits, sign1> &a, const double_int<bits, sign2> &b) noexcept { a.low |= b.low; a.high |= b.high; return a; }
	template<u64 bits, bool sign1, bool sign2> inline constexpr double_int<bits, sign1> &operator^=(double_int<bits, sign1> &a, const double_int<bits, sign2> &b) noexcept { a.low ^= b.low; a.high ^= b.high; return a; }

	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator<<=(double_int<bits, sign> &val, u64 count) noexcept
	{
		// count masked to bits limit - no-op on zero
		if (count &= bits - 1)
		{
			if (count >= bits / 2)
			{
				val.high = val.low;
				val.high <<= count - bits / 2;
				val.low = 0;
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

	template<u64 bits> inline constexpr double_int<bits, false> &operator>>=(double_int<bits, false> &val, u64 count) noexcept
	{
		// count masked to bits limit - no-op on zero
		if (count &= bits - 1)
		{
			if (count >= bits / 2)
			{
				val.low = val.high;
				val.low >>= count - bits / 2;
				val.high = 0;
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
				val.high = bit_test(val.high, bits / 2 - 1) ? ~(decltype(val.high))0 : (decltype(val.high))0;
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

	template<u64 bits> inline constexpr double_int<bits, true> operator-(const double_int<bits, true> &a) noexcept { double_int<bits, true> res = ~a; ++res; return res; }

	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> operator~(const double_int<bits, sign> &a) noexcept { return {~a.high, ~a.low}; }

	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> operator+(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { double_int<bits, sign> res = a; res += b; return res; }
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> operator-(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { double_int<bits, sign> res = a; res -= b; return res; }

	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> operator&(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { double_int<bits, sign> res = a; res &= b; return res; }
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> operator|(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { double_int<bits, sign> res = a; res |= b; return res; }
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> operator^(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { double_int<bits, sign> res = a; res ^= b; return res; }

	template<u64 bits, bool sign, typename T> inline constexpr double_int<bits, sign> operator<<(const double_int<bits, sign> &a, const T &count) noexcept { double_int<bits, sign> res = a; res <<= count; return res; }
	template<u64 bits, bool sign, typename T> inline constexpr double_int<bits, sign> operator>>(const double_int<bits, sign> &a, const T &count) noexcept { double_int<bits, sign> res = a; res >>= count; return res; }

	// -- multiplication -- //

	template<u64 bits> inline constexpr double_int<bits, false> operator*(double_int<bits, false> a, const double_int<bits, false> &b) noexcept
	{
		double_int<bits, false> res = 0;
		for (u64 bit = 0; a; ++bit, a <<= 1)
			if (bit_test(b, bit)) res += a;
		return res;
	}
	template<u64 bits> inline constexpr double_int<bits, true> operator*(double_int<bits, true> a, const double_int<bits, true> &b) noexcept
	{
		// we'll do signed divmod in terms of unsigned

		double_int<bits, false> &ua = a;
		double_int<bits, false> &ub = b;

		bool neg = false;
		if (bit_test(unum, bits - 1)) { make_neg(ua); neg = !neg; }
		if (bit_test(uden, bits - 1)) { make_neg(ub); neg = !neg; }

		ua = ua * ub;

		if (neg) make_neg(ua);

		return ua;
	}

	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator*=(double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { a = a * b; return a; }

	// -- division -- //

	inline constexpr std::pair<std::uint8_t, std::uint8_t> divmod(std::uint8_t num, std::uint8_t den) { return {num / den, num % den}; }
	inline constexpr std::pair<std::uint16_t, std::uint16_t> divmod(std::uint16_t num, std::uint16_t den) { return {num / den, num % den}; }
	inline constexpr std::pair<std::uint32_t, std::uint32_t> divmod(std::uint32_t num, std::uint32_t den) { return {num / den, num % den}; }
	inline constexpr std::pair<std::uint64_t, std::uint64_t> divmod(std::uint64_t num, std::uint64_t den) { return {num / den, num % den}; }
	template<u64 bits> inline constexpr std::pair<double_int<bits, false>, double_int<bits, false>> divmod(const double_int<bits, false> &num, const double_int<bits, false> &den)
	{
		if (!den) throw std::domain_error("divide by zero");

		auto res = std::make_pair<double_int<bits, false>, double_int<bits, false>>(0, 0);
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
	
	inline constexpr std::pair<std::int8_t, std::int8_t> divmod(std::int8_t num, std::int8_t den) { return {num / den, num % den}; }
	inline constexpr std::pair<std::int16_t, std::int16_t> divmod(std::int16_t num, std::int16_t den) { return {num / den, num % den}; }
	inline constexpr std::pair<std::int32_t, std::int32_t> divmod(std::int32_t num, std::int32_t den) { return {num / den, num % den}; }
	inline constexpr std::pair<std::int64_t, std::int64_t> divmod(std::int64_t num, std::int64_t den) { return {num / den, num % den}; }
	template<u64 bits> inline constexpr std::pair<double_int<bits, true>, double_int<bits, true>> divmod(const double_int<bits, true> &num, const double_int<bits, true> &den)
	{
		// we'll do signed divmod in terms of unsigned

		double_int<bits, false> unum = num;
		double_int<bits, false> uden = num;

		bool neg = false;
		if (bit_test(unum, bits - 1)) { make_neg(unum); neg = !neg; }
		if (bit_test(uden, bits - 1)) { make_neg(uden); neg = !neg; }

		auto res = divmod(unum, uden);
		
		if (neg)
		{
			make_neg(res.first);
			make_neg(res.second);
		}

		return std::make_pair<double_int<bits, true>, double_int<bits, true>>(res.first, res.second);
	}

	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator/=(double_int<bits, sign> &num, const double_int<bits, sign> &den) { num = divmod(num, den).first; return num; }
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> &operator%=(double_int<bits, sign> &num, const double_int<bits, sign> &den) { num = divmod(num, den).second; return num; }

	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> operator/(const double_int<bits, sign> &num, const double_int<bits, sign> &den) { return divmod(num, den).first; }
	template<u64 bits, bool sign> inline constexpr double_int<bits, sign> operator%(const double_int<bits, sign> &num, const double_int<bits, sign> &den) { return divmod(num, den).second; }

	// -- comparison -- //

	// cmp() will be equivalent to <=> but works for any version of C++

	template<u64 bits> inline constexpr int cmp(const double_int<bits, false> &a, const double_int<bits, false> &b) noexcept { return a.high < b.high ? -1 : a.high > b.high ? 1 : a.low < b.low ? -1 : a.low > b.low ? 1 : 0; }
	template<u64 bits> inline constexpr int cmp(const double_int<bits, true> &a, const double_int<bits, true> &b) noexcept { int ucmp = a.high < b.high ? -1 : a.high > b.high ? 1 : a.low < b.low ? -1 : a.low > b.low ? 1 : 0; return bit_test(a, bits - 1) ^ bit_test(b, bits - 1) ? -ucmp : ucmp; }
	
	template<u64 bits, bool sign> inline constexpr bool operator==(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { return cmp(a, b) == 0; }
	template<u64 bits, bool sign> inline constexpr bool operator!=(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { return cmp(a, b) != 0; }
	template<u64 bits, bool sign> inline constexpr bool operator<(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { return cmp(a, b) < 0; }
	template<u64 bits, bool sign> inline constexpr bool operator<=(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { return cmp(a, b) <= 0; }
	template<u64 bits, bool sign> inline constexpr bool operator>(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { return cmp(a, b) > 0; }
	template<u64 bits, bool sign> inline constexpr bool operator>=(const double_int<bits, sign> &a, const double_int<bits, sign> &b) noexcept { return cmp(a, b) >= 0; }

	// -- io -- //

	template<u64 bits> inline std::ostream &operator<<(std::ostream &ostr, double_int<bits, false> val)
	{
		std::string str;
		int digit;
		
		// build the string
		if (ostr.flags() & std::ios::oct)
		{
			static constexpr double_int<bits, false> mask = 7;
			do
			{
				digit = (int)(val & mask);
				val >>= 3;
				str.push_back('0' + digit);
			}
			while (val);
			if (ostr.flags() & std::ios::showbase) ostr.put('0');
		}
		else if(ostr.flags() & std::ios::hex)
		{
			static constexpr double_int<bits, false> mask = 15;
			const char hex_alpha = ostr.flags() & std::ios::uppercase ? 'A' : 'a';
			do
			{
				digit = (int)(val & mask);
				val >>= 4;
				str.push_back(digit < 10 ? '0' + digit : hex_alpha + digit - 10);
			}
			while (val);
			if (ostr.flags() & std::ios::showbase)
			{
				ostr.put('0'); // uses put() to make sure we don't clobber ostr.width()
				ostr.put('x');
			}
		}
		else // default to dec mode
		{
			static constexpr double_int<bits, false> base = 10;
			std::pair<double_int<bits, false>, double_int<bits, false>> temp;
			do
			{
				temp = divmod(val, base);
				digit = (int)temp.second;
				val = temp.first;
				str.push_back('0' + digit);
			}
			while (val);
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
	template<u64 bits> inline std::ostream &operator<<(std::ostream &ostr, double_int<bits, true> val)
	{
		// we'll do signed io in terms of unsigned

		if (bit_test(val, bits - 1))
		{
			make_neg(val);
			if (ostr.flags() & std::ios::showpos) ostr.put('-');
		}
		else if (ostr.flags() & std::ios::showpos) ostr.put('+');

		// refer to val as unsigned
		ostr << *(double_int<bits, false>*)&val;

		return ostr;
	}

	// -- extra utilities -- //

	inline constexpr void make_not(std::uint64_t &val) noexcept { val = ~val; }
	template<u64 bits, bool sign> inline constexpr void make_not(double_int<bits, sign> &val) noexcept { make_not(val.low); make_not(val.high); }

	template<u64 bits, bool sign> inline constexpr void make_neg(double_int<bits, sign> &val) noexcept { make_not(val); ++val; }

	inline constexpr bool bit_test(std::uint8_t val, u64 bit) noexcept { return (val >> bit) & 1; }
	inline constexpr bool bit_test(std::uint16_t val, u64 bit) noexcept { return (val >> bit) & 1; }
	inline constexpr bool bit_test(std::uint32_t val, u64 bit) noexcept { return (val >> bit) & 1; }
	inline constexpr bool bit_test(std::uint64_t val, u64 bit) noexcept { return (val >> bit) & 1; }
	template<u64 bits, bool sign> inline constexpr bool bit_test(const double_int<bits, sign> &val, u64 bit) noexcept { bit &= bits - 1; return bit >= bits / 2 ? bit_test(val.high, bit - bits / 2) : bit_test(val.low, bit); }

	inline constexpr void bit_set(std::uint8_t &val, u64 bit) noexcept { val |= (std::uint8_t)1 << bit; }
	inline constexpr void bit_set(std::uint16_t &val, u64 bit) noexcept { val |= (std::uint16_t)1 << bit; }
	inline constexpr void bit_set(std::uint32_t &val, u64 bit) noexcept { val |= (std::uint32_t)1 << bit; }
	inline constexpr void bit_set(std::uint64_t &val, u64 bit) noexcept { val |= (std::uint64_t)1 << bit; }
	template<u64 bits, bool sign> inline constexpr void bit_set(double_int<bits, sign> &val, u64 bit) noexcept { bit &= bits - 1; if (bit >= bits / 2) bit_set(val.high, bit - bits / 2); else bit_set(val.low, bit); }
}

#endif
