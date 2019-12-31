#include <cctype>
#include <sstream>
#include <algorithm>

#include "BiggerInts.h"

using namespace BiggerInts;

void detail::_collapse(bigint &a)
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

bool detail::is_neg(const bigint &val) noexcept { return !val.blocks.empty() && (val.blocks.back() & 0x8000000000000000); }

std::size_t detail::highest_set_bit(const bigint &val) noexcept
{
	if (val.blocks.empty()) return 0;
	if (val.blocks.back()) return (val.blocks.size() - 1) * 64 + highest_set_bit(val.blocks.back());
	else return (val.blocks.size() - 2) * 64 + highest_set_bit(val.blocks[val.blocks.size() - 2]);
}
void detail::make_not(bigint &a)
{
	if (a.blocks.empty()) a.blocks.push_back(0xffffffffffffffffull);
	else
	{
		for (std::uint64_t &v : a.blocks) v = ~v;
		_collapse(a);
	}
}
void detail::make_neg(bigint &a) { make_not(a); ++a; }

bool detail::bit_test(const bigint &val, std::uint64_t bit) noexcept
{
	std::uint64_t block = bit / 64;
	if (block >= val.blocks.size()) return detail::is_neg(val); // lookup in ghost blocks
	return (val.blocks[block] >> (bit % 64)) & 1; // lookup in declared range
}
bool detail::bit_test_in_bounds_nonzero(const bigint &val, std::uint64_t bit) noexcept
{
	return (val.blocks[bit / 64] >> (bit % 64)) & 1;
}
void detail::bit_set(bigint &val, std::uint64_t bit) noexcept
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

bool detail::cmp_less_non_negative(const bigint &a, const bigint &b) noexcept // same as cmp() but requires they both be non-negative
{
	if (a.blocks.size() > b.blocks.size()) return false;
	if (a.blocks.size() < b.blocks.size()) return true;
	for (std::size_t i = a.blocks.size(); i-- > 0; ) if (a.blocks[i] != b.blocks[i]) return a.blocks[i] < b.blocks[i];
	return false;
}
int detail::cmp(const bigint &a, const bigint &b) noexcept
{
	const bool a_neg = detail::is_neg(a);
	const bool b_neg = detail::is_neg(b);
	if (a_neg ^ b_neg) return a_neg ? -1 : 1; // if they have different signs we know right away

	if (a.blocks.size() > b.blocks.size()) return a_neg ? -1 : 1; // if a has higher magnitude than b we can tell from the sign
	if (a.blocks.size() < b.blocks.size()) return a_neg ? 1 : -1; // similarly

	for (std::size_t i = a.blocks.size(); i-- > 0; ) if (a.blocks[i] != b.blocks[i]) return a.blocks[i] < b.blocks[i] ? -1 : 1; // otherwise base decision on first different block

	return 0; // otherwise they're equal
}
int detail::cmp_bigint_builtin(const bigint &a, long long val) noexcept
{
	const bool a_neg = detail::is_neg(a);
	if (a_neg != (val < 0)) return a_neg ? -1 : 1; // if they have different signs we know right away

	if (a.blocks.size() > 1) return a_neg ? -1 : 1; // if a has higher magnitude than b we can tell from the sign
	else if (a.blocks.size() == 1)
	{
		if (a.blocks[0] == (unsigned long long)val) return 0;
		else return a.blocks[0] < (unsigned long long)val ? -1 : 1;
	}
	else return val < 0 ? 1 : val > 0 ? -1 : 0; // otherwise a == 0
}
int detail::cmp_bigint_builtin(const bigint &a, unsigned long long val) noexcept
{
	if (detail::is_neg(a)) return -1; // if a is negative we know it's < immediately

	if (a.blocks.size() > 2) return 1; // if a has higher magnitude, it's greater
	else if (a.blocks.size() == 2 && a.blocks[1]) return 1; // if second block is significant, a is greater
	else if (a.blocks.size() == 0) return val < 0 ? 1 : val > 0 ? -1 : 0; // if a == 0
	else return a.blocks[0] < val ? -1 : a.blocks[0] > val ? 1 : 0; // otherwise same magnitude
}

void bigint::zero_extend(detail::const_fixed_int_wrapper val)
{
	blocks.assign(val.blocks, val.blocks + val.blocks_n);
	std::size_t s = blocks.size();
	for (; s > 0; --s) if (blocks[s - 1]) break;
	blocks.resize(s);
	if (!blocks.empty() && (blocks.back() & 0x8000000000000000ull)) blocks.push_back(0ull);
}
void bigint::sign_extend(detail::const_fixed_int_wrapper val)
{
	blocks.assign(val.blocks, val.blocks + val.blocks_n);
	_collapse(*this);
}

bigint &bigint::operator++()
{
	for (std::size_t i = 1; i < blocks.size(); ++i) if (++blocks[i - 1]) { _collapse(*this); return *this; }

	if (blocks.empty()) blocks.push_back(1ull);
	else
	{
		std::uint64_t high = ++blocks.back();
		if (high == 0) blocks.clear();
		else if (high == 0x8000000000000000ull) blocks.push_back(0ull);
		else _collapse(*this);
	}

	return *this;
}
bigint bigint::operator++(int) { bigint cpy = *this; ++*this; return cpy; }

bigint &bigint::operator--()
{
	for (std::size_t i = 1; i < blocks.size(); ++i) if (blocks[i - 1]--) { _collapse(*this); return *this; }

	if (blocks.empty()) blocks.push_back(0xffffffffffffffffull);
	else
	{
		std::uint64_t high = blocks.back()--;
		/* high == 0 case cannot happen because 0 is represented by empty array and is therefore handled above */
		if (high == 0x8000000000000000ull) blocks.push_back(0xffffffffffffffffull);
		else _collapse(*this);
	}

	return *this;
}
bigint bigint::operator--(int) { bigint cpy = *this; --*this; return cpy; }

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

bigint &bigint::operator+=(const bigint &b) { return _add<false>(*this, b); }
bigint &bigint::operator-=(const bigint &b) { return _add<true>(*this, b); }

bigint detail::operator+(const bigint &a, const bigint &b) { bigint cpy = a; cpy += b; return cpy; }
bigint detail::operator+(bigint &&a, const bigint &b) { bigint cpy = std::move(a); cpy += b; return cpy; }
bigint detail::operator+(const bigint &a, bigint &&b) { bigint cpy = std::move(b); cpy += a; return cpy; }
bigint detail::operator+(bigint &&a, bigint &&b)
{
	if (a.blocks.size() >= b.blocks.size()) { bigint cpy = std::move(a); cpy += b; return cpy; }
	else { bigint cpy = std::move(b); cpy += a; return cpy; }
}

bigint detail::operator-(const bigint &a, const bigint &b) { bigint cpy = a; cpy -= b; return cpy; }
bigint detail::operator-(bigint &&a, const bigint &b) { bigint cpy = std::move(a); cpy -= b; return cpy; }
bigint detail::operator-(const bigint &a, bigint &&b) { bigint cpy = std::move(b); cpy -= a; detail::make_neg(cpy); return cpy; } // not sure if this is faster than just making a dynamically-allocated copy
bigint detail::operator-(bigint &&a, bigint &&b) { bigint cpy = std::move(a); cpy -= b; return cpy; }

bigint &bigint::operator&=(const bigint &b)
{
	const bool a_neg = detail::is_neg(*this);
	const bool b_neg = detail::is_neg(b);

	const std::size_t min = std::min(blocks.size(), b.blocks.size());

	for (std::size_t i = 0; i < min; ++i) blocks[i] &= b.blocks[i]; // perform on range declared by both

	if (blocks.size() < b.blocks.size())
	{
		// if a was positive then rest is all 0's and thus result is already computed
		if (a_neg) // otherwise a was negative, so rest is all 1's - just copy remaining blocks from b
		{
			blocks.resize(b.blocks.size());
			for (std::size_t i = min; i < b.blocks.size(); ++i) blocks[i] = b.blocks[i];
		}
	}
	else if (b.blocks.size() < blocks.size())
	{
		// if b was negative then rest is all 1's - we keep everything in a that was already there
		if (!b_neg) // otherwise b was positive - rest is all 0's and thus result is already computed - just truncate down to computed region and make sure result is positive
		{
			blocks.resize(min);
			if (detail::is_neg(*this)) blocks.push_back(0ull);
		}
	}

	_collapse(*this);
	return *this;
}

bigint detail::operator&(const bigint &a, const bigint &b) { bigint cpy = a; cpy &= b; return cpy; }
bigint detail::operator&(bigint &&a, const bigint &b) { bigint cpy = std::move(a); cpy &= b; return cpy; }
bigint detail::operator&(const bigint &a, bigint &&b) { bigint cpy = std::move(b); cpy &= a; return cpy; }
bigint detail::operator&(bigint &&a, bigint &&b) { bigint cpy = std::move(a); cpy &= b; return cpy; }

bigint &bigint::operator|=(const bigint &b)
{
	const bool a_neg = detail::is_neg(*this);
	const bool b_neg = detail::is_neg(b);

	const std::size_t min = std::min(blocks.size(), b.blocks.size());

	for (std::size_t i = 0; i < min; ++i) blocks[i] |= b.blocks[i]; // perform on range declared by both

	if (blocks.size() < b.blocks.size())
	{
		// if a was negative then rest is all 1's, and thus result is already computed
		if (!a_neg) // otherwise a positive, so rest is all 0's - just copy remaining blocks from b
		{
			blocks.resize(b.blocks.size());
			for (std::size_t i = min; i < b.blocks.size(); ++i) blocks[i] = b.blocks[i];
		}
	}
	else if (b.blocks.size() < blocks.size())
	{
		// if b was positive then rest is all 0's - we keep everything in a that was already there
		if (b_neg) // otherwise b was negative, so rest is all 1's - just truncate down to computed region and make sure result is negative
		{
			blocks.resize(min);
			if (!detail::is_neg(*this)) blocks.push_back(0xffffffffffffffffull);
		}
	}

	_collapse(*this);
	return *this;
}

bigint detail::operator|(const bigint &a, const bigint &b) { bigint cpy = a; cpy |= b; return cpy; }
bigint detail::operator|(bigint &&a, const bigint &b) { bigint cpy = std::move(a); cpy |= b; return cpy; }
bigint detail::operator|(const bigint &a, bigint &&b) { bigint cpy = std::move(b); cpy |= a; return cpy; }
bigint detail::operator|(bigint &&a, bigint &&b) { bigint cpy = std::move(a); cpy |= b; return cpy; }

bigint &bigint::operator^=(const bigint &b)
{
	const bool a_neg = detail::is_neg(*this);
	const bool b_neg = detail::is_neg(b);

	const std::size_t min = std::min(blocks.size(), b.blocks.size());

	for (std::size_t i = 0; i < min; ++i) blocks[i] ^= b.blocks[i]; // perform on range declared by both

	if (blocks.size() < b.blocks.size())
	{
		if (a_neg) // if a was negative then rest is all 1's - just copy the inverted bits of b
		{
			blocks.resize(b.blocks.size());
			for (std::size_t i = min; i < b.blocks.size(); ++i) blocks[i] = ~b.blocks[i];
		}
		else // otherwise a was positive, so rest is all 0's - just copy the bits of b
		{
			blocks.resize(b.blocks.size());
			for (std::size_t i = min; i < b.blocks.size(); ++i) blocks[i] = b.blocks[i];
		}
	}
	else if (b.blocks.size() < blocks.size())
	{
		// if b was positive then rest is all 0's - so we already have the result
		if (b_neg) // otherwise b negative, so rest is all 1's - just invert the bits of a
		{
			for (std::size_t i = min; i < blocks.size(); ++i) blocks[i] = ~blocks[i];
		}
	}

	_collapse(*this);
	return *this;
}

bigint detail::operator^(const bigint &a, const bigint &b) { bigint cpy = a; cpy ^= b; return cpy; }
bigint detail::operator^(bigint &&a, const bigint &b) { bigint cpy = std::move(a); cpy ^= b; return cpy; }
bigint detail::operator^(const bigint &a, bigint &&b) { bigint cpy = std::move(b); cpy ^= a; return cpy; }
bigint detail::operator^(bigint &&a, bigint &&b) { bigint cpy = std::move(a); cpy ^= b; return cpy; }

bigint &bigint::operator<<=(std::uint64_t count)
{
	if (!blocks.empty())
	{
		const std::uint64_t fill = detail::is_neg(*this) ? -1 : 0;
		const std::size_t full = count / 64;
		if (full)
		{
			blocks.resize(blocks.size() + full); // make room for the shifting blocks
			for (std::size_t i = blocks.size(); i-- > full; ) blocks[i] = blocks[i - full];
			for (std::size_t i = full; i-- > 0; ) blocks[i] = 0;
			count &= 63;
		}
		if (count)
		{
			std::uint64_t overflow_block = (fill << count) | (blocks.back() >> (64 - count)); // compute overflow block from high block (include ghost bits based on sign)
			for (std::size_t i = blocks.size() - 1; i > full; --i) blocks[i] = (blocks[i] << count) | (blocks[i - 1] >> (64 - count));
			blocks[full] <<= count;

			const std::uint64_t new_fill = detail::is_neg(*this) ? -1 : 0; // get new fill pattern
			if (overflow_block != new_fill) blocks.push_back(overflow_block); // if the overflow block isn't generated by the new fill pattern, we need to append it
		}
	}
	return *this;
}

bigint detail::operator<<(const bigint &val, std::uint64_t count) { bigint cpy = val; cpy <<= count; return cpy; }
bigint detail::operator<<(bigint &&val, std::uint64_t count) { bigint cpy = std::move(val); cpy <<= count; return cpy; }

bigint &bigint::operator>>=(std::uint64_t count)
{
	if (!blocks.empty())
	{
		const std::uint64_t fill = detail::is_neg(*this) ? -1 : 0;
		const std::size_t full = count / 64;
		if (full)
		{
			for (std::size_t i = 0; i < blocks.size() - full; ++i) blocks[i] = blocks[i + full];
			for (std::size_t i = blocks.size() - full; i < blocks.size(); ++i) blocks[i] = fill;
			count &= 63;
		}
		if (count)
		{
			for (std::size_t i = full; i < blocks.size() - 1; ++i) blocks[i] = (blocks[i] >> count) | (blocks[i + 1] << (64 - count));
			blocks[blocks.size() - 1] = (blocks[blocks.size() - 1] >> count) | (fill << (64 - count));
		}
		_collapse(*this);
	}
	return *this;
}

bigint detail::operator>>(const bigint &val, std::uint64_t count) { bigint cpy = val; cpy >>= count; return cpy; }
bigint detail::operator>>(bigint &&val, std::uint64_t count) { bigint cpy = std::move(val); cpy >>= count; return cpy; }

bigint detail::operator+(const bigint &a) { return a; }
bigint detail::operator+(bigint &&a) { return std::move(a); }

bigint detail::operator-(const bigint &a) { bigint res = a; detail::make_neg(res); return res; }
bigint detail::operator-(bigint &&a) { bigint res = std::move(a); detail::make_neg(res); return res; }

bigint detail::operator~(const bigint &a) { bigint res = a; detail::make_not(res); return res; }
bigint detail::operator~(bigint &&a) { bigint res = std::move(a); detail::make_not(res); return res; }

bigint _multiply_positive(const bigint &a, const bigint &b)
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

		return _multiply_positive(a_cpy, b_cpy);
	}
	else if (a_neg)
	{
		bigint a_cpy = std::forward<U>(a);
		detail::make_neg(a_cpy);

		bigint res = _multiply_positive(a_cpy, b);
		detail::make_neg(res);
		return res;
	}
	else if (b_neg)
	{
		bigint b_cpy = std::forward<V>(b);
		detail::make_neg(b_cpy);

		bigint res = _multiply_positive(a, b_cpy);
		detail::make_neg(res);
		return res;
	}
	else return _multiply_positive(a, b);
}

bigint detail::operator*(const bigint &a, const bigint &b) { return _multiply_unknown(a, b); }
bigint detail::operator*(bigint &&a, const bigint &b) { return _multiply_unknown(std::move(a), b); }
bigint detail::operator*(const bigint &a, bigint &&b) { return _multiply_unknown(a, std::move(b)); }
bigint detail::operator*(bigint &&a, bigint &&b) { return _multiply_unknown(std::move(a), std::move(b)); }

bigint &bigint::operator*=(const bigint &b)
{
	if (this != &b) *this = std::move(*this) * b;
	else *this = *this * b;
	return *this;
}
bigint &bigint::operator*=(bigint &&b) { *this = std::move(*this) * std::move(b); return *this; }

void detail::divmod_unchecked_positive_already_zero(bigint &quo, bigint &rem, const bigint &num, const bigint &den)
{
	// both positive, so if num < den then we know result immediately (this also handles num == 0 case, as we already know den != 0)
	if (num.blocks.size() < den.blocks.size()) { rem = num; return; }

	// pre-allocate space for the results
	quo.blocks.reserve(num.blocks.size() - den.blocks.size() + 1);
	rem.blocks.reserve(den.blocks.size());

	std::uint64_t bit = highest_set_bit(num);

	const std::size_t den_highest_bit = highest_set_bit(den);
	std::size_t res_second_highest_bit = 0;
	std::size_t shift_count = 0;

	while (true)
	{
		++shift_count;
		if (res_second_highest_bit + shift_count >= den_highest_bit)
		{
			rem <<= shift_count;
			res_second_highest_bit += shift_count;
			shift_count = 0;

			if (bit_test_in_bounds_nonzero(num, bit)) ++rem;

			if (!cmp_less_non_negative(rem, den))
			{
				rem -= den;
				res_second_highest_bit = highest_set_bit(rem);
				bit_set(quo, bit); // no unchecked equivalent because arbitrary precision
			}
		}
		else if (bit_test_in_bounds_nonzero(num, bit))
		{
			rem <<= shift_count;
			res_second_highest_bit += shift_count;
			shift_count = 0;
			++rem;
		}

		if (bit-- == 0) break;
	}

	rem <<= shift_count;
}

std::pair<bigint, bigint> detail::divmod_unchecked_positive(const bigint &num, const bigint &den)
{
	std::pair<bigint, bigint> res;
	detail::divmod_unchecked_positive_already_zero(res.first, res.second, num, den);
	return res;
}

template<typename U, typename V, std::enable_if_t<std::is_same_v<std::decay_t<U>, bigint> && std::is_same_v<std::decay_t<V>, bigint>, int> = 0>
std::pair<bigint, bigint> _divmod_unchecked(U &&a, V &&b)
{
	const bool a_neg = detail::is_neg(a);
	const bool b_neg = detail::is_neg(b);

	if (a_neg && b_neg)
	{
		bigint a_cpy = std::forward<U>(a);
		bigint b_cpy = std::forward<V>(b);
		detail::make_neg(a_cpy);
		detail::make_neg(b_cpy);

		auto res = detail::divmod_unchecked_positive(a_cpy, b_cpy);
		detail::make_neg(res.second);
		return res;
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
		return res;
	}
	else return detail::divmod_unchecked_positive(a, b);
}

template<typename U, typename V, std::enable_if_t<std::is_same_v<std::decay_t<U>, bigint> && std::is_same_v<std::decay_t<V>, bigint>, int> = 0>
std::pair<bigint, bigint> _divmod(U &&a, V &&b)
{
	if (!b) throw std::domain_error("divide by zero");
	return _divmod_unchecked(std::forward<U>(a), std::forward<V>(b));
}

std::pair<bigint, bigint> detail::divmod(const bigint &a, const bigint &b) { return _divmod(a, b); }
std::pair<bigint, bigint> detail::divmod(bigint &&a, const bigint &b) { return _divmod(std::move(a), b); }
std::pair<bigint, bigint> detail::divmod(const bigint &a, bigint &&b) { return _divmod(a, std::move(b)); }
std::pair<bigint, bigint> detail::divmod(bigint &&a, bigint &&b) { return _divmod(std::move(a), std::move(b)); }

bigint detail::operator/(const bigint &num, const bigint &den) { return detail::divmod(num, den).first; }
bigint detail::operator/(bigint &&num, const bigint &den) { return detail::divmod(std::move(num), den).first; }
bigint detail::operator/(const bigint &num, bigint &&den) { return detail::divmod(num, std::move(den)).first; }
bigint detail::operator/(bigint &&num, bigint &&den) { return detail::divmod(std::move(num), std::move(den)).first; }

bigint detail::operator%(const bigint &num, const bigint &den) { return detail::divmod(num, den).second; }
bigint detail::operator%(bigint &&num, const bigint &den) { return detail::divmod(std::move(num), den).second; }
bigint detail::operator%(const bigint &num, bigint &&den) { return detail::divmod(num, std::move(den)).second; }
bigint detail::operator%(bigint &&num, bigint &&den) { return detail::divmod(std::move(num), std::move(den)).second; }

bigint &bigint::operator/=(const bigint &den)
{
	if (this != &den) *this = std::move(*this) / den;
	else *this = *this / den;
	return *this;
}
bigint &bigint::operator/=(bigint &&den) { *this = std::move(*this) / std::move(den); return *this; }

bigint &bigint::operator%=(const bigint &den)
{
	if (this != &den) *this = std::move(*this) % den;
	else *this = *this % den;
	return *this;
}
bigint &bigint::operator%=(bigint &&den) { *this = std::move(*this) % std::move(den); return *this; }

bigint bigint::pow(bigint a, const bigint &b) // pass by value is intentional
{
	if (detail::is_neg(b)) return {}; // if exponent is negative just return 0

	bigint res = 1;
	std::uint64_t high_bit = highest_set_bit(b);

	for (std::uint64_t bit = 0; bit <= high_bit + 1; ++bit)
	{
		if (detail::bit_test(b, bit)) res *= a;
		a *= a;
	}

	return res;
}

static constexpr std::uint64_t _factorial_lookups[] = {
	1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200,
	1307674368000, 20922789888000, 355687428096000, 6402373705728000, 121645100408832000, 2432902008176640000
};
static constexpr std::uint64_t _factorial_lookups_count = sizeof(_factorial_lookups) / sizeof(*_factorial_lookups);
bigint bigint::factorial(std::uint64_t v)
{
	if (v < _factorial_lookups_count) return _factorial_lookups[v];
	bigint res = _factorial_lookups[_factorial_lookups_count - 1];
	for (std::uint64_t i = _factorial_lookups_count; i <= v; ++i) res *= i;
	return res;
}

// given a hex character, converts it to an integer [0, 15] - returns true if it was a valid hex digit. ch is only meaningful on success.
constexpr bool ext_hex(int &ch) noexcept
{
	if (ch >= '0' && ch <= '9') { ch -= '0'; return true; }
	ch |= 32;
	if (ch >= 'a' && ch <= 'f') { ch = ch - 'a' + 10; return true; }
	return false;
}

// prints the value - interpreted as non-negative - sign_ch (if not null) is appended to the front of the printed value
template<typename T, std::enable_if_t<std::is_same_v<T, detail::fixed_int_wrapper> || std::is_same_v<T, bigint>, int> = 0>
std::ostream &print_positive(std::ostream &ostr, T val, char sign_ch)
{
	std::ostream::sentry sentry(ostr);
	if (!sentry) return ostr;

	std::string str;
	int digit, dcount;
	std::uint64_t block;

	// build the string
	if (ostr.flags() & std::ios::hex)
	{
		const bool sign = detail::is_neg(val);
		const char hex_alpha = ostr.flags() & std::ios::uppercase ? 'A' : 'a';

		(void)sign; // suppress unused warnings
		if constexpr (std::is_same_v<T, bigint>)
		{
			if (sign) val.blocks.push_back(0ull); // fix the issue of using arithmetic right shifts for bigint (negative would be infinite loop)
		}

		while (true)
		{
			// get a block
			if constexpr (std::is_same_v<T, bigint>)
			{
				block = val.blocks.empty() ? 0ull : val.blocks[0];
				val >>= 64;
			}
			else
			{
				block = val.blocks[0];
				shr(val, 64);
			}
			dcount = 0;

			do // write the block - do-while to ensure 0 gets printed
			{
				digit = block & 15;
				str.push_back(digit < 10 ? '0' + digit : hex_alpha + digit - 10);
				block >>= 4;
				++dcount;
			} while (block);

			if (detail::nonzero(val)) { for (; dcount < 16; ++dcount) str.push_back('0'); } // if there's still stuff, pad with zeroes and continue
			else break; // otherwise we're done
		}
		if constexpr (std::is_same_v<T, bigint>)
		{
			if (sign)
			{
				const char f = ostr.flags() & std::ios::uppercase ? 'F' : 'f';

				std::size_t s; // now we truncate all leading f's
				for (s = str.size(); s > 0; --s) if (str[s - 1] != f) break;
				str.resize(s); // chop off everything we don't want
				if (str.empty() || str.back() <= '7') str.push_back(f); // if the string is now empty or no longer sign extends to negative, add a (single) f back on
			}
			else if (str.back() >= '8') str.push_back('0'); // if val was positive but ends in an 8 or higher, add a 0 to prevent sign extending to negative on parsing
		}
		if (ostr.flags() & std::ios::showbase)
		{
			str.push_back('x'); // uses put() to make sure we don't clobber ostr.width()
			str.push_back('0');
		}
	}
	else if (ostr.flags() & std::ios::oct)
	{
		const bool sign = detail::is_neg(val);

		(void)sign; // suppress unused warnings
		if constexpr (std::is_same_v<T, bigint>)
		{
			if (sign) val.blocks.push_back(0ull); // fix the issue of using arithmetic right shifts for bigint (negative would be infinite loop)
		}

		while (true)
		{
			// get a block
			if constexpr (std::is_same_v<T, bigint>)
			{
				block = val.blocks.empty() ? 0ull : val.blocks[0];
				val >>= 63;
			}
			else
			{
				block = val.blocks[0];
				shr(val, 63);
			}
			block &= 0777777777777777777777ull;
			dcount = 0;

			// write the block - do-while to ensure 0 gets printed
			do
			{
				str.push_back('0' + (block & 7));
				block >>= 3;
				++dcount;
			} while (block);

			if (detail::nonzero(val)) { for (; dcount < 21; ++dcount) str.push_back('0'); } // if there's still stuff, pad with zeroes and continue
			else break; // otherwise we're done
		}
		if constexpr (std::is_same_v<T, bigint>)
		{
			if (sign)
			{
				// if it was negative we need to ensure that the high char will sign extend to negative
				if (str.back() == '1') str.back() = '7'; // 001 -> 111
				else if (str.back() <= '3') str.back() += 4; // 010, 011 -> 110, 111
				// otherwise is at least 4 and therefore already good to go

				std::size_t s; // now we truncate all leading 7's
				for (s = str.size(); s > 0; --s) if (str[s - 1] != '7') break;
				str.resize(s); // chop off everything we don't want
				if (str.empty() || str.back() <= '3') str.push_back('7'); // if the string is now empty or no longer sign extends to negative, add a (single) 7 back on
			}
			else if (str.back() >= '4') str.push_back('0'); // if val was positive but ends in a 4 or higher, add a 0 to prevent sign extending to negative on parsing
		}
		if (ostr.flags() & std::ios::showbase) str.push_back('0');
	}
	else // default to dec mode
	{
		std::conditional_t<std::is_same_v<T, detail::fixed_int_wrapper>, std::vector<std::uint64_t>, int> buffer; // buffer area for use with fixed_int
		std::pair<T, T> temp;
		T base;
		if constexpr (std::is_same_v<T, detail::fixed_int_wrapper>)
		{
			buffer.resize(val.blocks_n * 3, 0); // allocate space for temp.first, temp.second, and base
			temp.first = { &buffer[0], val.blocks_n };
			temp.second = { &buffer[val.blocks_n], val.blocks_n };
			base = { &buffer[2 * val.blocks_n], val.blocks_n };
			base.blocks[0] = 10000000000000000000ull; // initialize base
		}
		else
		{
			(void)buffer; // suppress unused warnings (will just be optimized away)
			base = 10000000000000000000ull; // initialize base
		}

		while (true)
		{
			// get a block
			if constexpr (std::is_same_v<T, detail::fixed_int_wrapper>)
			{
				for (std::size_t i = 0; i < val.blocks_n * 2; ++i) buffer[i] = 0; // zero temp for the next step
				detail::divmod_unchecked_positive_already_zero(temp.first, temp.second, val, base);
				block = temp.second.blocks[0];
				for (std::size_t i = 0; i < val.blocks_n; ++i) val.blocks[i] = buffer[i]; // copy quotient back to val
			}
			else
			{
				using std::swap; // import std::swap for adl

				temp.first.blocks.clear(); // zero temp for the next step
				temp.second.blocks.clear();
				detail::divmod_unchecked_positive_already_zero(temp.first, temp.second, val, base);
				block = temp.second.blocks.empty() ? 0ull : temp.second.blocks[0];
				swap(val, temp.first); // copy quotient back to val (swap is so we don't reallocate a vector each time we loop)
			}
			dcount = 0;

			// write the block - do-while to ensure 0 gets printed
			do
			{
				str.push_back('0' + block % 10);
				block /= 10;
				++dcount;
			} while (block);

			if (detail::nonzero(val)) { for (; dcount < 19; ++dcount) str.push_back('0'); } // if there's still stuff, pad with zeroes and continue
			else break; // otherwise we're done
		}
	}

	if (sign_ch) str.push_back(sign_ch); // append the sign character if specified

	std::reverse(str.begin(), str.end()); // reverse the string for printing
	ostr << str;                          // print the string (this applies the width/fill)

	return ostr;
}
template<typename T, std::enable_if_t<std::is_same_v<T, detail::fixed_int_wrapper> || std::is_same_v<T, bigint>, int> = 0>
std::ostream &print_signed(std::ostream &ostr, T val)
{
	char sign_ch = 0;

	if ((ostr.flags() & std::ios::basefield) == std::ios::dec) // decimal is the only base where signs come into play
	{
		if (detail::is_neg(val))
		{
			detail::make_neg(val);
			sign_ch = '-';
		}
		else if (ostr.flags() & std::ios::showpos) sign_ch = '+';
	}

	print_positive(ostr, val, sign_ch);
	return ostr;
}

void detail::print_positive_core(std::ostream &ostr, const_fixed_int_wrapper val)
{
	std::vector<std::uint64_t> buffer(val.blocks_n);
	for (std::size_t i = 0; i < val.blocks_n; ++i) buffer[i] = val.blocks[i]; // copy to a location we can modify
	print_positive<fixed_int_wrapper>(ostr, { &buffer[0], buffer.size() }, 0);
}
void detail::print_signed_core(std::ostream &ostr, const_fixed_int_wrapper val)
{
	std::vector<std::uint64_t> buffer(val.blocks_n);
	for (std::size_t i = 0; i < val.blocks_n; ++i) buffer[i] = val.blocks[i]; // copy to a location we can modify
	print_signed<fixed_int_wrapper>(ostr, { &buffer[0], buffer.size() });
}

std::ostream &detail::operator<<(std::ostream &ostr, const bigint &val) { return print_signed<bigint>(ostr, val); }
std::ostream &detail::operator<<(std::ostream &ostr, bigint &&val) { return print_signed<bigint>(ostr, std::move(val)); }

template<typename T, std::enable_if_t<std::is_same_v<T, detail::fixed_int_wrapper> || std::is_same_v<T, bigint>, int> = 0>
std::istream &parse_positive_hex(std::istream &istr, T &val, bool noskipws = false)
{
	std::istream::sentry sentry(istr, noskipws);
	if (!sentry) return istr;

	// start by zeroing value
	if constexpr (std::is_same_v<T, bigint>) val.blocks.clear();
	else for (std::size_t i = 0; i < val.blocks_n; ++i) val.blocks[i] = 0;

	int first_digit, digit, num_digits;
	std::uint64_t block;

	// first char must be a hex digit - don't extract
	if ((first_digit = istr.peek()) == EOF) { istr.setstate(std::ios::failbit | std::ios::eofbit); return istr; }
	if (!ext_hex(first_digit)) { istr.setstate(std::ios::failbit); return istr; }

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
			if constexpr (std::is_same_v<T, detail::fixed_int_wrapper>)
			{
				// detect overflow
				if (val.blocks[val.blocks_n - 1] >> (std::uint64_t)(64 - 4 * num_digits)) { istr.setstate(std::ios::failbit); return istr; }

				// incorporate it into value
				shl(val, (std::uint64_t)(4 * num_digits));
				val.blocks[0] |= block;
			}
			else
			{
				// incorporate it into value
				val <<= (std::uint64_t)(4 * num_digits);
				val |= block;
			}
		}

		if (num_digits < 16) break;
	}

	if constexpr (std::is_same_v<T, bigint>)
	{
		static_assert(detail::highest_set_bit(0) == 0); // sanity check for the behavior below

		if (first_digit >= 8) // if the value should sign extend to negative
		{
			val.blocks.back() |= ~((1ull << detail::highest_set_bit(val.blocks.back())) - 1); // set all bits after the high bit to one
			detail::_collapse(val);
		}
		else // otherwise should sign extend to positive
		{
			if (detail::is_neg(val)) val.blocks.push_back(0ull); // if currently negative, toss on a 0 block to make it positive
		}
	}

	return istr;
}
template<typename T, std::enable_if_t<std::is_same_v<T, detail::fixed_int_wrapper> || std::is_same_v<T, bigint>, int> = 0>
std::istream &parse_positive_oct(std::istream &istr, T &val, bool noskipws = false)
{
	std::istream::sentry sentry(istr, noskipws);
	if (!sentry) return istr;

	// start by zeroing value
	if constexpr (std::is_same_v<T, bigint>) val.blocks.clear();
	else for (std::size_t i = 0; i < val.blocks_n; ++i) val.blocks[i] = 0;

	int first_digit, digit, num_digits;
	std::uint64_t block;

	// first char must be an oct digit - don't extract
	if ((first_digit = istr.peek()) == EOF) { istr.setstate(std::ios::failbit | std::ios::eofbit); return istr; }
	if (first_digit < '0' || first_digit > '7') { istr.setstate(std::ios::failbit); return istr; }
	first_digit -= '0';

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
			if constexpr (std::is_same_v<T, detail::fixed_int_wrapper>)
			{
				// detect overflow
				if (val.blocks[val.blocks_n - 1] >> (std::uint64_t)(64 - 3 * num_digits)) { istr.setstate(std::ios::failbit); return istr; }

				// incorporate it into value
				shl(val, (std::uint64_t)(3 * num_digits));
				val.blocks[0] |= block;
			}
			else
			{
				// incorporate it into value
				val <<= (std::uint64_t)(3 * num_digits);
				val |= block;
			}
		}

		if (num_digits < 21) break;
	}

	if constexpr (std::is_same_v<T, bigint>)
	{
		static_assert(detail::highest_set_bit(0) == 0); // sanity check for the behavior below

		if (first_digit >= 4) // if the value should sign extend to negative
		{
			val.blocks.back() |= ~((1ull << detail::highest_set_bit(val.blocks.back())) - 1); // set all bits after the high bit to one
			detail::_collapse(val);
		}
		else // otherwise should sign extend to positive
		{
			if (detail::is_neg(val)) val.blocks.push_back(0ull); // if currently negative, toss on a 0 block to make it positive
		}
	}

	return istr;
}
template<typename T, std::enable_if_t<std::is_same_v<T, detail::fixed_int_wrapper> || std::is_same_v<T, bigint>, int> = 0>
std::istream &parse_positive_dec(std::istream &istr, T &val, bool noskipws = false)
{
	std::istream::sentry sentry(istr, noskipws);
	if (!sentry) return istr;

	// start by zeroing value
	if constexpr (std::is_same_v<T, bigint>) val.blocks.clear();
	else for (std::size_t i = 0; i < val.blocks_n; ++i) val.blocks[i] = 0;

	int digit, num_digits;
	std::uint64_t block;

	std::conditional_t<std::is_same_v<T, detail::fixed_int_wrapper>, std::vector<std::uint64_t>, int> buffer; // buffer area (only used for fixed_int)
	if constexpr (std::is_same_v<T, detail::fixed_int_wrapper>)
	{
		buffer.resize(val.blocks_n);
	}
	else (void)buffer; // suppress unused warnings in the case we don't use it (it'll just be optimized away)

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
			if constexpr (std::is_same_v<T, detail::fixed_int_wrapper>)
			{
				// copy val to the buffer area and zero val
				for (std::size_t i = 0; i < val.blocks_n; ++i) buffer[i] = val.blocks[i];
				for (std::size_t i = 0; i < val.blocks_n; ++i) val.blocks[i] = 0;

				// incorporate it into val
				std::uint64_t overflow = multiply_u64_already_zero(val, { &buffer[0], buffer.size() }, scale);
				overflow += add_u64(val, block);

				// detect overflow
				if (overflow) { istr.setstate(std::ios::failbit); return istr; }
			}
			else
			{
				// incorporate it into val
				val *= scale;
				val += block;
			}
		}

		if (num_digits < 19) break;
	}

	return istr;
}
template<typename T, std::enable_if_t<std::is_same_v<T, detail::fixed_int_wrapper> || std::is_same_v<T, bigint>, int> = 0>
std::istream &parse_positive(std::istream &istr, T &val, bool noskipws = false)
{
	switch ((int)(istr.flags() & std::ios::basefield))
	{
	case (int)std::ios::hex: return parse_positive_hex(istr, val, noskipws);
	case (int)std::ios::oct: return parse_positive_oct(istr, val, noskipws);
	case (int)std::ios::dec: return parse_positive_dec(istr, val, noskipws);
	case 0:
	{
		std::istream::sentry sentry(istr, noskipws);
		if (!sentry) return istr;

		int digit;
		if ((digit = istr.peek()) == EOF) { istr.setstate(std::ios::failbit | std::ios::eofbit); return istr; }
		if (digit == '0') // if first digit is a zero we have a prefix
		{
			istr.get();
			if ((digit = istr.peek()) == EOF) return istr; // get the next character - if there is none, it's a valid decimal zero

			if (digit == 'x') { istr.get(); return parse_positive_hex(istr, val, true); }            // if second char is 'x', extract it and parse a hex value
			else if (std::isdigit((unsigned char)digit)) return parse_positive_oct(istr, val, true); // otherwise if it's a digit parse as oct
			else return istr;                                                                        // otherwise is a valid decimal zero
		}
		else return parse_positive_dec(istr, val, true); // no prefix is decimal
	}
	default: istr.setstate(std::ios::failbit); return istr; // if there are multiple base flags set, just fail (ambiguous)
	}
}

template<typename T, std::enable_if_t<std::is_same_v<T, detail::fixed_int_wrapper> || std::is_same_v<T, bigint>, int> = 0>
std::istream &parse_signed(std::istream &istr, T &val)
{
	std::istream::sentry sentry(istr);
	if (!sentry) return istr;

	int ch;
	bool neg = false;
	std::istream &(*parser)(std::istream&, T&, bool) = parse_positive; // the parser function to use (default to universal version)

	if ((istr.flags() & std::ios::basefield) == std::ios::dec) // decimal is the only base where signs come into play
	{
		if ((ch = istr.peek()) == EOF) { istr.setstate(std::ios::failbit | std::ios::eofbit); return istr; }
		if (ch == '-' || ch == '+')
		{
			istr.get();      // extract the sign char
			neg = ch == '-'; // update sign flag
		}
	}
	else if ((istr.flags() & std::ios::basefield) == 0) // otherwise if we're supposed to figure out the format and we see a sign char then it's decimal
	{
		if ((ch = istr.peek()) == EOF) { istr.setstate(std::ios::failbit | std::ios::eofbit); return istr; }
		if (ch == '-' || ch == '+')
		{
			istr.get();                  // extract the sign char
			neg = ch == '-';             // update sign flag
			parser = parse_positive_dec; // set to dec parser (rather than universal parser)
		}
	}

	if (!parser(istr, val, true)) return istr; // don't skip white space (we already did that above)

	if (neg) detail::make_neg(val); // account for sign in result
	if constexpr (!std::is_same_v<T, bigint>)
	{
		if (detail::is_neg(val) != neg) { istr.setstate(std::ios::failbit); return istr; } // check for signed overflow
	}

	return istr;
}

void detail::parse_positive_core(std::istream &istr, detail::fixed_int_wrapper val) { parse_positive(istr, val); }
void detail::parse_signed_core(std::istream &istr, detail::fixed_int_wrapper val) { parse_signed(istr, val); }

std::istream &detail::operator>>(std::istream &istr, bigint &val) { parse_signed(istr, val); return istr; }

template<bool sign, typename T>
bool _try_parse(T &res, std::string_view str, int base)
{
	std::istringstream ss(std::string(str.begin(), str.end())); // unfortunately istringstream requires a copy because stdlib is dumb

	switch (base)
	{
	case 10: break;
	case 16: ss.setf(std::ios::hex, std::ios::basefield); break;
	case 8: ss.setf(std::ios::oct, std::ios::basefield); break;
	case 0: ss.unsetf(std::ios::basefield); break;
	default: throw std::invalid_argument("unrecognized base specified");
	}

	if constexpr (sign) parse_signed(ss, res); else parse_positive(ss, res); // run the parser
	if (!ss) return false;                                                   // parse needs to succeed
	while (std::isspace((unsigned char)ss.peek())) ss.get();                 // consume trailing white space
	return ss.eof();                                                         // we need to have parsed the entire string
}

bool detail::try_parse_unsigned(fixed_int_wrapper res, std::string_view str, int base) { return _try_parse<false>(res, str, base); }
bool detail::try_parse_signed(fixed_int_wrapper res, std::string_view str, int base) { return _try_parse<true>(res, str, base); }

bool bigint::try_parse(bigint &res, std::string_view str, int base) { return _try_parse<true>(res, str, base); }
void bigint::parse(bigint &res, std::string_view str, int base) { if (!try_parse(res, str, base)) throw std::invalid_argument("failed to parse string"); }
bigint bigint::parse(std::string_view str, int base) { bigint res;  parse(res, str, base); return res; }
