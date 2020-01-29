#include <cctype>
#include <algorithm>

#include "BiggerInts.h"

using namespace BiggerInts;

void detail::demote(detail::fixed_int_wrapper a, const bigint &val) noexcept
{
	if (a.blocks_n <= val.blocks.size())
	{
		for (std::size_t i = 0; i < a.blocks_n; ++i) a.blocks[i] = val.blocks[i];
	}
	else
	{
		for (std::size_t i = 0; i < val.blocks.size(); ++i) a.blocks[i] = val.blocks[i];
		std::uint64_t fill = detail::is_neg(val) ? -1 : 0;
		for (std::size_t i = val.blocks.size(); i < a.blocks_n; ++i) a.blocks[i] = fill;
	}
}

// collapses the value in the given bigint to make it satisfy the requirements of the blocks array
void collapse(bigint &a)
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
		collapse(a);
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

	collapse(val);
}

bool detail::cmp_less_non_negative(const bigint &a, const bigint &b) noexcept // same as cmp() but requires they both be non-negative
{
	if (a.blocks.size() > b.blocks.size()) return false;
	if (a.blocks.size() < b.blocks.size()) return true;
	for (std::size_t i = a.blocks.size(); i-- > 0; ) if (a.blocks[i] != b.blocks[i]) return a.blocks[i] < b.blocks[i];
	return false;
}
int detail::cmp_raw(const bigint &a, const bigint &b) noexcept
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
	else if (a.blocks.size() == 0) return val ? -1 : 0; // if a == 0
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
	collapse(*this);
}

bigint &bigint::operator++()
{
	for (std::size_t i = 1; i < blocks.size(); ++i) if (++blocks[i - 1]) { collapse(*this); return *this; }

	if (blocks.empty()) blocks.push_back(1ull);
	else
	{
		std::uint64_t high = ++blocks.back();
		if (high == 0) blocks.clear();
		else if (high == 0x8000000000000000ull) blocks.push_back(0ull);
		else collapse(*this);
	}

	return *this;
}
bigint bigint::operator++(int) { bigint cpy = *this; ++*this; return cpy; }

bigint &bigint::operator--()
{
	for (std::size_t i = 1; i < blocks.size(); ++i) if (blocks[i - 1]--) { collapse(*this); return *this; }

	if (blocks.empty()) blocks.push_back(0xffffffffffffffffull);
	else
	{
		std::uint64_t high = blocks.back()--;
		/* high == 0 case cannot happen because 0 is represented by empty array and is therefore handled above */
		if (high == 0x8000000000000000ull) blocks.push_back(0xffffffffffffffffull);
		else collapse(*this);
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

	collapse(a);
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

	collapse(*this);
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

	collapse(*this);
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

	collapse(*this);
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
		collapse(*this);
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

	collapse(res); // we need to perform one collapse operation on the finished result to put it into a valid state
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

std::pair<bigint, bigint> BiggerInts::divmod(const bigint &a, const bigint &b) { return _divmod(a, b); }
std::pair<bigint, bigint> BiggerInts::divmod(bigint &&a, const bigint &b) { return _divmod(std::move(a), b); }
std::pair<bigint, bigint> BiggerInts::divmod(const bigint &a, bigint &&b) { return _divmod(a, std::move(b)); }
std::pair<bigint, bigint> BiggerInts::divmod(bigint &&a, bigint &&b) { return _divmod(std::move(a), std::move(b)); }

bigint detail::operator/(const bigint &num, const bigint &den) { return BiggerInts::divmod(num, den).first; }
bigint detail::operator/(bigint &&num, const bigint &den) { return BiggerInts::divmod(std::move(num), den).first; }
bigint detail::operator/(const bigint &num, bigint &&den) { return BiggerInts::divmod(num, std::move(den)).first; }
bigint detail::operator/(bigint &&num, bigint &&den) { return BiggerInts::divmod(std::move(num), std::move(den)).first; }

bigint detail::operator%(const bigint &num, const bigint &den) { return BiggerInts::divmod(num, den).second; }
bigint detail::operator%(bigint &&num, const bigint &den) { return BiggerInts::divmod(std::move(num), den).second; }
bigint detail::operator%(const bigint &num, bigint &&den) { return BiggerInts::divmod(num, std::move(den)).second; }
bigint detail::operator%(bigint &&num, bigint &&den) { return BiggerInts::divmod(std::move(num), std::move(den)).second; }

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
bigint bigint::factorial(const bigint &big_v)
{
	if (detail::is_neg(big_v)) throw std::domain_error("negative value passed to bigint::factorial()");
	if (big_v.blocks.size() == 0) return 1; // 0! = 1
	if (big_v.blocks.size() > 1) throw std::domain_error("value passed to bigint::factorial() was too large"); // don't even attempt to calc in this case (very very large)

	std::uint64_t v = big_v.blocks[0];
	if (v < _factorial_lookups_count) return _factorial_lookups[v];

	bigint res = _factorial_lookups[_factorial_lookups_count - 1];
	for (bigint i = _factorial_lookups_count; i <= v; ++i) res *= i;

	return res;
}

bigint bigint::permutations(const bigint &n, const bigint &k)
{
	bigint p = n - k;
	if (detail::is_neg(k) || detail::is_neg(p)) throw std::domain_error("invalid input to bigint::permutations()");

	bigint res = 1;
	while (p < n)
	{
		++p;
		res *= p;
	}

	return res;
}

bigint bigint::combinations(const bigint &big_n, const bigint &big_k)
{
	int _cmp = cmp(big_n, big_k);
	if (detail::is_neg(big_k) || _cmp < 0) throw std::domain_error("invalid input to bigint::combinations()");
	if (big_k.blocks.size() == 0 || _cmp == 0) return 1; // C(n, 0) = C(n, n) = 1
	if (big_n.blocks.size() > 1) throw std::domain_error("value passed to bigint::combinations() was too large"); // don't even attempt to calc in this case (very very long/big computation/value)

	// from above: n > k > 0 and #n blocks = #k blocks = 1
	// convert n and k to builtin types for speed
	std::uint64_t n = big_n.blocks[0];
	std::uint64_t k = big_k.blocks[0];

	std::vector<bigint> c(k + 1); // set aside space for our dynamic programming solution. source: https://www.geeksforgeeks.org/binomial-coefficient-dp-9/
								  // k + 1 can't overflow because big_k >= 0 and has only 1 block, which is signed (so 63-bit unsigned val at most)
	c[0] = 1;
	for (std::uint64_t i = 1; i <= n; ++i)
	{
		for (std::uint64_t j = detail::min(i, k); j > 0; --j) c[j] += c[j - 1];
	}
	return c[k];
}

fibonacci_generator &fibonacci_generator::next(std::size_t n) &
{
	using std::swap;

	for (std::size_t half = n / 2; half-- > 0; )
	{
		a += b;
		b += a;
	}
	if (n & 1)
	{
		a += b;
		swap(a, b);
	}
	return *this;
}

static std::uint64_t _fib[] = {
	0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181, 6765, 10946, 17711, 28657, 46368, 75025, 121393, 196418, 317811, 514229, 832040,
	1346269, 2178309, 3524578, 5702887, 9227465, 14930352, 24157817, 39088169, 63245986, 102334155, 165580141, 267914296, 433494437, 701408733, 1134903170, 1836311903,
	2971215073, 4807526976, 7778742049, 12586269025, 20365011074, 32951280099, 53316291173, 86267571272, 139583862445, 225851433717, 365435296162, 591286729879, 956722026041,
	1548008755920, 2504730781961, 4052739537881, 6557470319842, 10610209857723, 17167680177565, 27777890035288, 44945570212853, 72723460248141, 117669030460994, 190392490709135,
	308061521170129, 498454011879264, 806515533049393, 1304969544928657, 2111485077978050, 3416454622906707, 5527939700884757, 8944394323791464, 14472334024676221, 23416728348467685,
	37889062373143906, 61305790721611591, 99194853094755497, 160500643816367088, 259695496911122585, 420196140727489673, 679891637638612258, 1100087778366101931, 1779979416004714189, 
	2880067194370816120, 4660046610375530309, 7540113804746346429, 12200160415121876738,
};
static constexpr std::size_t _fib_n = sizeof(_fib) / sizeof(*_fib);
bigint bigint::fibonacci(const bigint &big_n)
{
	if (detail::is_neg(big_n)) throw std::domain_error("argument to bigint::fibonacci() too large");
	if (big_n.blocks.size() == 0) return {};
	if (big_n.blocks.size() > 1) throw std::domain_error("argument to bigint::fibonacci() was too large");

	// get the index of the fibonacci number to generate
	std::uint64_t n = big_n.blocks[0];
	if constexpr (sizeof(std::uint64_t) != sizeof(std::size_t))
	{
		if (n != (std::size_t)n) throw std::domain_error("argument to bigint::fibonacci() was too large"); // this would result in gigabytes of bigint data anyway...
	}

	if (n < _fib_n) return _fib[n];
	else {
		fibonacci_generator g{ _fib[_fib_n - 2], _fib[_fib_n - 1] };
		g.next(n - _fib_n + 1);
		return std::move(g.b); // by returning b (instead of a like current() would do) we save 1 bigint addition (hence why line above isn't +2)
	}
}

bigint bigint::catalan(bigint n)
{
	bigint n2 = n;
	n2 <<= 1;

	bigint c = bigint::combinations(n2, n);
	++n;
	c /= n;

	return c;
}

// given a hex character, converts it to an integer [0, 15] - returns true if it was a valid hex digit. ch is only meaningful on success.
bool ext_hex(int &ch) noexcept
{
	if (ch >= '0' && ch <= '9') { ch -= '0'; return true; }
	ch |= 32;
	if (ch >= 'a' && ch <= 'f') { ch = ch - 'a' + 10; return true; }
	return false;
}

int parse_base(const std::ios_base &stream)
{
	switch ((int)(stream.flags() & std::ios::basefield))
	{
	case 0: return 0; // if no base flags are set, use prefix mode
	case (int)std::ios::dec: return 10;
	case (int)std::ios::hex: return 16;
	case (int)std::ios::oct: return 8;
	default: throw std::invalid_argument("multiple base flags were set for stream object");
	}
}

parse_fmt::parse_fmt(const std::ios_base &stream) : base(parse_base(stream)) {}
tostr_fmt::tostr_fmt(const std::ios_base &stream)
{
	switch ((int)(stream.flags() & std::ios::basefield))
	{
	case 0: // if no base flags are set, use decimal as a default
	case (int)std::ios::dec: base = 10; break;
	case (int)std::ios::hex: base = 16; break;
	case (int)std::ios::oct: base = 8; break;
	default: throw std::invalid_argument("multiple base flags were set for stream object");
	}

	showpos = stream.flags() & std::ios::showpos;
	showbase = stream.flags() & std::ios::showbase;
	uppercase = stream.flags() & std::ios::uppercase;
}

template<typename T, std::enable_if_t<std::is_same_v<T, detail::fixed_int_wrapper> || std::is_same_v<T, bigint&>, int> = 0>
std::string tostr_positive_hex(T val, const tostr_fmt &fmt, char sign_ch)
{
	constexpr bool is_bigint = std::is_same_v<T, bigint&>;

	std::string str;
	int digit, dcount;
	std::uint64_t block;

	const bool sign = detail::is_neg(val);
	const char hex_alpha = fmt.uppercase ? 'A' : 'a';

	if constexpr (is_bigint)
	{
		if (sign) val.blocks.push_back(0ull); // fix the issue of using arithmetic right shifts for bigint (negative would be infinite loop)
		str.reserve(val.blocks.size() * 16 + 4);
	}
	else
	{
		(void)sign; // suppress unused warnings
		str.reserve(val.blocks_n * 16 + 4);
	}

	while (true)
	{
		// get a block
		if constexpr (is_bigint)
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
	if constexpr (is_bigint)
	{
		if (sign)
		{
			const char f = fmt.uppercase ? 'F' : 'f';

			std::size_t s; // now we truncate all leading f's
			for (s = str.size(); s > 0; --s) if (str[s - 1] != f) break;
			str.resize(s); // chop off everything we don't want
			if (str.empty() || str.back() <= '7') str.push_back(f); // if the string is now empty or no longer sign extends to negative, add a (single) f back on
		}
		else if (str.back() >= '8') str.push_back('0'); // if val was positive but ends in an 8 or higher, add a 0 to prevent sign extending to negative on parsing
	}
	if (fmt.showbase)
	{
		str.push_back('x');
		str.push_back('0');
	}

	if (sign_ch) str.push_back(sign_ch);  // append the sign character if specified
	std::reverse(str.begin(), str.end()); // reverse the string for printing

	return str;
}
template<typename T, std::enable_if_t<std::is_same_v<T, detail::fixed_int_wrapper> || std::is_same_v<T, bigint&>, int> = 0>
std::string tostr_positive_bin(T val, const tostr_fmt &fmt, char sign_ch)
{
	constexpr bool is_bigint = std::is_same_v<T, bigint&>;

	std::string str;
	int digit, dcount;
	std::uint64_t block;

	const bool sign = detail::is_neg(val);

	if constexpr (is_bigint)
	{
		if (sign) val.blocks.push_back(0ull); // fix the issue of using arithmetic right shifts for bigint (negative would be infinite loop)
		str.reserve(val.blocks.size() * 64 + 4);
	}
	else
	{
		(void)sign; // suppress unused warnings
		str.reserve(val.blocks_n * 64 + 4);
	}

	while (true)
	{
		// get a block
		if constexpr (is_bigint)
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
			digit = block & 1;
			str.push_back('0' + digit);
			block >>= 1;
			++dcount;
		} while (block);

		if (detail::nonzero(val)) { for (; dcount < 64; ++dcount) str.push_back('0'); } // if there's still stuff, pad with zeroes and continue
		else break; // otherwise we're done
	}
	if constexpr (is_bigint)
	{
		if (sign)
		{
			std::size_t s; // now we truncate all leading 1's
			for (s = str.size(); s > 0; --s) if (str[s - 1] != '1') break;
			str.resize(s); // chop off everything we don't want
			if (str.empty() || str.back() == '0') str.push_back('1'); // if the string is now empty or no longer sign extends to negative, add a (single) 1 back on
		}
		else if (str.back() == '1') str.push_back('0'); // if val was positive but ends in a 1, add a 0 to prevent sign extending to negative on parsing
	}
	if (fmt.showbase)
	{
		str.push_back('b');
		str.push_back('0');
	}

	if (sign_ch) str.push_back(sign_ch);  // append the sign character if specified
	std::reverse(str.begin(), str.end()); // reverse the string for printing

	return str;
}
template<typename T, std::enable_if_t<std::is_same_v<T, detail::fixed_int_wrapper> || std::is_same_v<T, bigint&>, int> = 0>
std::string tostr_positive_oct(T val, const tostr_fmt &fmt, char sign_ch)
{
	constexpr bool is_bigint = std::is_same_v<T, bigint&>;

	std::string str;
	int dcount;
	std::uint64_t block;

	const bool sign = detail::is_neg(val);

	if constexpr (is_bigint)
	{
		if (sign) val.blocks.push_back(0ull); // fix the issue of using arithmetic right shifts for bigint (negative would be infinite loop)
		str.reserve(val.blocks.size() * 22 + 4);
	}
	else
	{
		(void)sign; // suppress unused warnings
		str.reserve(val.blocks_n * 22 + 4);
	}

	while (true)
	{
		// get a block
		if constexpr (is_bigint)
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
	if constexpr (is_bigint)
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
	if (fmt.showbase)
	{
		str.push_back('o');
		str.push_back('0');
	}

	if (sign_ch) str.push_back(sign_ch);  // append the sign character if specified
	std::reverse(str.begin(), str.end()); // reverse the string for printing

	return str;
}
template<typename T, std::enable_if_t<std::is_same_v<T, detail::fixed_int_wrapper> || std::is_same_v<T, bigint&>, int> = 0>
std::string tostr_positive_dec(T val, const tostr_fmt &/*fmt*/, char sign_ch)
{
	constexpr bool is_bigint = std::is_same_v<T, bigint&>;

	std::string str;
	int dcount;
	std::uint64_t block;

	typedef std::remove_reference_t<T> T_v;
	std::conditional_t<!is_bigint, std::vector<std::uint64_t>, int> buffer; // buffer area for use with fixed_int
	std::pair<T_v, T_v> temp;
	T_v base;
	if constexpr (!is_bigint)
	{
		buffer.resize(val.blocks_n * 3, 0); // allocate space for temp.first, temp.second, and base
		temp.first = { &buffer[0], val.blocks_n };
		temp.second = { &buffer[val.blocks_n], val.blocks_n };
		base = { &buffer[2 * val.blocks_n], val.blocks_n };
		base.blocks[0] = 10000000000000000000ull; // initialize base

		str.reserve(val.blocks_n * 20 + 4);
	}
	else
	{
		(void)buffer; // suppress unused warnings (will just be optimized away)
		base = 10000000000000000000ull; // initialize base

		str.reserve(val.blocks.size() * 20 + 4);
	}

	while (true)
	{
		// get a block
		if constexpr (!is_bigint)
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

	if (sign_ch) str.push_back(sign_ch);  // append the sign character if specified
	std::reverse(str.begin(), str.end()); // reverse the string for printing

	return str;
}
template<typename T, std::enable_if_t<std::is_same_v<T, detail::fixed_int_wrapper> || std::is_same_v<T, bigint&>, int> = 0>
std::string tostr_positive(T val, const tostr_fmt &fmt, char sign_ch)
{
	switch (fmt.base)
	{
	case 10: return tostr_positive_dec<T>(val, fmt, sign_ch);
	case 16: return tostr_positive_hex<T>(val, fmt, sign_ch);
	case 8: return tostr_positive_oct<T>(val, fmt, sign_ch);
	case 2: return tostr_positive_bin<T>(val, fmt, sign_ch);
	default: throw std::invalid_argument("unsupported base specified for conversion to string");
	}
}

template<typename T, std::enable_if_t<std::is_same_v<T, detail::fixed_int_wrapper> || std::is_same_v<T, bigint&>, int> = 0>
std::string tostr_signed_raw(T val, const tostr_fmt &fmt)
{
	char sign_ch = 0;
	if (fmt.base == 10) // decimal is the only base where signs come into play
	{
		if (detail::is_neg(val))
		{
			detail::make_neg(val);
			sign_ch = '-';
		}
		else if (fmt.showpos) sign_ch = '+';
	}
	return tostr_positive<T>(val, fmt, sign_ch);
}

// this is so we can reuse the stream logic for parsing values to also parse strings (without resorting to string stream)
struct string_reader
{
	const char *pos, *end;
	string_reader(std::string_view str) : pos(str.data()), end(str.data() + str.size()) {}

	inline int peek() { return pos != end ? *pos : EOF; }
	inline int get() { return pos != end ? *pos++ : EOF; }
};

void skip_whitespace(string_reader &str)
{
	while (str.pos != str.end && std::isspace((unsigned char)*str.pos)) ++str.pos;
}
void skip_whitespace(std::istream &istr)
{
	for (int ch; (ch = istr.peek()) != EOF && std::isspace((unsigned char)ch); istr.get());
}

// for sanity checking we impose this restriction on the impl templates
template<typename S, typename T>
inline constexpr bool parse_allowed = (std::is_same_v<S, std::istream&> || std::is_same_v<S, string_reader>) && (std::is_same_v<T, detail::fixed_int_wrapper> || std::is_same_v<T, bigint&>);

template<typename S, typename T, std::enable_if_t<parse_allowed<S, T>, int> = 0>
bool try_parse_positive_hex(S src, T val, int /*base*/, bool skipws = true)
{
	constexpr bool is_string = std::is_same_v<S, string_reader>;
	constexpr bool is_bigint = std::is_same_v<T, bigint&>;

	if (skipws) skip_whitespace(src);

	// start by zeroing value
	if constexpr (is_bigint) val.blocks.clear();
	else for (std::size_t i = 0; i < val.blocks_n; ++i) val.blocks[i] = 0;

	int first_digit, digit, num_digits;
	std::uint64_t block;

	// first char must be a hex digit - don't extract
	if ((first_digit = src.peek()) == EOF) return false;
	if (!ext_hex(first_digit)) return false;

	while (true)
	{
		block = 0; // read a block of 16 hex digits
		for (num_digits = 0; num_digits < 16; ++num_digits)
		{
			// get the digit and make sure it's in range - only extract it if it's good
			if ((digit = src.peek()) == EOF) break;
			if (!ext_hex(digit)) break;
			src.get();

			block <<= 4;
			block |= digit;
		}
		if (num_digits != 0)
		{
			if constexpr (!is_bigint)
			{
				// detect overflow
				if (val.blocks[val.blocks_n - 1] >> (std::uint64_t)(64 - 4 * num_digits))
				{
					if constexpr (!is_string)
					{
						for (int ch; (ch = src.peek()) != EOF && ext_hex(ch); src.get()); // if reading from a stream consume extra valid chars on overflow failure
					}
					return false;
				}

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

	if constexpr (is_bigint)
	{
		static_assert(detail::highest_set_bit(0) == 0); // sanity check for the behavior below

		if (first_digit >= 8) // if the value should sign extend to negative
		{
			val.blocks.back() |= ~((1ull << detail::highest_set_bit(val.blocks.back())) - 1); // set all bits after the high bit to one
			collapse(val);
		}
		else // otherwise should sign extend to positive
		{
			if (detail::is_neg(val)) val.blocks.push_back(0ull); // if currently negative, toss on a 0 block to make it positive
		}
	}

	if constexpr (is_string)
	{
		skip_whitespace(src);
		return src.pos == src.end; // we need to have processed the entire string
	}
	else return true;
}
template<typename S, typename T, std::enable_if_t<parse_allowed<S, T>, int> = 0>
bool try_parse_positive_bin(S src, T val, int /*base*/, bool skipws = true)
{
	constexpr bool is_string = std::is_same_v<S, string_reader>;
	constexpr bool is_bigint = std::is_same_v<T, bigint&>;

	if (skipws) skip_whitespace(src);

	// start by zeroing value
	if constexpr (is_bigint) val.blocks.clear();
	else for (std::size_t i = 0; i < val.blocks_n; ++i) val.blocks[i] = 0;

	int first_digit, digit, num_digits;
	std::uint64_t block;

	// first char must be a hex digit - don't extract
	if ((first_digit = src.peek()) == EOF) return false;
	if (first_digit < '0' || first_digit > '1') return false;
	first_digit -= '0';

	while (true)
	{
		block = 0; // read a block of 64 bin digits
		for (num_digits = 0; num_digits < 64; ++num_digits)
		{
			// get the digit and make sure it's in range - only extract it if it's good
			if ((digit = src.peek()) == EOF) break;
			if (digit < '0' || digit > '1') break;
			digit -= '0';
			src.get();

			block <<= 1;
			block |= digit;
		}
		if (num_digits != 0)
		{
			if constexpr (!is_bigint)
			{
				// detect overflow
				if (val.blocks[val.blocks_n - 1] >> (std::uint64_t)(64 - num_digits))
				{
					if constexpr (!is_string)
					{
						for (int ch; (ch = src.peek()) != EOF && ch >= '0' && ch <= '1'; src.get()); // if reading from a stream consume extra valid chars on overflow failure
					}
					return false;
				}

				// incorporate it into value
				shl(val, (std::uint64_t)(num_digits));
				val.blocks[0] |= block;
			}
			else
			{
				// incorporate it into value
				val <<= (std::uint64_t)(num_digits);
				val |= block;
			}
		}

		if (num_digits < 64) break;
	}

	if constexpr (is_bigint)
	{
		static_assert(detail::highest_set_bit(0) == 0); // sanity check for the behavior below

		if (first_digit == 1) // if the value should sign extend to negative
		{
			val.blocks.back() |= ~((1ull << detail::highest_set_bit(val.blocks.back())) - 1); // set all bits after the high bit to one
			collapse(val);
		}
		else // otherwise should sign extend to positive
		{
			if (detail::is_neg(val)) val.blocks.push_back(0ull); // if currently negative, toss on a 0 block to make it positive
		}
	}

	if constexpr (is_string)
	{
		skip_whitespace(src);
		return src.pos == src.end; // we need to have processed the entire string
	}
	else return true;
}
template<typename S, typename T, std::enable_if_t<parse_allowed<S, T>, int> = 0>
bool try_parse_positive_oct(S src, T val, int /*base*/, bool skipws = true)
{
	constexpr bool is_string = std::is_same_v<S, string_reader>;
	constexpr bool is_bigint = std::is_same_v<T, bigint&>;

	if (skipws) skip_whitespace(src);

	// start by zeroing value
	if constexpr (is_bigint) val.blocks.clear();
	else for (std::size_t i = 0; i < val.blocks_n; ++i) val.blocks[i] = 0;

	int first_digit, digit, num_digits;
	std::uint64_t block;

	// first char must be an oct digit - don't extract
	if ((first_digit = src.peek()) == EOF) return false;
	if (first_digit < '0' || first_digit > '7') return false;
	first_digit -= '0';

	while (true)
	{
		block = 0; // read a block of 21 oct digits
		for (num_digits = 0; num_digits < 21; ++num_digits)
		{
			// get the digit and make sure it's in range - only extract it if it's good
			if ((digit = src.peek()) == EOF) break;
			if (digit < '0' || digit > '7') break;
			src.get();

			block <<= 3;
			block |= digit - '0';
		}
		if (num_digits != 0)
		{
			if constexpr (!is_bigint)
			{
				// detect overflow
				if (val.blocks[val.blocks_n - 1] >> (std::uint64_t)(64 - 3 * num_digits))
				{
					if constexpr (!is_string)
					{
						for (int ch; (ch = src.peek()) != EOF && ch >= '0' && ch <= '7'; src.get()); // if reading from a stream consume extra valid chars on overflow failure
					}
					return false;
				}

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

	if constexpr (is_bigint)
	{
		static_assert(detail::highest_set_bit(0) == 0); // sanity check for the behavior below

		if (first_digit >= 4) // if the value should sign extend to negative
		{
			val.blocks.back() |= ~((1ull << detail::highest_set_bit(val.blocks.back())) - 1); // set all bits after the high bit to one
			collapse(val);
		}
		else // otherwise should sign extend to positive
		{
			if (detail::is_neg(val)) val.blocks.push_back(0ull); // if currently negative, toss on a 0 block to make it positive
		}
	}

	if constexpr (is_string)
	{
		skip_whitespace(src);
		return src.pos == src.end; // we need to have processed the entire string
	}
	else return true;
}
template<typename S, typename T, std::enable_if_t<parse_allowed<S, T>, int> = 0>
bool try_parse_positive_dec(S src, T val, int /*base*/, bool skipws = true)
{
	constexpr bool is_string = std::is_same_v<S, string_reader>;
	constexpr bool is_bigint = std::is_same_v<T, bigint&>;

	if (skipws) skip_whitespace(src);

	// start by zeroing value
	if constexpr (is_bigint) val.blocks.clear();
	else for (std::size_t i = 0; i < val.blocks_n; ++i) val.blocks[i] = 0;

	int digit, num_digits;
	std::uint64_t block;

	std::conditional_t<!is_bigint, std::vector<std::uint64_t>, int> buffer; // buffer area (only used for fixed_int)
	if constexpr (!is_bigint)
	{
		buffer.resize(val.blocks_n);
	}
	else (void)buffer; // suppress unused warnings in the case we don't use it (it'll just be optimized away)

	// first char must be a dec digit - don't extract
	if ((digit = src.peek()) == EOF) return false;
	if (digit < '0' || digit > '9') return false;

	while (true)
	{
		block = 0;               // read a block of 19 dec digits
		std::uint64_t scale = 1; // scale factor to apply to val for block incorporation (see below)
		for (num_digits = 0; num_digits < 19; ++num_digits)
		{
			// get the digit and make sure it's in range - only extract it if it's good
			if ((digit = src.peek()) == EOF) break;
			if (digit < '0' || digit > '9') break;
			src.get();

			scale *= 10;
			block *= 10;
			block += digit - '0';
		}
		if (num_digits != 0)
		{
			if constexpr (!is_bigint)
			{
				// copy val to the buffer area and zero val
				for (std::size_t i = 0; i < val.blocks_n; ++i) buffer[i] = val.blocks[i];
				for (std::size_t i = 0; i < val.blocks_n; ++i) val.blocks[i] = 0;

				// incorporate it into val
				std::uint64_t overflow = multiply_u64_already_zero(val, { &buffer[0], buffer.size() }, scale);
				overflow += add_u64(val, block);

				// detect overflow
				if (overflow)
				{
					if constexpr (!is_string)
					{
						for (int ch; (ch = src.peek()) != EOF && ch >= '0' && ch <= '9'; src.get()); // if reading from a stream consume extra valid chars on overflow failure
					}
					return false;
				}
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

	if constexpr (is_string)
	{
		skip_whitespace(src);
		return src.pos == src.end; // we need to have processed the entire string
	}
	else return true;
}
template<typename S, typename T, std::enable_if_t<parse_allowed<S, T>, int> = 0>
bool try_parse_positive(S src, T val, int base, bool skipws = true)
{
	constexpr bool is_string = std::is_same_v<S, string_reader>;

	switch (base)
	{
	case 16: return try_parse_positive_hex<S, T>(src, val, base, skipws);
	case 10: return try_parse_positive_dec<S, T>(src, val, base, skipws);
	case 8: return try_parse_positive_oct<S, T>(src, val, base, skipws);
	case 2: return try_parse_positive_bin<S, T>(src, val, base, skipws);
	case 0:
	{
		if (skipws) skip_whitespace(src);

		int digit;
		if ((digit = src.peek()) == EOF) return false;
		if (digit == '0') // if first digit is a zero we have a prefix
		{
			src.get();
			if ((digit = src.peek()) == EOF) return true; // get the next character - if there is none, it's a valid decimal zero
			digit = std::tolower((unsigned char)digit);

			if (digit == 'x') { src.get(); return try_parse_positive_hex<S, T>(src, val, base, false); }             // if second char is 'x', extract it and parse a hex value
			else if (digit == 'o') { src.get(); return try_parse_positive_oct<S, T>(src, val, base, false); }        // if second char is 'o', extract it and parse an oct value
			else if (digit == 'b') { src.get(); return try_parse_positive_bin<S, T>(src, val, base, false); }        // if second char is 'b', extract it and parse a bin value
			else if (std::isdigit((unsigned char)digit)) return try_parse_positive_dec<S, T>(src, val, base, false); // if second char is a digit, parse as dec (no prefix)
			else // otherwise invalid char
			{
				if constexpr (is_string) // if we're in string mode, we need to consume the entire string
				{
					skip_whitespace(src);
					return src.pos == src.end; // otherwise if the rest of string is white space, is valid, otherwise contains illegal chars
				}
				else return true; // otherwise we're good to go
			}
		}
		else return try_parse_positive_dec<S, T>(src, val, base, false); // no prefix is decimal
	}
	default: throw std::invalid_argument("unrecognized parse base specified");
	}
}

template<typename S, typename T, std::enable_if_t<parse_allowed<S, T>, int> = 0>
bool try_parse_signed(S src, T val, int base)
{
	constexpr bool is_bigint = std::is_same_v<T, bigint&>;

	skip_whitespace(src);

	int ch;
	bool neg = false;
	bool (*parser)(S, T, int, bool) = try_parse_positive; // the parser function to use (default to universal version)

	if (base == 10) // decimal is the only base where signs come into play
	{
		if ((ch = src.peek()) == EOF) return false;
		if (ch == '-' || ch == '+')
		{
			src.get();       // extract the sign char
			neg = ch == '-'; // update sign flag
		}
	}
	else if (base == 0) // otherwise if we're supposed to figure out the format and we see a sign char then it's decimal
	{
		if ((ch = src.peek()) == EOF) return false;
		if (ch == '-' || ch == '+')
		{
			src.get();                       // extract the sign char
			neg = ch == '-';                 // update sign flag
			parser = try_parse_positive_dec; // set to dec parser (rather than universal parser)
		}
	}

	if (!parser(src, val, base, false)) return false; // don't skip white space (we already did that above)

	if (neg) detail::make_neg(val); // account for sign in result
	if constexpr (!is_bigint)
	{
		if (detail::is_neg(val) != neg) return false; // check for signed overflow
	}

	return true;
}

std::string detail::tostr_unsigned(detail::fixed_int_wrapper val, const tostr_fmt &fmt) { return tostr_positive<detail::fixed_int_wrapper>(val, fmt, 0); }
std::string detail::tostr_signed(detail::fixed_int_wrapper val, const tostr_fmt &fmt) { return tostr_signed_raw<detail::fixed_int_wrapper>(val, fmt); }

std::string tostr_fmt::operator()(const bigint &val) const
{
	auto cpy = val;
	return tostr_signed_raw<bigint&>(cpy, *this);
}
std::string tostr_fmt::operator()(bigint &&val) const { return tostr_signed_raw<bigint&>(val, *this); }

std::ostream &detail::print_unsigned(std::ostream &ostr, detail::fixed_int_wrapper val) { return ostr << detail::tostr_unsigned(val, tostr_fmt{ ostr }); }
std::ostream &detail::print_signed(std::ostream &ostr, detail::fixed_int_wrapper val) { return ostr << detail::tostr_signed(val, tostr_fmt{ ostr }); }

std::ostream &detail::operator<<(std::ostream &ostr, const bigint &val) { return ostr << tostr_fmt{ ostr }(val); }
std::ostream &detail::operator<<(std::ostream &ostr, bigint &&val) { return ostr << tostr_fmt{ ostr }(std::move(val)); }

bool detail::try_parse_unsigned(detail::fixed_int_wrapper val, std::string_view str, int base) { return ::try_parse_positive<string_reader, fixed_int_wrapper>(str, val, base); }
bool detail::try_parse_signed(detail::fixed_int_wrapper val, std::string_view str, int base) { return ::try_parse_signed<string_reader, fixed_int_wrapper>(str, val, base); }

std::istream &detail::extract_unsigned(detail::fixed_int_wrapper res, std::istream &istr, int base)
{
	if (!::try_parse_positive<std::istream&, fixed_int_wrapper>(istr, res, base)) istr.setstate(std::ios::failbit);
	return istr;
}
std::istream &detail::extract_signed(detail::fixed_int_wrapper res, std::istream &istr, int base)
{
	if (!::try_parse_signed<std::istream&, fixed_int_wrapper>(istr, res, base)) istr.setstate(std::ios::failbit);
	return istr;
}

bool parse_fmt::operator()(bigint &res, std::string_view str) const { return ::try_parse_signed<string_reader, bigint&>(str, res, base); }
std::istream &parse_fmt::operator()(bigint &res, std::istream &istr) const { if (!::try_parse_signed<std::istream&, bigint&>(istr, res, parse_base(istr))) istr.setstate(std::ios::failbit); return istr; }

bigint parse_fmt::operator()(std::string_view str) const { detail::bigint res; if (!operator()(res, str)) throw std::invalid_argument("failed to parse string"); return res; }