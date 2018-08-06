#include <iostream>
#include <iomanip>
#include <type_traits>
#include <limits>

#include "BiggerInts.h"

int main(int argc, const char *argv[])
{
	using namespace BiggerInts;
	
	//uint_t<32> a;

	//std::uint32_t *b = &a;
	//std::cout << std::hex;
	double_int<512, false> base = 10000000000000;

	std::cout << "base: " << base << "\n\n\n";

	/*
	u64 num_high = 5;
	u64 num_low = 4;
	u64 denom = 2;

	auto dm = divmod(build_uint<128>(num_high, num_low), (uint_t<128>)denom);
	*/
	
	double_int<8192, false> big = 24;

	std::is_trivial<decltype(big)>::value;
	std::is_integral<decltype(big)>::value;

	std::is_signed<uint_t<512>>::value;
	std::is_signed<int_t<512>>::value;

	std::is_unsigned<uint_t<512>>::value;
	std::is_unsigned<int_t<512>>::value;
	
	std::make_signed_t<decltype(big)> s_big;
	std::make_unsigned_t<decltype(big)> us_big;

	std::numeric_limits<decltype(big)>::is_specialized;

	int_t<512> _big_;
	_big_.high = 1;
	_big_.low = 2;

	_big_ <<= 2;
	_big_ >>= 2;

	_big_ *= _big_;
	_big_ /= (decltype(_big_))23;

	//std::cout << std::hex;
	std::cout << "built value: " << _big_ << "\n\n";
	std::cout << "built value x-77: " << (_big_ * -(decltype(_big_))77) << "\n\n";

	//big = 23;
	//big = big & big | big;

	//bool thing = big;
	//thing = big;

	//if (big) big += 6;

	//std::cout << big << '\n';

	//uint_t<128> less_big = big;

	//big = 0xf4da218bcffe;
	//big = big;
	big = 4156475324365;
	//less_big = ++big;
	//less_big--;

	//less_big = 493;

	//big = 21;
	big *= (decltype(big))0xf4da218bcffe;
	big = big >> 5;

	int_t<8192> &as_signed = *(int_t<8192>*)&big;
	std::cout << "signed: " << as_signed << '\n';

	std::cout << "\n\n";

	std::cout << "min: " << std::numeric_limits<uint_t<512>>::min() << '\n';
	std::cout << "max: " << std::numeric_limits<uint_t<512>>::max() << '\n';
	std::cout << "low: " << std::numeric_limits<uint_t<512>>::lowest() << '\n';
	std::cout << "eps: " << std::numeric_limits<uint_t<512>>::epsilon() << '\n';

	std::cout << '\n';

	std::cout << "min: " << std::numeric_limits<int_t<512>>::min() << '\n';
	std::cout << "max: " << std::numeric_limits<int_t<512>>::max() << '\n';
	std::cout << "low: " << std::numeric_limits<int_t<512>>::lowest() << '\n';
	std::cout << "eps: " << std::numeric_limits<int_t<512>>::epsilon() << '\n';

	std::cout << "\n\n";

	decltype(big) den = 0x4d5837abcdeab75;
	//divmod(4, 5);
	auto thing = divmod(big, den);
	std::cout << big << " / " << den << " = " << thing.first << " R " << thing.second << '\n';
	std::cout << "patching: " << ((thing.first * den + thing.second) == big) << '\n';

	if (big) std::cout << "nonzero\n";



	uint_t<12> masked;

	masked = 17;
	masked += 56;
	//--masked;
	
	//++masked;
	//masked = 12;
	//masked ^= masked;
	//masked = ~masked;

	std::cout << masked << '\n';

	//bool thing = false;
	//++masked;

	//thing = masked <= masked;

	//if (thing) std::cout << "nonzero\n";
	//else std::cout << "zero\n";
	
	std::cout << std::hex;

	int_t<512> a;// , b, c, d;
	int_t<128> b, f, g;

	//asgn(a, -1);
	//asgn(b, a);

	b = -1;
	//b = a;
	a = b;
	b = (decltype(b))a;

	int i_thing;

	i_thing = (int)a;

	//asgn(b, (u64)-1);
	//asgn(c, (u64)-1);
	//d = -1;

	//asgn(e, -1);
	//asgn(f, (signed int)(-1));
	//asgn(g, (unsigned int)(-1));
	//h = -1

	std::cout << std::numeric_limits<decltype(a)>::max();

	std::cout << "\n\n--------------------------------\n\n";

	std::cout << a << '\n' << b << '\n';// << b << '\n' << c << '\n' << d << "\n\n";

	//std::cout << e << '\n' << f << '\n' << g << '\n' << h << "\n\n";



	std::cin.get();
	std::system("pause");
	return 0;
}
