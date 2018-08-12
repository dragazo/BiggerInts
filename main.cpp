#include <iostream>
#include <iomanip>
#include <type_traits>
#include <limits>
#include <chrono>

#include "BiggerInts.h"

#define COMMA ,

#define test_fail(expr, expected, name) { std::cerr << "Failed Test: " << (name) << " " expr " -> " << res << " != " << (expected) << '\n'; }
#define test(expr, expected, name) { auto res = (expr); if(res != (expected)) test_fail(expr, expected, name); }

#define t_start auto _t = std::chrono::high_resolution_clock::now();
#define t_end auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - _t).count();
#define t_print(expr, count) { std::cout << #expr " -> " << (expr) << " : " << dt << "ms\n"; }

#define t_test(init, expr, count) { init; t_start; for(int i = 0; i < count; ++i) (expr); t_end; t_print(expr, count); }

int main(int argc, const char *argv[])
{
	using namespace BiggerInts;
	
	//int_t<512> _readme = 0;

	#define __count__ 900000
	#define __bits__ 8192 * 8
	//#define __bits__ 1024

	//std::cout << std::hex;

	//t_test(uint_t<128> a = 74652465 COMMA b = 947563412 COMMA c, c = ++a, __count__);
	//t_test(uint_t<128> a = 74652465 COMMA b = 947563412 COMMA c; a <<= 64; b <<= 64, c = ++a, __count__);

	//t_test(uint_t<__bits__> a = __count__ + 2 COMMA c, c = (a << 175) >> 175, __count__);
	//t_test(uint_t<__bits__> a = __count__ + 3 COMMA c, c = (a << 175) >> 175, __count__);
	//t_test(uint_t<__bits__> a = __count__ + 4 COMMA c, c = (a << 175) >> 175, __count__);
	//t_test(uint_t<__bits__> a = __count__ + 5 COMMA c, c = (a << 175) >> 175, __count__);

	t_test(uint_t<__bits__> a = 74652465 COMMA b = 947563412 COMMA c, c = a ^ b, __count__);
	t_test(uint_t<__bits__> a = 74652465 COMMA b = 947563412 COMMA c; a <<= 256; b <<= 256, c = a ^ b, __count__);
	t_test(uint_t<__bits__> a = 74652465 COMMA b = 947563412 COMMA c; b <<= 256, c = a ^ b, __count__);
	t_test(uint_t<__bits__> a = 74652465 COMMA b = 947563412 COMMA c; a <<= 256, c = a ^ b, __count__);

	//t_test(uint_t<__bits__> a = 74652465 COMMA b = 947563412 COMMA c, c = ++a, __count__);
	//t_test(uint_t<__bits__> a = 74652465 COMMA b = 947563412 COMMA c; a <<= 256; b <<= 256, c = ++a, __count__);
	//t_test(uint_t<__bits__> a = 74652465 COMMA b = 947563412 COMMA c; b <<= 256, c = ++a, __count__);
	//t_test(uint_t<__bits__> a = 74652465 COMMA b = 947563412 COMMA c; a <<= 256, c = ++a, __count__);

	std::cin.get();
	return 0;


	/*
	std::cin >> std::setbase(0);
	std::cout << std::boolalpha;

	//std::cin >> _readme;
	std::cout << "got: " << _readme << '\n';
	std::cout << (bool)std::cin << "\n\n\n";
	//std::cin.clear();
	//std::cout << '\'' << (char)std::cin.get() << "\'\n";

	//uint_t<32> a;
	//std::cout << "offset: " << (std::size_t)(&((double_int<512, false>*)0)->low) << '\n';
	//auto off = offsetof(decltype(double_int<512, false>()), low);

	typedef int_t<23> t_t;

	//std::cout << std::boolalpha << std::hex;

	std::cout << "min: " << std::numeric_limits<t_t>::min() << '\n';
	std::cout << "max: " << std::numeric_limits<t_t>::max() << '\n';
	std::cout << "bits: " << std::numeric_limits<t_t>::digits << '\n';
	std::cout << "dig: " << std::numeric_limits<t_t>::digits10 << '\n';
	std::cout << "signed: " << std::numeric_limits<t_t>::is_signed << '\n';
	std::cout << "eps: " << std::numeric_limits<t_t>::epsilon() << '\n';

	std::cout << "\n\n\n";

	//std::uint32_t *b = &a;
	//std::cout << std::hex;
	double_int<512, false> base = 10000000000000;

	std::cout << "base: " << base << "\n\n\n";

	
	u64 num_high = 5;
	u64 num_low = 4;
	u64 denom = 2;

	auto dm = divmod(build_uint<128>(num_high, num_low), (uint_t<128>)denom);
	
	
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
	int_t<2048> _bigger_ = -46;

	_big_.high = 1;
	_big_.low = 2;

	_big_ <<= 2;
	_big_ >>= 2;

	if (base == (unsigned long long)10000000000000) std::cout << "ALL GOOD!!!\n";
	if (_big_ > _bigger_) std::cout << "less\n";
	if (_bigger_ != 65) std::cout << "less equal\n";

	//(const double &)val = 7;

	_big_ *= _big_;
	_big_ /= 23;

	_big_ *= 53;

	_big_ += 75;
	_big_ += big;
	big += _big_;
	big += 7;

	big |= 7;
	big &= 7;
	big ^= 7;

	int int_thing = 8;

	+int_thing;
	+_big_;

	//std::cout << std::hex;
	std::cout << "built value: " << _big_ << "\n\n";
	std::cout << "built value x-77: " << (_big_ % -77) << "\n\n";

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
	
	//std::cout << std::hex;

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
	*/
}
