#include <iostream>
#include <iomanip>
#include <type_traits>
#include <limits>
#include <chrono>
#include <string>
#include <sstream>
#include <utility>
#include <cassert>

#include "BiggerInts.h"

#define COMMA ,

#define test_fail(expr, expected, name) { std::cerr << "Failed Test: " << (name) << " " expr " -> " << res << " != " << (expected) << '\n'; }
#define test(expr, expected, name) { auto res = (expr); if(res != (expected)) test_fail(expr, expected, name); }

#define t_start auto _t = std::chrono::high_resolution_clock::now();
#define t_end auto _dt = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - _t).count();
#define t_print(expr, count) { std::cout << #expr " -> " << (expr) << " : " << _dt << "ms\n"; }

#define t_test(init, expr, count) { init; t_start; for(int i = 0; i < count; ++i) (expr); t_end; t_print(expr, count); }

template<typename ...Args>
std::string tostr(Args &&...args)
{
	std::ostringstream ostr;
	(ostr << ... << std::forward<Args>(args));
	return ostr.str();
}

int main(int argc, const char *argv[])
{
	using namespace BiggerInts;
	
	int_t<512> _readme = 0;

	#define __count__ 9000
	#define __bits__ 8192 * 8
	//#define __bits__ 1024

	//std::cout << std::hex;

	t_test(uint_t<128> a = 74652465 COMMA b = 947563412 COMMA c, c = ++a, __count__);
	t_test(uint_t<128> a = 74652465 COMMA b = 947563412 COMMA c; a <<= 64; b <<= 64, c = ++a, __count__);

	t_test(int_t<__bits__> a = -(__count__ + 2) COMMA c, c = (a << 175) >> 175, __count__);
	t_test(int_t<__bits__> a = -(__count__ + 3) COMMA c, c = (a << 175) >> 175, __count__);
	t_test(int_t<__bits__> a = -(__count__ + 4) COMMA c, c = (a << 175) >> 175, __count__);
	t_test(int_t<__bits__> a = -(__count__ + 5) COMMA c, c = (a << 175) >> 175, __count__);
	
	uint_t<__bits__> big_mul = 947563412;
	assert(tostr(big_mul) == "947563412");

	big_mul *= big_mul;
	assert(tostr(big_mul) == "897876419761081744");

	big_mul *= big_mul;
	assert(tostr(big_mul) == "806182065162978263317234893050081536");

	big_mul *= big_mul;
	assert(tostr(big_mul) == "649929522190444530768966266652239714986927778859800606430979456248119296");

	big_mul += 17;
	assert(tostr(big_mul) == "649929522190444530768966266652239714986927778859800606430979456248119313");

	big_mul <<= 621;
	assert(tostr(big_mul) == "56557938587827109863786197009098171333880971620065133028223237982646344014826583116417132"
		"577306321974342507460953540245890763571403111747364084775326663465806403623083715574"
		"23307676007538288423649789830139543353463939819407429243267044228077146543261877272576");

	big_mul &= big_mul;
	assert(tostr(big_mul) == "56557938587827109863786197009098171333880971620065133028223237982646344014826583116417132"
		"577306321974342507460953540245890763571403111747364084775326663465806403623083715574"
		"23307676007538288423649789830139543353463939819407429243267044228077146543261877272576");

	big_mul |= big_mul;
	assert(tostr(big_mul) == "56557938587827109863786197009098171333880971620065133028223237982646344014826583116417132"
		"577306321974342507460953540245890763571403111747364084775326663465806403623083715574"
		"23307676007538288423649789830139543353463939819407429243267044228077146543261877272576");

	auto post_dec = big_mul--;
	assert(tostr(post_dec) == "56557938587827109863786197009098171333880971620065133028223237982646344014826583116417132"
		"577306321974342507460953540245890763571403111747364084775326663465806403623083715574"
		"23307676007538288423649789830139543353463939819407429243267044228077146543261877272576");
	assert(tostr(big_mul) == "56557938587827109863786197009098171333880971620065133028223237982646344014826583116417132"
		"577306321974342507460953540245890763571403111747364084775326663465806403623083715574"
		"23307676007538288423649789830139543353463939819407429243267044228077146543261877272575");

	auto post_inc = big_mul++;
	assert(tostr(post_inc) == "56557938587827109863786197009098171333880971620065133028223237982646344014826583116417132"
		"577306321974342507460953540245890763571403111747364084775326663465806403623083715574"
		"23307676007538288423649789830139543353463939819407429243267044228077146543261877272575");
	assert(tostr(big_mul) == "56557938587827109863786197009098171333880971620065133028223237982646344014826583116417132"
		"577306321974342507460953540245890763571403111747364084775326663465806403623083715574"
		"23307676007538288423649789830139543353463939819407429243267044228077146543261877272576");

	auto pre_inc = ++big_mul;
	assert(tostr(big_mul) == "56557938587827109863786197009098171333880971620065133028223237982646344014826583116417132"
		"577306321974342507460953540245890763571403111747364084775326663465806403623083715574"
		"23307676007538288423649789830139543353463939819407429243267044228077146543261877272577");
	assert(tostr(big_mul) == "56557938587827109863786197009098171333880971620065133028223237982646344014826583116417132"
		"577306321974342507460953540245890763571403111747364084775326663465806403623083715574"
		"23307676007538288423649789830139543353463939819407429243267044228077146543261877272577");

	auto pre_dec = --big_mul;
	assert(tostr(big_mul) == "56557938587827109863786197009098171333880971620065133028223237982646344014826583116417132"
		"577306321974342507460953540245890763571403111747364084775326663465806403623083715574"
		"23307676007538288423649789830139543353463939819407429243267044228077146543261877272576");
	assert(tostr(big_mul) == "56557938587827109863786197009098171333880971620065133028223237982646344014826583116417132"
		"577306321974342507460953540245890763571403111747364084775326663465806403623083715574"
		"23307676007538288423649789830139543353463939819407429243267044228077146543261877272576");

	big_mul >>= 191;
	assert(tostr(big_mul) == "1802039889491867190810467664715451904918476332004551151301726717494487639814495660570778575783656"
		"106649698951732305310318830222989557694740573768908873663803263415298739647066704572678996885278044454912");

	big_mul = 465;
	assert(tostr(big_mul) == "465");

	//assert(tostr(big_mul + 75u) == "540");
	tostr(big_mul ^ big_mul);

	std::cout << "outside again!\n";

	{
		t_start;

		t_test(int_t<__bits__> a = 74652465 COMMA b = 947563412 COMMA c, c = -(a ^ b), __count__);
		t_test(int_t<__bits__> a = 74652465 COMMA b = 947563412 COMMA c; a <<= 256; b <<= 256, c = -(a ^ b), __count__);
		t_test(int_t<__bits__> a = 74652465 COMMA b = 947563412 COMMA c; b <<= 256, c = -(a ^ b), __count__);
		t_test(int_t<__bits__> a = 74652465 COMMA b = 947563412 COMMA c; a <<= 256, c = -(a ^ b), __count__);

		std::cout << big_mul << "\n\n";

		t_end;
		t_print("completed", 1);
	}

	t_test(uint_t<__bits__> a = 74652465 COMMA b = 947563412 COMMA c, c = ++a, __count__);
	t_test(uint_t<__bits__> a = 74652465 COMMA b = 947563412 COMMA c; a <<= 256; b <<= 256, c = ++a, __count__);
	t_test(uint_t<__bits__> a = 74652465 COMMA b = 947563412 COMMA c; b <<= 256, c = ++a, __count__);
	t_test(uint_t<__bits__> a = 74652465 COMMA b = 947563412 COMMA c; a <<= 256, c = ++a, __count__);


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

	std::cout << "int_t<22>\n";
	std::cout << "min: " << std::numeric_limits<int_t<23>>::min() << '\n';
	std::cout << "max: " << std::numeric_limits<int_t<23>>::max() << '\n';
	std::cout << "bits: " << std::numeric_limits<int_t<23>>::digits << '\n';
	std::cout << "dig: " << std::numeric_limits<int_t<23>>::digits10 << '\n';
	std::cout << "signed: " << std::numeric_limits<int_t<23>>::is_signed << '\n';
	std::cout << "eps: " << std::numeric_limits<int_t<23>>::epsilon() << '\n';

	std::cout << "\n\n\n";

	//std::uint32_t *b = &a;
	//std::cout << std::hex;
	double_int<512, false> base = 10000000000000;

	std::cout << "base: " << base << "\n\n\n";

	
	u64 num_high = 5;
	u64 num_low = 4;
	u64 denom = 2;
	
	
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

	if (base == 10000000000000ULL) std::cout << "ALL GOOD!!!\n";
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

	std::cout << "uint_t<512>:\n";
	std::cout << "min: " << std::numeric_limits<uint_t<512>>::min() << '\n';
	std::cout << "max: " << std::numeric_limits<uint_t<512>>::max() << '\n';
	std::cout << "low: " << std::numeric_limits<uint_t<512>>::lowest() << '\n';
	std::cout << "eps: " << std::numeric_limits<uint_t<512>>::epsilon() << "\n\n";

	std::cout << "int_t<512>:\n";
	std::cout << "min: " << std::numeric_limits<int_t<512>>::min() << '\n';
	std::cout << "max: " << std::numeric_limits<int_t<512>>::max() << '\n';
	std::cout << "low: " << std::numeric_limits<int_t<512>>::lowest() << '\n';
	std::cout << "eps: " << std::numeric_limits<int_t<512>>::epsilon() << "\n\n";

	decltype(big) den = 0x4d5837abcdeab75;
	//divmod(4, 5);
	auto dmod = divmod(big, den);
	std::cout << big << " / " << den << " = " << dmod.first << " R " << dmod.second << '\n';
	std::cout << "patching: " << ((dmod.first * den + dmod.second) == big) << '\n';

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

	std::cerr << "\n\nall tests completed\n";
	return 0;
}
