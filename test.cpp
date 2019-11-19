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

	assert(tostr(uint_t<128>(2) * uint_t<128>(13)) == "26");
	assert(tostr(uint_t<128>(2) + uint_t<128>(13)) == "15");
	assert(tostr(uint_t<128>(2) / uint_t<128>(13)) == "0");
	assert(tostr(uint_t<128>(2) % uint_t<128>(13)) == "2");

	assert(tostr(uint_t<128>(52) * uint_t<128>(3)) == "156");
	assert(tostr(uint_t<128>(52) + uint_t<128>(3)) == "55");
	assert(tostr(uint_t<128>(52) / uint_t<128>(3)) == "17");
	assert(tostr(uint_t<128>(52) % uint_t<128>(3)) == "1");

	assert(tostr(int_t<128>(2) * int_t<128>(-13)) == "-26");
	assert(tostr(int_t<128>(2) + int_t<128>(-13)) == "-11");
	assert(tostr(int_t<128>(2) / int_t<128>(-13)) == "0");
	assert(tostr(int_t<128>(2) % int_t<128>(13)) == "2");

	assert(tostr(int_t<128>(2) & int_t<128>(-13)) == "2");
	assert(tostr(int_t<128>(2) | int_t<128>(-13)) == "-13");
	assert(tostr(int_t<128>(2) ^ int_t<128>(-13)) == "-15");

	assert(tostr(int_t<128>(-13) << 1) == "-26");
	assert(tostr(int_t<128>(-13) << 2) == "-52");

	assert(tostr(int_t<128>(-13) >> 1) == "-7");
	assert(tostr(int_t<128>(-13) >> 2) == "-4");

	{
		static_assert(std::is_same<decltype(uint_t<128>::high), std::uint64_t>::value);
		static_assert(std::is_same<decltype(uint_t<128>::low), std::uint64_t>::value);

		static_assert(std::is_same<decltype(int_t<128>::high), std::uint64_t>::value);
		static_assert(std::is_same<decltype(int_t<128>::low), std::uint64_t>::value);
	}

	uint_t<8192> big_mul = 947563412;
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
	assert(tostr(pre_inc) == "56557938587827109863786197009098171333880971620065133028223237982646344014826583116417132"
		"577306321974342507460953540245890763571403111747364084775326663465806403623083715574"
		"23307676007538288423649789830139543353463939819407429243267044228077146543261877272577");
	assert(tostr(big_mul) == "56557938587827109863786197009098171333880971620065133028223237982646344014826583116417132"
		"577306321974342507460953540245890763571403111747364084775326663465806403623083715574"
		"23307676007538288423649789830139543353463939819407429243267044228077146543261877272577");
	assert(pre_inc == big_mul);
	assert(pre_inc <= big_mul);
	assert(pre_inc >= big_mul);

	auto pre_dec = --big_mul;
	assert(tostr(pre_dec) == "56557938587827109863786197009098171333880971620065133028223237982646344014826583116417132"
		"577306321974342507460953540245890763571403111747364084775326663465806403623083715574"
		"23307676007538288423649789830139543353463939819407429243267044228077146543261877272576");
	assert(tostr(big_mul) == "56557938587827109863786197009098171333880971620065133028223237982646344014826583116417132"
		"577306321974342507460953540245890763571403111747364084775326663465806403623083715574"
		"23307676007538288423649789830139543353463939819407429243267044228077146543261877272576");
	assert(pre_dec == big_mul);
	assert(pre_dec <= big_mul);
	assert(pre_dec >= big_mul);

	big_mul >>= 191;
	assert(tostr(big_mul) == "1802039889491867190810467664715451904918476332004551151301726717494487639814495660570778575783656"
		"106649698951732305310318830222989557694740573768908873663803263415298739647066704572678996885278044454912");

	big_mul = 465;
	assert(tostr(big_mul) == "465");

	assert(tostr(big_mul + 75u) == "540");
	assert(tostr(big_mul ^ big_mul) == "0");

	// -- type equivalence tests -- //

	static_assert(std::is_same<uint_t<8>, std::uint8_t>::value, "type equivalence violation");
	static_assert(std::is_same<uint_t<16>, std::uint16_t>::value, "type equivalence violation");
	static_assert(std::is_same<uint_t<32>, std::uint32_t>::value, "type equivalence violation");
	static_assert(std::is_same<uint_t<64>, std::uint64_t>::value, "type equivalence violation");

	static_assert(std::is_same<int_t<8>, std::int8_t>::value, "type equivalence violation");
	static_assert(std::is_same<int_t<16>, std::int16_t>::value, "type equivalence violation");
	static_assert(std::is_same<int_t<32>, std::int32_t>::value, "type equivalence violation");
	static_assert(std::is_same<int_t<64>, std::int64_t>::value, "type equivalence violation");

	// -- numeric limits tests -- //

	assert(tostr(std::numeric_limits<uint_t<512>>::min()) == "0");
	assert(tostr(std::numeric_limits<uint_t<512>>::max()) == "13407807929942597099574024998205846127479365820592393377723561443721764030073546976801874298166903427690031858186486050853753882811946569946433649006084095");
	assert(tostr(std::numeric_limits<uint_t<512>>::lowest()) == "0");

	assert(tostr(std::numeric_limits<int_t<512>>::min()) == "-6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042048");
	assert(tostr(std::numeric_limits<int_t<512>>::max()) == "6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042047");
	assert(tostr(std::numeric_limits<int_t<512>>::lowest()) == "-6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042048");

	assert(tostr(std::numeric_limits<uint_t<128>>::min()) == "0");
	assert(tostr(std::numeric_limits<uint_t<128>>::max()) == "340282366920938463463374607431768211455");
	assert(tostr(std::numeric_limits<uint_t<128>>::lowest()) == "0");

	assert(tostr(std::numeric_limits<int_t<128>>::min()) == "-170141183460469231731687303715884105728");
	assert(tostr(std::numeric_limits<int_t<128>>::max()) == "170141183460469231731687303715884105727");
	assert(tostr(std::numeric_limits<int_t<128>>::lowest()) == "-170141183460469231731687303715884105728");

	assert(tostr(std::numeric_limits<uint_t<64>>::min()) == "0");
	assert(tostr(std::numeric_limits<uint_t<64>>::max()) == "18446744073709551615");
	assert(tostr(std::numeric_limits<uint_t<64>>::lowest()) == "0");

	assert(tostr(std::numeric_limits<int_t<64>>::min()) == "-9223372036854775808");
	assert(tostr(std::numeric_limits<int_t<64>>::max()) == "9223372036854775807");
	assert(tostr(std::numeric_limits<int_t<64>>::lowest()) == "-9223372036854775808");

	// -- all tests completed -- //

	std::cerr << "\n\nall tests completed successfully\n";
	return 0;
}
