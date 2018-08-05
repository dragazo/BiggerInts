#include <iostream>
#include <iomanip>

#include "BiggerInts.h"

int main(int argc, const char *argv[])
{
	using namespace BiggerInts;
	
	//uint_t<32> a;

	//std::uint32_t *b = &a;

	/*
	u64 num_high = 5;
	u64 num_low = 4;
	u64 denom = 2;

	auto dm = divmod(build_uint<128>(num_high, num_low), (uint_t<128>)denom);
	*/
	
	double_int<8192, false> big = 24;

	auto _big_ = uint_t<512>(1, 2);
	std::cout << std::hex;
	std::cout << "built value: " << _big_ << "\n\n";

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

	decltype(big) den = 0x4d5837abcdeab75;
	//divmod(4, 5);
	auto thing = divmod(big, den);
	std::cout << big << " / " << den << " = " << thing.first << " R " << thing.second << '\n';
	std::cout << "patching: " << ((thing.first * den + thing.second) == big) << '\n';

	std::cout << "plz work :3 .... " << std::setw(sizeof(big) * 2) << std::setfill('0') << std::showbase << big << '\n';
	std::cout << "width: " << std::cout.width() << '\n';
	std::cout << "bit test: " << bit_test(big, 4006) << '\n';

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
	
	std::cin.get();
	
	return 0;
}
