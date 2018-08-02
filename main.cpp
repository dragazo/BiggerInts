#include <iostream>

#include "BiggerInts.h"

int main(int argc, const char *argv[])
{
	using namespace BiggerInts;

	uint_t<32> a;

	std::uint32_t *b = &a;

	

	uint_t<128> big = 24;

	big = 23;
	big = big & big | big;

	uint_t<128> less_big = big;

	big = big;

	less_big = big;

	uint_t<12> masked;

	masked = 17;
	masked += 56;
	--masked;
	
	++masked;
	masked = 12;
	//masked ^= masked;
	//masked = ~masked;

	std::cout << masked << '\n';

	bool thing = false;
	//++masked;

	thing = masked <= masked;

	if (thing) std::cout << "nonzero\n";
	else std::cout << "zero\n";

	std::cin.get();
	return 0;
}
