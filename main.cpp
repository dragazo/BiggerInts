#include <iostream>
#include <iomanip>

#include "BiggerInts.h"

int main(int argc, const char *argv[])
{
	using namespace BiggerInts;

	//uint_t<32> a;

	//std::uint32_t *b = &a;

	

	uint_t<8192> big = 24;



	//big = 23;
	//big = big & big | big;

	//bool thing = big;
	//thing = big;

	//if (big) big += 6;

	//std::cout << big << '\n';

	//uint_t<128> less_big = big;

	//big = big;
	big = 4156475324365;
	big <<= 4000;
	//less_big = ++big;
	//less_big--;

	//less_big = 493;

	std::cout << "plz work :3 .... " << std::setw(sizeof(big) * 2) << std::setfill('0') << std::showbase << big << '\n';
	std::cout << "width: " << std::cout.width() << '\n';

	uint_t<12> masked;

	masked = 17;
	masked += 56;
	//--masked;
	
	//++masked;
	//masked = 12;
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
