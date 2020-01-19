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

using namespace BiggerInts;

#define COMMA ,

#define assert_throws(except, expr) try { expr; std::cerr << "didn't throw\n"; assert(false); } catch (const except&) {} catch (...) { std::cerr << "threw wrong type\n"; assert(false); }

#define test_fail(expr, expected, name) { std::cerr << "Failed Test: " << (name) << " " expr " -> " << res << " != " << (expected) << '\n'; }
#define test(expr, expected, name) { auto res = (expr); if(res != (expected)) test_fail(expr, expected, name); }

#define t_start auto _t = std::chrono::high_resolution_clock::now()
#define t_end(name, cnt) { auto _dt = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - _t).count(); std::cerr << name ": " << _dt << " ns (" << (cnt) << " iterations) avg: " << (_dt / (cnt)) << " ns\n"; }

#define t_test(init, expr, count) { init; t_start; for(int i = 0; i < count; ++i) (expr); t_end; t_print(expr, count); }

template<typename ...Args>
std::string sstostr(const Args &...args)
{
	std::ostringstream ostr;
	(ostr << ... << args);
	return ostr.str();
}

constexpr bool benchmark_do_uint = true;
constexpr bool benchmark_do_int = true;
constexpr bool benchmark_do_bigint = true;

volatile std::uint64_t dont_optimize_this_away = 0; // Kelek's breath, optimizers can really make benchmarking annoying

void do_something(const BiggerInts::bigint &v) { dont_optimize_this_away += v.blocks.empty() ? 0ull : v.blocks[0]; }
template<std::uint64_t bits, bool sign>
void do_something(const BiggerInts::detail::fixed_int<bits, sign> &v) { dont_optimize_this_away += v.blocks[0]; }
template<typename T>
void do_something(const std::pair<T, T> &v) { do_something(v.first); do_something(v.second); }
void do_something(const std::string &a) { dont_optimize_this_away += (std::uint64_t)a.size(); }

template<typename BinaryFunction>
void benchmark_binary(const char *name, std::size_t count_1, std::size_t count_2, std::size_t count_3, std::size_t count_4, BinaryFunction f)
{
	if constexpr (benchmark_do_uint)
	{
		uint_t<128> vals[] =
		{
			uint_t<128>::parse("164953864532146553885641326541323564132"),
			uint_t<128>::parse("68573958787958483948687979685743988795"),
			uint_t<128>::parse("47785319425432789132452387667540040594"),
			uint_t<128>::parse("62641345234867594511231243576678624112"),
			uint_t<128>::parse("12642723264275675264211264237286172641"),
			uint_t<128>::parse("41672765971459452836324591059442640942"),
			uint_t<128>::parse("46176276219426832272141594256267216412"),
			uint_t<128>::parse("16126745353512798654214787864567825452"),
			uint_t<128>::parse("13786541237897643215487675429731824613"),
			uint_t<128>::parse("45389783215348783514264316267678795642"),
			uint_t<128>::parse("26461343345768768645949213234357678665"),
			uint_t<128>::parse("45676597896543534326162612343537367866"),
			uint_t<128>::parse("16435738689456523313264567898883765642"),
			uint_t<128>::parse("46567383388656261334357573838366565421"),
		};

		uint_t<128> val2 = uint_t<128>::parse("165798541236525214");
		decltype(f(vals[0], vals[0])) dd;
		
		t_start;
		for (std::size_t i = 0; i < count_1; ++i)
		{
			for (const auto &v : vals)
			{
				dd = f(v, val2); do_something(dd);
				dd = f(val2, v); do_something(dd);
				dd = f(v, v); do_something(dd);
				dd = f(val2, val2); do_something(dd);
			}
		}
		t_end("1) " << name << " unsigned 128", 4 * count_1);
	}
	
	if constexpr (benchmark_do_int)
	{
		int_t<128> vals[] =
		{
			int_t<128>::parse("164953864532146553885641326541323564132"),
			int_t<128>::parse("68573958787958483948687979685743988795"),
			int_t<128>::parse("47785319425432789132452387667540040594"),
			int_t<128>::parse("62641345234867594511231243576678624112"),
			int_t<128>::parse("12642723264275675264211264237286172641"),
			int_t<128>::parse("41672765971459452836324591059442640942"),
			int_t<128>::parse("46176276219426832272141594256267216412"),
			int_t<128>::parse("16126745353512798654214787864567825452"),
			int_t<128>::parse("13786541237897643215487675429731824613"),
			int_t<128>::parse("45389783215348783514264316267678795642"),
			int_t<128>::parse("26461343345768768645949213234357678665"),
			int_t<128>::parse("45676597896543534326162612343537367866"),
			int_t<128>::parse("16435738689456523313264567898883765642"),
			int_t<128>::parse("46567383388656261334357573838366565421"),
		};

		int_t<128> val2 = int_t<128>::parse("165798541236525214");
		decltype(f(vals[0], vals[0])) dd;

		t_start;
		for (std::size_t i = 0; i < count_1; ++i)
		{
			for (const auto &v : vals)
			{
				dd = f(v, val2); do_something(dd);
				dd = f(val2, v); do_something(dd);
				dd = f(v, v); do_something(dd);
				dd = f(val2, val2); do_something(dd);
			}
		}
		t_end("1) " << name << " signed 128", 4 * count_1);
	}

	if constexpr (benchmark_do_bigint)
	{
		bigint vals[] =
		{
			bigint::parse("164953864532146553885641326541323564132"),
			bigint::parse("68573958787958483948687979685743988795"),
			bigint::parse("47785319425432789132452387667540040594"),
			bigint::parse("62641345234867594511231243576678624112"),
			bigint::parse("12642723264275675264211264237286172641"),
			bigint::parse("41672765971459452836324591059442640942"),
			bigint::parse("46176276219426832272141594256267216412"),
			bigint::parse("16126745353512798654214787864567825452"),
			bigint::parse("13786541237897643215487675429731824613"),
			bigint::parse("45389783215348783514264316267678795642"),
			bigint::parse("26461343345768768645949213234357678665"),
			bigint::parse("45676597896543534326162612343537367866"),
			bigint::parse("16435738689456523313264567898883765642"),
			bigint::parse("46567383388656261334357573838366565421"),
		};

		bigint val2 = bigint::parse("165798541236525214");
		decltype(f(vals[0], vals[0])) dd;

		t_start;
		for (std::size_t i = 0; i < count_1; ++i)
		{
			for (const auto &v : vals)
			{
				dd = f(v, val2); do_something(dd);
				dd = f(val2, v); do_something(dd);
				dd = f(v, v); do_something(dd);
				dd = f(val2, val2); do_something(dd);
			}
		}
		t_end("1) " << name << " bigint", 4 * count_1);
	}
	
	if constexpr (benchmark_do_uint)
	{
		uint_t<256> vals[] =
		{
			uint_t<256>::parse("1649538645321465538856413265789453212354678989785653212354878653424132356413"),
			uint_t<256>::parse("6857395878795848394868797968574345687978925312545679532123453789786514908879"),
			uint_t<256>::parse("4778531942543278913245238766235456868653123123456789899868652525252175004059"),
			uint_t<256>::parse("6264134523486759565378935212111104040045605451123124323546474868757667862411"),
			uint_t<256>::parse("1264272326427560643579865532334260603130262002462678317526421126423728617264"),
			uint_t<256>::parse("4167276597145945283632459105940634353766783723761163234768566090999144264094"),
			uint_t<256>::parse("4617627621942683227214159425626700616234303760377685879004613264267270421641"),
			uint_t<256>::parse("1612674535351270645679076866462442652656767925402042569865421478786456782545"),
			uint_t<256>::parse("1378654123789764321548066643503435773738979646034503430343045767542973182461"),
			uint_t<256>::parse("4538978321534878351406064967605678378337866654650644303133226431626767879564"),
			uint_t<256>::parse("2646134334576876864594921324653303750264246256728549413640672653423435767866"),
			uint_t<256>::parse("4567659789654353432616261234350664572678024624362312312645273860795403736786"),
			uint_t<256>::parse("1643573868945652331326456064672738264260264256425372867294025645789888376564"),
			uint_t<256>::parse("4656738338865626133435757383064650762860724950414629645642567597260836656542"),
		};

		uint_t<256> val2 = uint_t<256>::parse("16579862451765276262462456641541236525214");
		decltype(f(vals[0], vals[0])) dd;

		t_start;
		for (std::size_t i = 0; i < count_2; ++i)
		{
			for (const auto &v : vals)
			{
				dd = f(v, val2); do_something(dd);
				dd = f(val2, v); do_something(dd);
				dd = f(v, v); do_something(dd);
				dd = f(val2, val2); do_something(dd);
			}
		}
		t_end("2) " << name << " unsigned 256", 4 * count_2);
	}

	if constexpr (benchmark_do_int)
	{
		int_t<256> vals[] =
		{
			int_t<256>::parse("1649538645321465538856413265789453212354678989785653212354878653424132356413"),
			int_t<256>::parse("6857395878795848394868797968574345687978925312545679532123453789786514908879"),
			int_t<256>::parse("4778531942543278913245238766235456868653123123456789899868652525252175004059"),
			int_t<256>::parse("6264134523486759565378935212111104040045605451123124323546474868757667862411"),
			int_t<256>::parse("1264272326427560643579865532334260603130262002462678317526421126423728617264"),
			int_t<256>::parse("4167276597145945283632459105940634353766783723761163234768566090999144264094"),
			int_t<256>::parse("4617627621942683227214159425626700616234303760377685879004613264267270421641"),
			int_t<256>::parse("1612674535351270645679076866462442652656767925402042569865421478786456782545"),
			int_t<256>::parse("1378654123789764321548066643503435773738979646034503430343045767542973182461"),
			int_t<256>::parse("4538978321534878351406064967605678378337866654650644303133226431626767879564"),
			int_t<256>::parse("2646134334576876864594921324653303750264246256728549413640672653423435767866"),
			int_t<256>::parse("4567659789654353432616261234350664572678024624362312312645273860795403736786"),
			int_t<256>::parse("1643573868945652331326456064672738264260264256425372867294025645789888376564"),
			int_t<256>::parse("4656738338865626133435757383064650762860724950414629645642567597260836656542"),
		};

		int_t<256> val2 = int_t<256>::parse("16579862451765276262462456641541236525214");
		decltype(f(vals[0], vals[0])) dd;

		t_start;
		for (std::size_t i = 0; i < count_2; ++i)
		{
			for (const auto &v : vals)
			{
				dd = f(v, val2); do_something(dd);
				dd = f(val2, v); do_something(dd);
				dd = f(v, v); do_something(dd);
				dd = f(val2, val2); do_something(dd);
			}
		}
		t_end("2) " << name << " signed 256", 4 * count_2);
	}

	if constexpr (benchmark_do_uint)
	{
		uint_t<8192> vals[] =
		{
			uint_t<8192>::parse("1649538645321465538856413265789453212354678989785653212354878653424132356413"),
			uint_t<8192>::parse("6857395878795848394868797968574345687978925312545679532123453789786514908879"),
			uint_t<8192>::parse("4778531942543278913245238766235456868653123123456789899868652525252175004059"),
			uint_t<8192>::parse("6264134523486759565378935212111104040045605451123124323546474868757667862411"),
			uint_t<8192>::parse("1264272326427560643579865532334260603130262002462678317526421126423728617264"),
			uint_t<8192>::parse("4167276597145945283632459105940634353766783723761163234768566090999144264094"),
			uint_t<8192>::parse("4617627621942683227214159425626700616234303760377685879004613264267270421641"),
			uint_t<8192>::parse("1612674535351270645679076866462442652656767925402042569865421478786456782545"),
			uint_t<8192>::parse("1378654123789764321548066643503435773738979646034503430343045767542973182461"),
			uint_t<8192>::parse("4538978321534878351406064967605678378337866654650644303133226431626767879564"),
			uint_t<8192>::parse("2646134334576876864594921324653303750264246256728549413640672653423435767866"),
			uint_t<8192>::parse("4567659789654353432616261234350664572678024624362312312645273860795403736786"),
			uint_t<8192>::parse("1643573868945652331326456064672738264260264256425372867294025645789888376564"),
			uint_t<8192>::parse("4656738338865626133435757383064650762860724950414629645642567597260836656542"),
		};

		uint_t<8192> val2 = uint_t<8192>::parse("16579862451765276262462456641541236525214");
		decltype(f(vals[0], vals[0])) dd;

		t_start;
		for (std::size_t i = 0; i < count_3; ++i)
		{
			for (const auto &v : vals)
			{
				dd = f(v, val2); do_something(dd);
				dd = f(val2, v); do_something(dd);
				dd = f(v, v); do_something(dd);
				dd = f(val2, val2); do_something(dd);
			}
		}
		t_end("2) " << name << " unsigned 8192 (small)", 4 * count_3);
	}

	if constexpr (benchmark_do_int)
	{
		int_t<8192> vals[] =
		{
			int_t<8192>::parse("1649538645321465538856413265789453212354678989785653212354878653424132356413"),
			int_t<8192>::parse("6857395878795848394868797968574345687978925312545679532123453789786514908879"),
			int_t<8192>::parse("4778531942543278913245238766235456868653123123456789899868652525252175004059"),
			int_t<8192>::parse("6264134523486759565378935212111104040045605451123124323546474868757667862411"),
			int_t<8192>::parse("1264272326427560643579865532334260603130262002462678317526421126423728617264"),
			int_t<8192>::parse("4167276597145945283632459105940634353766783723761163234768566090999144264094"),
			int_t<8192>::parse("4617627621942683227214159425626700616234303760377685879004613264267270421641"),
			int_t<8192>::parse("1612674535351270645679076866462442652656767925402042569865421478786456782545"),
			int_t<8192>::parse("1378654123789764321548066643503435773738979646034503430343045767542973182461"),
			int_t<8192>::parse("4538978321534878351406064967605678378337866654650644303133226431626767879564"),
			int_t<8192>::parse("2646134334576876864594921324653303750264246256728549413640672653423435767866"),
			int_t<8192>::parse("4567659789654353432616261234350664572678024624362312312645273860795403736786"),
			int_t<8192>::parse("1643573868945652331326456064672738264260264256425372867294025645789888376564"),
			int_t<8192>::parse("4656738338865626133435757383064650762860724950414629645642567597260836656542"),
		};

		int_t<8192> val2 = int_t<8192>::parse("16579862451765276262462456641541236525214");
		decltype(f(vals[0], vals[0])) dd;

		t_start;
		for (std::size_t i = 0; i < count_3; ++i)
		{
			for (const auto &v : vals)
			{
				dd = f(v, val2); do_something(dd);
				dd = f(val2, v); do_something(dd);
				dd = f(v, v); do_something(dd);
				dd = f(val2, val2); do_something(dd);
			}
		}
		t_end("2) " << name << " signed 8192 (small)", 4 * count_3);
	}
	
	if constexpr (benchmark_do_bigint)
	{
		bigint vals[] =
		{
			bigint::parse("1649538645321465538856413265789453212354678989785653212354878653424132356413"),
			bigint::parse("6857395878795848394868797968574345687978925312545679532123453789786514908879"),
			bigint::parse("4778531942543278913245238766235456868653123123456789899868652525252175004059"),
			bigint::parse("6264134523486759565378935212111104040045605451123124323546474868757667862411"),
			bigint::parse("1264272326427560643579865532334260603130262002462678317526421126423728617264"),
			bigint::parse("4167276597145945283632459105940634353766783723761163234768566090999144264094"),
			bigint::parse("4617627621942683227214159425626700616234303760377685879004613264267270421641"),
			bigint::parse("1612674535351270645679076866462442652656767925402042569865421478786456782545"),
			bigint::parse("1378654123789764321548066643503435773738979646034503430343045767542973182461"),
			bigint::parse("4538978321534878351406064967605678378337866654650644303133226431626767879564"),
			bigint::parse("2646134334576876864594921324653303750264246256728549413640672653423435767866"),
			bigint::parse("4567659789654353432616261234350664572678024624362312312645273860795403736786"),
			bigint::parse("1643573868945652331326456064672738264260264256425372867294025645789888376564"),
			bigint::parse("4656738338865626133435757383064650762860724950414629645642567597260836656542"),
		};

		bigint val2 = bigint::parse("16579862451765276262462456641541236525214");
		decltype(f(vals[0], vals[0])) dd;

		t_start;
		for (std::size_t i = 0; i < count_4; ++i)
		{
			for (const auto &v : vals)
			{
				dd = f(v, val2); do_something(dd);
				dd = f(val2, v); do_something(dd);
				dd = f(v, v); do_something(dd);
				dd = f(val2, val2); do_something(dd);
			}
		}
		t_end("2) " << name << " bigint (small)", 4 * count_4);
	}

	std::cerr << '\n';
}

void benchmark_string(std::size_t count_1, std::size_t count_2, std::size_t count_3)
{
	const char *strs_1[] = {
		"164953864532146553885641326541323564132",
		"68573958787958483948687979685743988795",
		"47785319425432789132452387667540040594",
		"62641345234867594511231243576678624112",
		"12642723264275675264211264237286172641",
		"41672765971459452836324591059442640942",
		"46176276219426832272141594256267216412",
		"16126745353512798654214787864567825452",
		"13786541237897643215487675429731824613",
		"45389783215348783514264316267678795642",
		"26461343345768768645949213234357678665",
		"45676597896543534326162612343537367866",
		"16435738689456523313264567898883765642",
		"46567383388656261334357573838366565421",
	};
	const char *strs_2[] = {
		"164953864532146553323564132",
		"685739587879584839743988795",
		"477853194254327891540040594",
		"626413452348675945678624112",
		"126427232642756752286172641",
		"416727659714594528442640942",
		"461762762194268322267216412",
		"161267453535127986567825452",
		"137865412378976432731824613",
		"453897832153487835678795642",
		"264613433457687686357678665",
		"456765978965435343537367866",
		"164357386894565233883765642",
		"465673833886562613366565421",
	};
	const char *strs_3[] = {
		"10010111010101011101010101010101010101010011011001010101111110",
		"11010100100001111101010101010101001111010011100001111100010111",
		"01010100111010010101001101011010110001011010100101010111111010",
		"10100101010111010010100110101010101010011001010101100101011101",
		"10010100110101010110001101100010110101010111101100011110000110",
		"10101011010100101000101111111001000011110001110001110001110001",
		"10010100111001010101010100101000010101001000000000000000000000",
		"01010101010011111100000101011001110100110100101010101010101010",
		"00000000000000000000000000000000000000000000000000000000000000",
		"11111111111111111111111111111111111111111111111111111111111111",
		"10101010010101011010011010101010101010101010101010010101010101",
		"00110011001100110011001100110011001100110011001100110011001100",
		"11001100110011001100110011001100110011001100110011001100110011",
		"10011100011010101100101010101010001010100101010111010010101001",
	};

	{
		uint_t<128> dd;

		t_start;
		for (std::size_t i = 0; i < count_1; ++i)
		{
			for (const char *v : strs_1)
			{
				uint_t<128>::parse(dd, v, 10); do_something(dd);
			}
		}
		t_end("1) parse unsigned 128 dec", count_1);
	}
	{
		uint_t<128> dd;

		t_start;
		for (std::size_t i = 0; i < count_2; ++i)
		{
			for (const char *v : strs_2)
			{
				uint_t<128>::parse(dd, v, 16); do_something(dd);
			}
		}
		t_end("2) parse unsigned 128 hex", count_2);
	}
	{
		uint_t<128> dd;

		t_start;
		for (std::size_t i = 0; i < count_3; ++i)
		{
			for (const char *v : strs_3)
			{
				uint_t<128>::parse(dd, v, 2); do_something(dd);
			}
		}
		t_end("2) parse unsigned 128 bin", count_3);
	}
	{
		bigint dd;

		t_start;
		for (std::size_t i = 0; i < count_1; ++i)
		{
			for (const char *v : strs_1)
			{
				bigint::parse(dd, v, 10); do_something(dd);
			}
		}
		t_end("1) parse bigint dec", count_1);
	}
	{
		bigint dd;

		t_start;
		for (std::size_t i = 0; i < count_2; ++i)
		{
			for (const char *v : strs_2)
			{
				bigint::parse(dd, v, 16); do_something(dd);
			}
		}
		t_end("2) parse bigint hex", count_2);
	}
	{
		bigint dd;

		t_start;
		for (std::size_t i = 0; i < count_3; ++i)
		{
			for (const char *v : strs_3)
			{
				bigint::parse(dd, v, 2); do_something(dd);
			}
		}
		t_end("2) parse bigint bin", count_3);
	}

	std::cerr << '\n';
}

void detail_tests()
{
	assert(detail::_mul_u64(0 COMMA 0) == std::make_pair<std::uint64_t COMMA std::uint64_t>(0ull COMMA 0ull));
	assert(detail::_mul_u64(34346 COMMA 0) == std::make_pair<std::uint64_t COMMA std::uint64_t>(0ull COMMA 0ull));
	assert(detail::_mul_u64(0 COMMA 4564657) == std::make_pair<std::uint64_t COMMA std::uint64_t>(0ull COMMA 0ull));
	assert(detail::_mul_u64(3546482324346ull COMMA 45654513154657ull) == std::make_pair<std::uint64_t COMMA std::uint64_t>(3576792889872885050ull COMMA 8777317ull));
	assert(detail::_mul_u64(18446744073709551615ull COMMA 18446744073709551615ull) == std::make_pair<std::uint64_t COMMA std::uint64_t>(1ull COMMA 18446744073709551614ull));

	assert(detail::highest_set_bit(0) == 0);
	for (std::size_t i = 0; i < 64; ++i)
	{
		assert(detail::highest_set_bit((std::uint64_t)1 << i) == i);
	}
	assert(detail::highest_set_bit(uint_t<512>(0)) == 0);
}
void smol_tests()
{
	assert(sstostr(uint_t<128>(2) * uint_t<128>(13)) == "26");
	assert(sstostr(uint_t<128>(2) + uint_t<128>(13)) == "15");
	assert(sstostr(uint_t<128>(20) - uint_t<128>(13)) == "7");
	assert(sstostr(uint_t<128>(2) / uint_t<128>(13)) == "0");
	assert(sstostr(uint_t<128>(2) % uint_t<128>(13)) == "2");

	assert(sstostr(2u * uint_t<128>(13)) == "26");
	assert(sstostr(2u + uint_t<128>(13)) == "15");
	assert(sstostr(20u - uint_t<128>(13)) == "7");
	assert(sstostr(uint_t<128>(2) / 13u) == "0");
	assert(sstostr(uint_t<128>(2) % 13u) == "2");

	assert(sstostr(uint_t<128>(52) * uint_t<128>(3)) == "156");
	assert(sstostr(uint_t<128>(52) + uint_t<128>(3)) == "55");
	assert(sstostr(uint_t<128>(52) / uint_t<128>(3)) == "17");
	assert(sstostr(uint_t<128>(52) % uint_t<128>(3)) == "1");

	assert(sstostr(int_t<128>(2) * int_t<128>(-13)) == "-26");
	assert(sstostr(int_t<128>(2) + int_t<128>(-13)) == "-11");
	assert(sstostr(int_t<128>(2) / int_t<128>(-13)) == "0");
	assert(sstostr(int_t<128>(2) % int_t<128>(13)) == "2");

	assert(sstostr(int_t<128>(2) & int_t<128>(-13)) == "2");
	assert(sstostr(int_t<128>(2) | int_t<128>(-13)) == "-13");
	assert(sstostr(int_t<128>(2) ^ int_t<128>(-13)) == "-15");

	assert(sstostr(int_t<128>(-13) << 1) == "-26");
	assert(sstostr(int_t<128>(-13) << 2) == "-52");

	assert(sstostr(int_t<128>(-13) >> 1) == "-7");
	assert(sstostr(int_t<128>(-13) >> 2) == "-4");
}
void bigboi_tests()
{
	uint_t<8192> big_mul = 947563412;
	assert(sstostr(big_mul) == "947563412");

	big_mul *= big_mul;
	assert(sstostr(big_mul) == "897876419761081744");

	big_mul *= big_mul;
	assert(sstostr(big_mul) == "806182065162978263317234893050081536");

	big_mul *= big_mul;
	assert(sstostr(big_mul) == "649929522190444530768966266652239714986927778859800606430979456248119296");

	big_mul += 17;
	big_mul += (decltype(big_mul))0;
	big_mul = big_mul + (decltype(big_mul))0;
	big_mul = (decltype(big_mul))0 + big_mul;
	big_mul = big_mul + 0;
	big_mul = 0 + big_mul;
	assert(sstostr(big_mul) == "649929522190444530768966266652239714986927778859800606430979456248119313");
	big_mul -= 0;
	big_mul = big_mul - (decltype(big_mul))0;
	big_mul = (decltype(big_mul))0 - big_mul;
	big_mul = big_mul - 0;
	big_mul = 0 - big_mul;
	assert(sstostr(big_mul) == "649929522190444530768966266652239714986927778859800606430979456248119313");
	big_mul &= (decltype(big_mul))-1;
	big_mul &= -1;
	big_mul = big_mul & (decltype(big_mul))-1;
	big_mul = big_mul & -1;
	big_mul = (decltype(big_mul))-1 & big_mul;
	big_mul = -1 & big_mul;
	assert(sstostr(big_mul) == "649929522190444530768966266652239714986927778859800606430979456248119313");
	big_mul |= (decltype(big_mul))0;
	big_mul |= 0;
	big_mul = big_mul | (decltype(big_mul))0;
	big_mul = big_mul | 0;
	big_mul = (decltype(big_mul))0 | big_mul;
	big_mul = 0 | big_mul;
	assert(sstostr(big_mul) == "649929522190444530768966266652239714986927778859800606430979456248119313");
	big_mul ^= (decltype(big_mul))0;
	big_mul ^= 0;
	big_mul = big_mul ^ (decltype(big_mul))0;
	big_mul = big_mul ^ 0;
	big_mul = (decltype(big_mul))0 ^ big_mul;
	big_mul = 0 ^ big_mul;
	assert(sstostr(big_mul) == "649929522190444530768966266652239714986927778859800606430979456248119313");

	big_mul <<= 621;
	assert(sstostr(big_mul) == "56557938587827109863786197009098171333880971620065133028223237982646344014826583116417132"
		"577306321974342507460953540245890763571403111747364084775326663465806403623083715574"
		"23307676007538288423649789830139543353463939819407429243267044228077146543261877272576");

	big_mul &= big_mul;
	assert(sstostr(big_mul) == "56557938587827109863786197009098171333880971620065133028223237982646344014826583116417132"
		"577306321974342507460953540245890763571403111747364084775326663465806403623083715574"
		"23307676007538288423649789830139543353463939819407429243267044228077146543261877272576");

	big_mul |= big_mul;
	assert(sstostr(big_mul) == "56557938587827109863786197009098171333880971620065133028223237982646344014826583116417132"
		"577306321974342507460953540245890763571403111747364084775326663465806403623083715574"
		"23307676007538288423649789830139543353463939819407429243267044228077146543261877272576");

	auto post_dec = big_mul--;
	assert(sstostr(post_dec) == "56557938587827109863786197009098171333880971620065133028223237982646344014826583116417132"
		"577306321974342507460953540245890763571403111747364084775326663465806403623083715574"
		"23307676007538288423649789830139543353463939819407429243267044228077146543261877272576");
	assert(sstostr(big_mul) == "56557938587827109863786197009098171333880971620065133028223237982646344014826583116417132"
		"577306321974342507460953540245890763571403111747364084775326663465806403623083715574"
		"23307676007538288423649789830139543353463939819407429243267044228077146543261877272575");

	auto post_inc = big_mul++;
	assert(sstostr(post_inc) == "56557938587827109863786197009098171333880971620065133028223237982646344014826583116417132"
		"577306321974342507460953540245890763571403111747364084775326663465806403623083715574"
		"23307676007538288423649789830139543353463939819407429243267044228077146543261877272575");
	assert(sstostr(big_mul) == "56557938587827109863786197009098171333880971620065133028223237982646344014826583116417132"
		"577306321974342507460953540245890763571403111747364084775326663465806403623083715574"
		"23307676007538288423649789830139543353463939819407429243267044228077146543261877272576");

	auto pre_inc = ++big_mul;
	assert(sstostr(pre_inc) == "56557938587827109863786197009098171333880971620065133028223237982646344014826583116417132"
		"577306321974342507460953540245890763571403111747364084775326663465806403623083715574"
		"23307676007538288423649789830139543353463939819407429243267044228077146543261877272577");
	assert(sstostr(big_mul) == "56557938587827109863786197009098171333880971620065133028223237982646344014826583116417132"
		"577306321974342507460953540245890763571403111747364084775326663465806403623083715574"
		"23307676007538288423649789830139543353463939819407429243267044228077146543261877272577");
	assert(pre_inc == big_mul);
	assert(pre_inc <= big_mul);
	assert(pre_inc >= big_mul);

	auto pre_dec = --big_mul;
	assert(sstostr(pre_dec) == "56557938587827109863786197009098171333880971620065133028223237982646344014826583116417132"
		"577306321974342507460953540245890763571403111747364084775326663465806403623083715574"
		"23307676007538288423649789830139543353463939819407429243267044228077146543261877272576");
	assert(sstostr(big_mul) == "56557938587827109863786197009098171333880971620065133028223237982646344014826583116417132"
		"577306321974342507460953540245890763571403111747364084775326663465806403623083715574"
		"23307676007538288423649789830139543353463939819407429243267044228077146543261877272576");
	assert(pre_dec == big_mul);
	assert(pre_dec <= big_mul);
	assert(pre_dec >= big_mul);

	big_mul >>= 191;
	assert(sstostr(big_mul) == "1802039889491867190810467664715451904918476332004551151301726717494487639814495660570778575783656"
		"106649698951732305310318830222989557694740573768908873663803263415298739647066704572678996885278044454912");
	big_mul <<= 73;
	assert(sstostr(big_mul) == "1701978554986102599652944396520696900597046164730926970313320209887640881966509866088480298291994"
		"7794655263581165454491067120838946815593625140970325732225722954705464939198103933056536439508186725987781864851905868153749504");

	big_mul += big_mul;
	assert(sstostr(big_mul) == "34039571099722051993058887930413938011940923294618539406266404197752817639330197321769605965839895589"
		"310527162330908982134241677893631187250281940651464451445909410929878396207866113072879016373451975563729703811736307499008");

	big_mul <<= 1;
	assert(sstostr(big_mul) == "68079142199444103986117775860827876023881846589237078812532808395505635278660394643539211931679791178"
		"621054324661817964268483355787262374500563881302928902891818821859756792415732226145758032746903951127459407623472614998016");

	big_mul = 465;
	assert(sstostr(big_mul) == "465");

	assert(sstostr(big_mul + 75u) == "540");
	assert(sstostr(big_mul ^ big_mul) == "0");
}
void numeric_limits_tests()
{
	static_assert(std::is_reference_v<decltype(std::numeric_limits<uint_t<512>>::min())>, "large constexpr object returned by value");
	static_assert(std::is_reference_v<decltype(std::numeric_limits<int_t<512>>::max())>, "large constexpr object returned by value");
	
	static_assert(std::is_reference_v<decltype(uint_t<512>::min())>, "large constexpr object returned by value");
	static_assert(std::is_reference_v<decltype(int_t<512>::max())>, "large constexpr object returned by value");

	tostr_fmt dec;

	assert(dec(std::numeric_limits<uint_t<512>>::min()) == "0");
	assert(dec(std::numeric_limits<uint_t<512>>::max()) == "13407807929942597099574024998205846127479365820592393377723561443721764030073546976801874298166903427690031858186486050853753882811946569946433649006084095");
	assert(dec(std::numeric_limits<uint_t<512>>::lowest()) == "0");

	assert(dec(std::numeric_limits<int_t<512>>::min()) == "-6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042048");
	assert(dec(std::numeric_limits<int_t<512>>::max()) == "6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042047");
	assert(dec(std::numeric_limits<int_t<512>>::lowest()) == "-6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042048");

	assert(std::numeric_limits<uint_t<128>>::min().tostr(dec) == "0");
	assert(std::numeric_limits<uint_t<128>>::max().tostr(dec) == "340282366920938463463374607431768211455");
	assert(std::numeric_limits<uint_t<128>>::lowest().tostr(dec) == "0");

	assert(std::numeric_limits<int_t<128>>::min().tostr(dec) == "-170141183460469231731687303715884105728");
	assert(std::numeric_limits<int_t<128>>::max().tostr(dec) == "170141183460469231731687303715884105727");
	assert(std::numeric_limits<int_t<128>>::lowest().tostr(dec) == "-170141183460469231731687303715884105728");

	assert(dec(std::numeric_limits<uint_t<64>>::min()) == "0");
	assert(dec(std::numeric_limits<uint_t<64>>::max()) == "18446744073709551615");
	assert(dec(std::numeric_limits<uint_t<64>>::lowest()) == "0");

	assert(dec(std::numeric_limits<int_t<64>>::min()) == "-9223372036854775808");
	assert(dec(std::numeric_limits<int_t<64>>::max()) == "9223372036854775807");
	assert(dec(std::numeric_limits<int_t<64>>::lowest()) == "-9223372036854775808");

	assert(&uint_t<512>::min() == &std::numeric_limits<uint_t<512>>::min()); // we ensure reference equality to save space
	assert(&int_t<512>::max() == &std::numeric_limits<int_t<512>>::max());
}
void simple_parse_tests()
{
	assert(sstostr(uint_t<256>::parse("184738578493786576848368767654873647697975746763857664")) == "184738578493786576848368767654873647697975746763857664");
	assert(sstostr(uint_t<256>::parse("18")) == "18");
	assert(uint_t<256>::parse("18") == 18u);

	assert(sstostr(int_t<256>::parse("184738578493786576848368767654873647697975746763857664")) == "184738578493786576848368767654873647697975746763857664");
	assert(sstostr(int_t<256>::parse("18")) == "18");
	assert(int_t<256>::parse("18") == 18);
	assert(int_t<256>::parse("+18") == 18);

	assert(sstostr(int_t<256>::parse("-184738578493786576848368767654873647697975746763857664")) == "-184738578493786576848368767654873647697975746763857664");
	assert(sstostr(int_t<256>::parse("-18")) == "-18");
	assert(int_t<256>::parse("-18") == -18);

	assert(sstostr(int_t<256>::parse(" -18")) == "-18");
	assert(sstostr(int_t<256>::parse("-18 ")) == "-18");
	assert(sstostr(int_t<256>::parse(" -18 ")) == "-18");
	assert(int_t<256>::parse("-18") == -18);
	assert(int_t<256>::parse("-18 ") == -18);
	assert(int_t<256>::parse(" -18") == -18);

	uint_t<512> temp_u512;
	int_t<512> temp_i512;
	uint_t<256> temp_u256;
	int_t<256> temp_i256;

	assert_throws(std::invalid_argument, uint_t<512>::parse(temp_u512, "- 18"));
	assert_throws(std::invalid_argument, uint_t<512>::parse(temp_u512, "+ 18"));

	assert_throws(std::invalid_argument, uint_t<512>::parse(temp_u512, "abc"));
	assert_throws(std::invalid_argument, uint_t<512>::parse(temp_u512, " abc"));
	assert_throws(std::invalid_argument, uint_t<512>::parse(temp_u512, "abc "));
	assert_throws(std::invalid_argument, uint_t<512>::parse(temp_u512, " abc "));

	assert_throws(std::invalid_argument, uint_t<512>::parse(temp_u512, "a18"));
	assert_throws(std::invalid_argument, uint_t<512>::parse(temp_u512, "18a"));

	assert_throws(std::invalid_argument, uint_t<512>::parse(temp_u512, "a 18"));
	assert_throws(std::invalid_argument, uint_t<512>::parse(temp_u512, "18 a"));

	// ----------

	assert_throws(std::invalid_argument, int_t<512>::parse(temp_i512, "- 18"));
	assert_throws(std::invalid_argument, int_t<512>::parse(temp_i512, "+ 18"));

	assert_throws(std::invalid_argument, int_t<512>::parse(temp_i512, "abc"));
	assert_throws(std::invalid_argument, int_t<512>::parse(temp_i512, " abc"));
	assert_throws(std::invalid_argument, int_t<512>::parse(temp_i512, "abc "));
	assert_throws(std::invalid_argument, int_t<512>::parse(temp_i512, " abc "));

	assert_throws(std::invalid_argument, int_t<512>::parse(temp_i512, "a18"));
	assert_throws(std::invalid_argument, int_t<512>::parse(temp_i512, "18a"));

	assert_throws(std::invalid_argument, int_t<512>::parse(temp_i512, "a 18"));
	assert_throws(std::invalid_argument, int_t<512>::parse(temp_i512, "18 a"));

	// base parsing tests

	assert_throws(std::invalid_argument, uint_t<256>::parse(temp_u256, "0x764"));
	assert_throws(std::invalid_argument, uint_t<256>::parse(temp_u256, "0x764", 10));
	assert_throws(std::invalid_argument, uint_t<256>::parse(temp_u256, "0x764", 8));
	assert_throws(std::invalid_argument, uint_t<256>::parse(temp_u256, "0x764", 16));
	assert(sstostr(uint_t<256>::parse("0x764", 0)) == "1892");
	assert(sstostr(uint_t<256>::parse("0o764", 0)) == "500");
	assert(sstostr(uint_t<256>::parse("764", 0)) == "764");
	assert(sstostr(uint_t<256>::parse("0764", 0)) == "764");

	assert(sstostr(uint_t<256>::parse("167345562314657356422211643535261643535621", 10)) == "167345562314657356422211643535261643535621");
	assert(sstostr(uint_t<256>::parse("167345562314657356422211643535261643535621", 8)) == "19846815955032434075951735996930767761");
	assert(sstostr(uint_t<256>::parse("167345562314657356422211643535261643535621", 16)) == "32811116214936888653588305588611479218962394469921");

	assert(sstostr(int_t<256>::parse("167345562314657356422211643535261643535621", 10)) == "167345562314657356422211643535261643535621");
	assert(sstostr(int_t<256>::parse("+167345562314657356422211643535261643535621", 10)) == "167345562314657356422211643535261643535621");
	assert(sstostr(int_t<256>::parse("-167345562314657356422211643535261643535621", 10)) == "-167345562314657356422211643535261643535621");

	assert(sstostr(int_t<256>::parse("167345562314657356422211643535261643535621", 8)) == "19846815955032434075951735996930767761");
	assert(sstostr(int_t<256>::parse("167345562314657356422211643535261643535621", 16)) == "32811116214936888653588305588611479218962394469921");

	// zero parsing edge cases

	auto zero_tester = [](auto v) {
		std::conditional_t<std::is_signed_v<decltype(v)>, int, unsigned> zero = 0;

		assert(decltype(v)::parse("0") == zero);
		assert(decltype(v)::parse(" 0") == zero);
		assert(decltype(v)::parse("0 ") == zero);
		assert(decltype(v)::parse(" 0 ") == zero);

		assert(decltype(v)::parse("0", 10) == zero);
		assert(decltype(v)::parse(" 0", 10) == zero);
		assert(decltype(v)::parse("0 ", 10) == zero);
		assert(decltype(v)::parse(" 0 ", 10) == zero);

		assert(decltype(v)::parse("0", 8) == zero);
		assert(decltype(v)::parse(" 0", 8) == zero);
		assert(decltype(v)::parse("0 ", 8) == zero);
		assert(decltype(v)::parse(" 0 ", 8) == zero);

		assert(decltype(v)::parse("0", 16) == zero);
		assert(decltype(v)::parse(" 0", 16) == zero);
		assert(decltype(v)::parse("0 ", 16) == zero);
		assert(decltype(v)::parse(" 0 ", 16) == zero);

		assert(decltype(v)::parse("0", 0) == zero);
		assert(decltype(v)::parse(" 0", 0) == zero);
		assert(decltype(v)::parse("0 ", 0) == zero);
		assert(decltype(v)::parse(" 0 ", 0) == zero);

		assert(!decltype(v)::try_parse(v, "0x", 0));
		assert(!decltype(v)::try_parse(v, " 0x", 0));
		assert(!decltype(v)::try_parse(v, "0x ", 0));
		assert(!decltype(v)::try_parse(v, " 0x ", 0));

		assert(!decltype(v)::try_parse(v, "0 x", 0));
		assert(!decltype(v)::try_parse(v, " 0 x", 0));
		assert(!decltype(v)::try_parse(v, "0 x ", 0));
		assert(!decltype(v)::try_parse(v, " 0 x ", 0));

		assert(!decltype(v)::try_parse(v, "0 1", 0));
		assert(!decltype(v)::try_parse(v, " 0 1", 0));
		assert(!decltype(v)::try_parse(v, "0 1 ", 0));
		assert(!decltype(v)::try_parse(v, " 0 1 ", 0));
	};

	zero_tester(uint_t<256>{});
	zero_tester(int_t<256>{});
	zero_tester(bigint{});

	{
		std::istringstream ss("  543   ");
		uint_t<256> v;
		ss >> v;
		assert((bool)ss);
		assert(v == 543u);
		assert(ss.get() == ' ');
	}
	{
		std::istringstream ss("  543   ");
		int_t<256> v;
		ss >> v;
		assert((bool)ss);
		assert(v == 543);
		assert(ss.get() == ' ');
	}
	{
		std::istringstream ss("  543   ");
		bigint v;
		ss >> v;
		assert((bool)ss);
		assert(v == 543);
		assert(ss.get() == ' ');
	}
}
void big_parse_tests()
{
	{
		uint_t<256> temp;
		bigint big_temp;

		assert(!uint_t<256>::try_parse(temp, "0x764"));
		assert(!uint_t<256>::try_parse(temp, "0x764", 10));
		assert(!uint_t<256>::try_parse(temp, "0x764", 8));
		assert(!uint_t<256>::try_parse(temp, "0x764", 16));
		assert(uint_t<256>::try_parse(temp, "0x764", 0));

		assert(!bigint::try_parse(big_temp, "0x764"));
		assert(!bigint::try_parse(big_temp, "0x764", 10));
		assert(!bigint::try_parse(big_temp, "0x764", 8));
		assert(!bigint::try_parse(big_temp, "0x764", 16));
		assert(bigint::try_parse(big_temp, "0x764", 0));
	}

	{
		int_t<256> temp;
		assert(int_t<256>::try_parse(temp, "-167345562314657356422211643535261643535621", 10));
		assert(sstostr(temp) == "-167345562314657356422211643535261643535621");

		assert(int_t<256>::try_parse(temp, "167345562314657356422211643535261643535621", 8));
		assert(sstostr(temp) == "19846815955032434075951735996930767761");

		assert(int_t<256>::try_parse(temp, "167345562314657356422211643535261643535621", 16));
		assert(sstostr(temp) == "32811116214936888653588305588611479218962394469921");
	}

	tostr_fmt tobin(2);

	tostr_fmt tobin_showbase = tobin;
	tobin_showbase.showbase = true;

	tostr_fmt tobin_uppercase = tobin;
	tobin_uppercase.uppercase = true;

	parse_fmt frombin(2);
	parse_fmt fromunk(0);

	{
		bigint v = 453;
		assert(sstostr(v) == "453");
		assert(sstostr(std::setw(10), v) == "       453");
		assert(sstostr(std::showpos, v) == "+453");
		assert(sstostr(std::showpos, std::setw(10), v) == "      +453");
		assert(sstostr(std::dec, v) == "453");
		assert(sstostr(std::oct, v) == "0705");
		assert(sstostr(std::hex, v) == "1c5");
		assert(sstostr(std::dec, std::showpos, v) == "+453");
		assert(sstostr(std::oct, std::showpos, v) == "0705");
		assert(sstostr(std::hex, std::showpos, v) == "1c5");
		assert(sstostr(std::dec, std::showbase, v) == "453");
		assert(sstostr(std::oct, std::showbase, v) == "0o0705");
		assert(sstostr(std::hex, std::showbase, v) == "0x1c5");
		assert(sstostr(std::dec, std::showbase, std::showpos, v) == "+453");
		assert(sstostr(std::oct, std::showbase, std::showpos, v) == "0o0705");
		assert(sstostr(std::hex, std::showbase, std::showpos, v) == "0x1c5");
		assert(sstostr(std::dec, std::showbase, std::showpos, std::setw(10), v) == "      +453");
		assert(sstostr(std::oct, std::showbase, std::showpos, std::setw(10), v) == "    0o0705");
		assert(sstostr(std::hex, std::showbase, std::showpos, std::setw(10), v) == "     0x1c5");

		assert(tobin(v) == "0111000101");
		assert(tobin_showbase(v) == "0b0111000101");
		assert(tobin_uppercase(v) == "0111000101");

		assert(frombin("0111000101") == 453);
		assert(fromunk("0b0111000101") == 453);

		assert_throws(std::invalid_argument, (void)frombin("01110001012"));
		assert_throws(std::invalid_argument, (void)frombin("0b01110001012"));
		assert_throws(std::invalid_argument, (void)frombin("0b 0111000101"));

		assert(bigint::parse("453") == 453);
		assert(bigint::parse("       453") == 453);
		assert(bigint::parse("+453") == 453);
		assert(bigint::parse("      +453") == 453);
		assert(bigint::parse("453", 10) == 453);
		assert(bigint::parse("0705", 8) == 453);
		assert(bigint::parse("1c5", 16) == 453);
		assert(bigint::parse("+453", 10) == 453);
		assert(bigint::parse("0705", 8) == 453);
		assert(bigint::parse("1c5", 16) == 453);
		assert(bigint::parse("453", 0) = 453);
		assert(bigint::parse("0o0705", 0) == 453);
		assert(bigint::parse("0x1c5", 0) == 453);
		assert(bigint::parse("+453", 0) == 453);
		assert(bigint::parse("0o0705", 0) == 453);
		assert(bigint::parse("0x1c5", 0) == 453);
		assert(bigint::parse("      +453", 0) == 453);
		assert(bigint::parse("    0o0705", 0) == 453);
		assert(bigint::parse("     0x1c5", 0) == 453);
	}
	{
		bigint v = 1;
		assert(sstostr(v) == "1");
		assert(sstostr(std::setw(10), v) == "         1");
		assert(sstostr(std::showpos, v) == "+1");
		assert(sstostr(std::showpos, std::setw(10), v) == "        +1");
		assert(sstostr(std::dec, v) == "1");
		assert(sstostr(std::oct, v) == "1");
		assert(sstostr(std::hex, v) == "1");
		assert(sstostr(std::dec, std::showpos, v) == "+1");
		assert(sstostr(std::oct, std::showpos, v) == "1");
		assert(sstostr(std::hex, std::showpos, v) == "1");
		assert(sstostr(std::dec, std::showbase, v) == "1");
		assert(sstostr(std::oct, std::showbase, v) == "0o1");
		assert(sstostr(std::hex, std::showbase, v) == "0x1");
		assert(sstostr(std::dec, std::showbase, std::showpos, v) == "+1");
		assert(sstostr(std::oct, std::showbase, std::showpos, v) == "0o1");
		assert(sstostr(std::hex, std::showbase, std::showpos, v) == "0x1");
		assert(sstostr(std::dec, std::showbase, std::showpos, std::setw(10), v) == "        +1");
		assert(sstostr(std::oct, std::showbase, std::showpos, std::setw(10), v) == "       0o1");
		assert(sstostr(std::hex, std::showbase, std::showpos, std::setw(10), v) == "       0x1");

		assert(tobin(v) == "01");
		assert(tobin_showbase(v) == "0b01");
		assert(tobin_uppercase(v) == "01");

		assert(frombin("01") == 1);
		assert(fromunk("0b01") == 1);

		assert(bigint::parse("1") == 1);
		assert(bigint::parse("         1") == 1);
		assert(bigint::parse("+1") == 1);
		assert(bigint::parse("        +1") == 1);
		assert(bigint::parse("1", 10) == 1);
		assert(bigint::parse("1", 8) == 1);
		assert(bigint::parse("1", 16) == 1);
		assert(bigint::parse("+1", 10) == 1);
		assert(bigint::parse("1", 8) == 1);
		assert(bigint::parse("1", 16) == 1);
		assert(bigint::parse("1", 0) == 1);
		assert(bigint::parse("0o1", 0) == 1);
		assert(bigint::parse("0x1", 0) == 1);
		assert(bigint::parse("+1", 0) == 1);
		assert(bigint::parse("0o1", 0) == 1);
		assert(bigint::parse("0x1", 0) == 1);
		assert(bigint::parse("        +1", 0) == 1);
		assert(bigint::parse("       0o1", 0) == 1);
		assert(bigint::parse("       0x1", 0) == 1);
	}
	{
		bigint v = 4;
		assert(sstostr(v) == "4");
		assert(sstostr(std::setw(10), v) == "         4");
		assert(sstostr(std::showpos, v) == "+4");
		assert(sstostr(std::showpos, std::setw(10), v) == "        +4");
		assert(sstostr(std::dec, v) == "4");
		assert(sstostr(std::oct, v) == "04");
		assert(sstostr(std::hex, v) == "4");
		assert(sstostr(std::dec, std::showpos, v) == "+4");
		assert(sstostr(std::oct, std::showpos, v) == "04");
		assert(sstostr(std::hex, std::showpos, v) == "4");
		assert(sstostr(std::dec, std::showbase, v) == "4");
		assert(sstostr(std::oct, std::showbase, v) == "0o04");
		assert(sstostr(std::hex, std::showbase, v) == "0x4");
		assert(sstostr(std::dec, std::showbase, std::showpos, v) == "+4");
		assert(sstostr(std::oct, std::showbase, std::showpos, v) == "0o04");
		assert(sstostr(std::hex, std::showbase, std::showpos, v) == "0x4");
		assert(sstostr(std::dec, std::showbase, std::showpos, std::setw(10), v) == "        +4");
		assert(sstostr(std::oct, std::showbase, std::showpos, std::setw(10), v) == "      0o04");
		assert(sstostr(std::hex, std::showbase, std::showpos, std::setw(10), v) == "       0x4");

		assert(tobin(v) == "0100");
		assert(tobin_showbase(v) == "0b0100");
		assert(tobin_uppercase(v) == "0100");

		assert(frombin("0100") == 4);
		assert(fromunk("0b0100") == 4);

		assert(bigint::parse("4") == 4);
		assert(bigint::parse("         4") == 4);
		assert(bigint::parse("+4") = 4);
		assert(bigint::parse("        +4") == 4);
		assert(bigint::parse("4", 10) == 4);
		assert(bigint::parse("04", 8) == 4);
		assert(bigint::parse("4", 16) == 4);
		assert(bigint::parse("+4", 10) == 4);
		assert(bigint::parse("04", 8) == 4);
		assert(bigint::parse("4", 16) == 4);
		assert(bigint::parse("4", 0) == 4);
		assert(bigint::parse("0o04", 0) == 4);
		assert(bigint::parse("0x4", 0) == 4);
		assert(bigint::parse("+4", 0) == 4);
		assert(bigint::parse("0o04", 0) == 4);
		assert(bigint::parse("0x4", 0) == 4);
		assert(bigint::parse("        +4", 0) == 4);
		assert(bigint::parse("      0o04", 0) == 4);
		assert(bigint::parse("       0x4", 0) == 4);
	}
	{
		bigint v = 8;
		assert(sstostr(v) == "8");
		assert(sstostr(std::setw(10), v) == "         8");
		assert(sstostr(std::showpos, v) == "+8");
		assert(sstostr(std::showpos, std::setw(10), v) == "        +8");
		assert(sstostr(std::dec, v) == "8");
		assert(sstostr(std::oct, v) == "10");
		assert(sstostr(std::hex, v) == "08");
		assert(sstostr(std::dec, std::showpos, v) == "+8");
		assert(sstostr(std::oct, std::showpos, v) == "10");
		assert(sstostr(std::hex, std::showpos, v) == "08");
		assert(sstostr(std::dec, std::showbase, v) == "8");
		assert(sstostr(std::oct, std::showbase, v) == "0o10");
		assert(sstostr(std::hex, std::showbase, v) == "0x08");
		assert(sstostr(std::dec, std::showbase, std::showpos, v) == "+8");
		assert(sstostr(std::oct, std::showbase, std::showpos, v) == "0o10");
		assert(sstostr(std::hex, std::showbase, std::showpos, v) == "0x08");
		assert(sstostr(std::dec, std::showbase, std::showpos, std::setw(10), v) == "        +8");
		assert(sstostr(std::oct, std::showbase, std::showpos, std::setw(10), v) == "      0o10");
		assert(sstostr(std::hex, std::showbase, std::showpos, std::setw(10), v) == "      0x08");

		assert(tobin(v) == "01000");
		assert(tobin_showbase(v) == "0b01000");
		assert(tobin_uppercase(v) == "01000");

		assert(frombin("01000") == 8);
		assert(fromunk("0b01000") == 8);

		assert(bigint::parse("8") == 8);
		assert(bigint::parse("         8") == 8);
		assert(bigint::parse("+8") == 8);
		assert(bigint::parse("        +8") == 8);
		assert(bigint::parse("8", 10) == 8);
		assert(bigint::parse("10", 8) == 8);
		assert(bigint::parse("08", 16) == 8);
		assert(bigint::parse("+8", 10) == 8);
		assert(bigint::parse("10", 8) == 8);
		assert(bigint::parse("08", 16) == 8);
		assert(bigint::parse("8", 0) == 8);
		assert(bigint::parse("0o10", 0) == 8);
		assert(bigint::parse("0x08", 0) == 8);
		assert(bigint::parse("+8", 0) == 8);
		assert(bigint::parse("0o10", 0) == 8);
		assert(bigint::parse("0x08", 0) == 8);
		assert(bigint::parse("        +8", 0) == 8);
		assert(bigint::parse("      0o10", 0) == 8);
		assert(bigint::parse("      0x08", 0) == 8);
	}
	{
		bigint v = 0x966;
		assert(sstostr(v) == "2406");
		assert(sstostr(std::setw(10), v) == "      2406");
		assert(sstostr(std::showpos, v) == "+2406");
		assert(sstostr(std::showpos, std::setw(10), v) == "     +2406");
		assert(sstostr(std::dec, v) == "2406");
		assert(sstostr(std::oct, v) == "04546");
		assert(sstostr(std::hex, v) == "0966");
		assert(sstostr(std::dec, std::showpos, v) == "+2406");
		assert(sstostr(std::oct, std::showpos, v) == "04546");
		assert(sstostr(std::hex, std::showpos, v) == "0966");
		assert(sstostr(std::dec, std::showbase, v) == "2406");
		assert(sstostr(std::oct, std::showbase, v) == "0o04546");
		assert(sstostr(std::hex, std::showbase, v) == "0x0966");
		assert(sstostr(std::dec, std::showbase, std::showpos, v) == "+2406");
		assert(sstostr(std::oct, std::showbase, std::showpos, v) == "0o04546");
		assert(sstostr(std::hex, std::showbase, std::showpos, v) == "0x0966");
		assert(sstostr(std::dec, std::showbase, std::showpos, std::setw(10), v) == "     +2406");
		assert(sstostr(std::oct, std::showbase, std::showpos, std::setw(10), v) == "   0o04546");
		assert(sstostr(std::hex, std::showbase, std::showpos, std::setw(10), v) == "    0x0966");

		assert(tobin(v) == "0100101100110");
		assert(tobin_showbase(v) == "0b0100101100110");
		assert(tobin_uppercase(v) == "0100101100110");

		assert(frombin("0100101100110") == 2406);
		assert(fromunk("0b0100101100110") == 2406);

		assert(bigint::parse("2406") == 0x966);
		assert(bigint::parse("      2406") == 0x966);
		assert(bigint::parse("+2406") == 0x966);
		assert(bigint::parse("     +2406") == 0x966);
		assert(bigint::parse("2406", 10) == 0x966);
		assert(bigint::parse("04546", 8) == 0x966);
		assert(bigint::parse("0966", 16) == 0x966);
		assert(bigint::parse("+2406", 10) == 0x966);
		assert(bigint::parse("04546", 8) == 0x966);
		assert(bigint::parse("0966", 16) == 0x966);
		assert(bigint::parse("2406", 0) == 0x966);
		assert(bigint::parse("0o04546", 0) == 0x966);
		assert(bigint::parse("0x0966", 0) == 0x966);
		assert(bigint::parse("+2406", 0) == 0x966);
		assert(bigint::parse("0o04546", 0) == 0x966);
		assert(bigint::parse("0x0966", 0) == 0x966);
		assert(bigint::parse("     +2406", 0) == 0x966);
		assert(bigint::parse("   0o04546", 0) == 0x966);
		assert(bigint::parse("    0x0966", 0) == 0x966);
	}
	{
		bigint v = -1;
		assert(sstostr(v) == "-1");
		assert(sstostr(std::setw(10), v) == "        -1");
		assert(sstostr(std::showpos, v) == "-1");
		assert(sstostr(std::showpos, std::setw(10), v) == "        -1");
		assert(sstostr(std::dec, v) == "-1");
		assert(sstostr(std::oct, v) == "7");
		assert(sstostr(std::hex, v) == "f");
		assert(sstostr(std::dec, std::showpos, v) == "-1");
		assert(sstostr(std::oct, std::showpos, v) == "7");
		assert(sstostr(std::hex, std::showpos, v) == "f");
		assert(sstostr(std::dec, std::showbase, v) == "-1");
		assert(sstostr(std::oct, std::showbase, v) == "0o7");
		assert(sstostr(std::hex, std::showbase, v) == "0xf");
		assert(sstostr(std::dec, std::showbase, std::showpos, v) == "-1");
		assert(sstostr(std::oct, std::showbase, std::showpos, v) == "0o7");
		assert(sstostr(std::hex, std::showbase, std::showpos, v) == "0xf");
		assert(sstostr(std::dec, std::showbase, std::showpos, std::setw(10), v) == "        -1");
		assert(sstostr(std::oct, std::showbase, std::showpos, std::setw(10), v) == "       0o7");
		assert(sstostr(std::hex, std::showbase, std::showpos, std::setw(10), v) == "       0xf");

		assert(tobin(v) == "1");
		assert(tobin_showbase(v) == "0b1");
		assert(tobin_uppercase(v) == "1");

		assert(frombin("1") == -1);
		assert(fromunk("0b1") == -1);

		assert(bigint::parse("-1") == -1);
		assert(bigint::parse("        -1") == -1);
		assert(bigint::parse("-1") == -1);
		assert(bigint::parse("        -1") == -1);
		assert(bigint::parse("-1", 10) == -1);
		assert(bigint::parse("7", 8) == -1);
		assert(bigint::parse("f", 16) == -1);
		assert(bigint::parse("-1", 10) == -1);
		assert(bigint::parse("7", 8) == -1);
		assert(bigint::parse("f", 16) == -1);
		assert(bigint::parse("-1", 0) == -1);
		assert(bigint::parse("0o7", 0) == -1);
		assert(bigint::parse("0xf", 0) == -1);
		assert(bigint::parse("-1", 0) == -1);
		assert(bigint::parse("0o7", 0) == -1);
		assert(bigint::parse("0xf", 0) == -1);
		assert(bigint::parse("        -1", 0) == -1);
		assert(bigint::parse("       0o7", 0) == -1);
		assert(bigint::parse("       0xf", 0) == -1);
	}
	{
		bigint v = -2;
		assert(sstostr(v) == "-2");
		assert(sstostr(std::setw(10), v) == "        -2");
		assert(sstostr(std::showpos, v) == "-2");
		assert(sstostr(std::showpos, std::setw(10), v) == "        -2");
		assert(sstostr(std::dec, v) == "-2");
		assert(sstostr(std::oct, v) == "6");
		assert(sstostr(std::hex, v) == "e");
		assert(sstostr(std::dec, std::showpos, v) == "-2");
		assert(sstostr(std::oct, std::showpos, v) == "6");
		assert(sstostr(std::hex, std::showpos, v) == "e");
		assert(sstostr(std::dec, std::showbase, v) == "-2");
		assert(sstostr(std::oct, std::showbase, v) == "0o6");
		assert(sstostr(std::hex, std::showbase, v) == "0xe");
		assert(sstostr(std::dec, std::showbase, std::showpos, v) == "-2");
		assert(sstostr(std::oct, std::showbase, std::showpos, v) == "0o6");
		assert(sstostr(std::hex, std::showbase, std::showpos, v) == "0xe");
		assert(sstostr(std::dec, std::showbase, std::showpos, std::setw(10), v) == "        -2");
		assert(sstostr(std::oct, std::showbase, std::showpos, std::setw(10), v) == "       0o6");
		assert(sstostr(std::hex, std::showbase, std::showpos, std::setw(10), v) == "       0xe");

		assert(tobin(v) == "10");
		assert(tobin_showbase(v) == "0b10");
		assert(tobin_uppercase(v) == "10");

		assert(frombin("10") == -2);
		assert(fromunk("0b10") == -2);

		assert(bigint::parse("-2") == -2);
		assert(bigint::parse("        -2") == -2);
		assert(bigint::parse("-2") == -2);
		assert(bigint::parse("        -2") == -2);
		assert(bigint::parse("-2", 10) = -2);
		assert(bigint::parse("6", 8) == -2);
		assert(bigint::parse("e", 16) == -2);
		assert(bigint::parse("-2", 10) == -2);
		assert(bigint::parse("6", 8) == -2);
		assert(bigint::parse("e", 16) == -2);
		assert(bigint::parse("-2", 0) == -2);
		assert(bigint::parse("0o6", 0) == -2);
		assert(bigint::parse("0xe", 0) == -2);
		assert(bigint::parse("-2", 0) == -2);
		assert(bigint::parse("0o6", 0) == -2);
		assert(bigint::parse("0xe", 0) == -2);
		assert(bigint::parse("        -2", 0) == -2);
		assert(bigint::parse("        0o6", 0) == -2);
		assert(bigint::parse("       0xe", 0) == -2);
	}
	{
		bigint v = -4;
		assert(sstostr(v) == "-4");
		assert(sstostr(std::setw(10), v) == "        -4");
		assert(sstostr(std::showpos, v) == "-4");
		assert(sstostr(std::showpos, std::setw(10), v) == "        -4");
		assert(sstostr(std::dec, v) == "-4");
		assert(sstostr(std::oct, v) == "4");
		assert(sstostr(std::hex, v) == "c");
		assert(sstostr(std::dec, std::showpos, v) == "-4");
		assert(sstostr(std::oct, std::showpos, v) == "4");
		assert(sstostr(std::hex, std::showpos, v) == "c");
		assert(sstostr(std::dec, std::showbase, v) == "-4");
		assert(sstostr(std::oct, std::showbase, v) == "0o4");
		assert(sstostr(std::hex, std::showbase, v) == "0xc");
		assert(sstostr(std::dec, std::showbase, std::showpos, v) == "-4");
		assert(sstostr(std::oct, std::showbase, std::showpos, v) == "0o4");
		assert(sstostr(std::hex, std::showbase, std::showpos, v) == "0xc");
		assert(sstostr(std::dec, std::showbase, std::showpos, std::setw(10), v) == "        -4");
		assert(sstostr(std::oct, std::showbase, std::showpos, std::setw(10), v) == "       0o4");
		assert(sstostr(std::hex, std::showbase, std::showpos, std::setw(10), v) == "       0xc");

		assert(tobin(v) == "100");
		assert(tobin_showbase(v) == "0b100");
		assert(tobin_uppercase(v) == "100");

		assert(frombin("100") == -4);
		assert(fromunk("0b100") == -4);

		assert(bigint::parse("-4") == -4);
		assert(bigint::parse("        -4") == -4);
		assert(bigint::parse("-4") == -4);
		assert(bigint::parse("        -4") == -4);
		assert(bigint::parse("-4", 10) == -4);
		assert(bigint::parse("4", 8) == -4);
		assert(bigint::parse("c", 16) == -4);
		assert(bigint::parse("-4", 10) == -4);
		assert(bigint::parse("4", 8) == -4);
		assert(bigint::parse("c", 16) == -4);
		assert(bigint::parse("-4", 0) == -4);
		assert(bigint::parse("0o4", 0) == -4);
		assert(bigint::parse("0xc", 0) == -4);
		assert(bigint::parse("-4", 0) == -4);
		assert(bigint::parse("0o4", 0) == -4);
		assert(bigint::parse("0xc", 0) == -4);
		assert(bigint::parse("        -4", 0) == -4);
		assert(bigint::parse("       0o4", 0) == -4);
		assert(bigint::parse("       0xc", 0) == -4);
	}
	{
		bigint v = -5;
		assert(sstostr(v) == "-5");
		assert(sstostr(std::setw(10), v) == "        -5");
		assert(sstostr(std::showpos, v) == "-5");
		assert(sstostr(std::showpos, std::setw(10), v) == "        -5");
		assert(sstostr(std::dec, v) == "-5");
		assert(sstostr(std::oct, v) == "73");
		assert(sstostr(std::hex, v) == "b");
		assert(sstostr(std::dec, std::showpos, v) == "-5");
		assert(sstostr(std::oct, std::showpos, v) == "73");
		assert(sstostr(std::hex, std::showpos, v) == "b");
		assert(sstostr(std::dec, std::showbase, v) == "-5");
		assert(sstostr(std::oct, std::showbase, v) == "0o73");
		assert(sstostr(std::hex, std::showbase, v) == "0xb");
		assert(sstostr(std::dec, std::showbase, std::showpos, v) == "-5");
		assert(sstostr(std::oct, std::showbase, std::showpos, v) == "0o73");
		assert(sstostr(std::hex, std::showbase, std::showpos, v) == "0xb");
		assert(sstostr(std::dec, std::showbase, std::showpos, std::setw(10), v) == "        -5");
		assert(sstostr(std::oct, std::showbase, std::showpos, std::setw(10), v) == "      0o73");
		assert(sstostr(std::hex, std::showbase, std::showpos, std::setw(10), v) == "       0xb");

		assert(tobin(v) == "1011");
		assert(tobin_showbase(v) == "0b1011");
		assert(tobin_uppercase(v) == "1011");

		assert(frombin("1011") == -5);
		assert(fromunk("0b1011") == -5);

		assert(bigint::parse("-5") == -5);
		assert(bigint::parse("        -5") == -5);
		assert(bigint::parse("-5") == -5);
		assert(bigint::parse("        -5") == -5);
		assert(bigint::parse("-5", 10) == -5);
		assert(bigint::parse("73", 8) == -5);
		assert(bigint::parse("b", 16) == -5);
		assert(bigint::parse("-5", 10) == -5);
		assert(bigint::parse("73", 8) == -5);
		assert(bigint::parse("b", 16) == -5);
		assert(bigint::parse("-5", 0) == -5);
		assert(bigint::parse("0o73", 0) == -5);
		assert(bigint::parse("0xb", 0) == -5);
		assert(bigint::parse("-5", 0) == -5);
		assert(bigint::parse("0o73", 0) == -5);
		assert(bigint::parse("0xb", 0) == -5);
		assert(bigint::parse("        -5", 0) == -5);
		assert(bigint::parse("      0o73", 0) == -5);
		assert(bigint::parse("       0xb", 0) == -5);
	}
	{
		bigint v = -16;
		assert(sstostr(v) == "-16");
		assert(sstostr(std::setw(10), v) == "       -16");
		assert(sstostr(std::showpos, v) == "-16");
		assert(sstostr(std::showpos, std::setw(10), v) == "       -16");
		assert(sstostr(std::dec, v) == "-16");
		assert(sstostr(std::oct, v) == "60");
		assert(sstostr(std::hex, v) == "f0");
		assert(sstostr(std::dec, std::showpos, v) == "-16");
		assert(sstostr(std::oct, std::showpos, v) == "60");
		assert(sstostr(std::hex, std::showpos, v) == "f0");
		assert(sstostr(std::dec, std::showbase, v) == "-16");
		assert(sstostr(std::oct, std::showbase, v) == "0o60");
		assert(sstostr(std::hex, std::showbase, v) == "0xf0");
		assert(sstostr(std::dec, std::showbase, std::showpos, v) == "-16");
		assert(sstostr(std::oct, std::showbase, std::showpos, v) == "0o60");
		assert(sstostr(std::hex, std::showbase, std::showpos, v) == "0xf0");
		assert(sstostr(std::dec, std::showbase, std::showpos, std::setw(10), v) == "       -16");
		assert(sstostr(std::oct, std::showbase, std::showpos, std::setw(10), v) == "      0o60");
		assert(sstostr(std::hex, std::showbase, std::showpos, std::setw(10), v) == "      0xf0");

		assert(tobin(v) == "10000");
		assert(tobin_showbase(v) == "0b10000");
		assert(tobin_uppercase(v) == "10000");

		assert(frombin("10000") == -16);
		assert(fromunk("0b10000") == -16);

		assert(bigint::parse("-16") == -16);
		assert(bigint::parse("       -16") == -16);
		assert(bigint::parse("-16") == -16);
		assert(bigint::parse("       -16") == -16);
		assert(bigint::parse("-16", 10) == -16);
		assert(bigint::parse("60", 8) == -16);
		assert(bigint::parse("f0", 16) == -16);
		assert(bigint::parse("-16", 10) == -16);
		assert(bigint::parse("60", 8) == -16);
		assert(bigint::parse("f0", 16) == -16);
		assert(bigint::parse("-16", 0) == -16);
		assert(bigint::parse("0o60", 0) == -16);
		assert(bigint::parse("0xf0", 0) == -16);
		assert(bigint::parse("-16", 0) == -16);
		assert(bigint::parse("0o60", 0) == -16);
		assert(bigint::parse("0xf0", 0) == -16);
		assert(bigint::parse("       -16", 0) == -16);
		assert(bigint::parse("      0o60", 0) == -16);
		assert(bigint::parse("      0xf0", 0) == -16);
	}
	{
		bigint v = 0xffffffffffffffffull;
		assert(sstostr(v) == "18446744073709551615");
		assert(sstostr(std::setw(10), v) == "18446744073709551615");
		assert(sstostr(std::showpos, v) == "+18446744073709551615");
		assert(sstostr(std::showpos, std::setw(10), v) == "+18446744073709551615");
		assert(sstostr(std::dec, v) == "18446744073709551615");
		assert(sstostr(std::oct, v) == "1777777777777777777777");
		assert(sstostr(std::hex, v) == "0ffffffffffffffff");
		assert(sstostr(std::dec, std::showpos, v) == "+18446744073709551615");
		assert(sstostr(std::oct, std::showpos, v) == "1777777777777777777777");
		assert(sstostr(std::hex, std::showpos, v) == "0ffffffffffffffff");
		assert(sstostr(std::dec, std::showbase, v) == "18446744073709551615");
		assert(sstostr(std::oct, std::showbase, v) == "0o1777777777777777777777");
		assert(sstostr(std::hex, std::showbase, v) == "0x0ffffffffffffffff");
		assert(sstostr(std::dec, std::showbase, std::showpos, v) == "+18446744073709551615");
		assert(sstostr(std::oct, std::showbase, std::showpos, v) == "0o1777777777777777777777");
		assert(sstostr(std::hex, std::showbase, std::showpos, v) == "0x0ffffffffffffffff");
		assert(sstostr(std::dec, std::showbase, std::showpos, std::setw(10), v) == "+18446744073709551615");
		assert(sstostr(std::oct, std::showbase, std::showpos, std::setw(10), v) == "0o1777777777777777777777");
		assert(sstostr(std::hex, std::showbase, std::showpos, std::setw(10), v) == "0x0ffffffffffffffff");
		assert(sstostr(std::dec, std::showbase, std::showpos, std::uppercase, std::setw(10), v) == "+18446744073709551615");
		assert(sstostr(std::oct, std::showbase, std::showpos, std::uppercase, std::setw(10), v) == "0o1777777777777777777777");
		assert(sstostr(std::hex, std::showbase, std::showpos, std::uppercase, std::setw(10), v) == "0x0FFFFFFFFFFFFFFFF");

		assert(tobin(v) == "01111111111111111111111111111111111111111111111111111111111111111");
		assert(tobin_showbase(v) == "0b01111111111111111111111111111111111111111111111111111111111111111");
		assert(tobin_uppercase(v) == "01111111111111111111111111111111111111111111111111111111111111111");

		assert(frombin("01111111111111111111111111111111111111111111111111111111111111111") == v);
		assert(fromunk("0b01111111111111111111111111111111111111111111111111111111111111111") == v);

		assert(bigint::parse("18446744073709551615") == v);
		assert(bigint::parse("18446744073709551615") == v);
		assert(bigint::parse("+18446744073709551615") == v);
		assert(bigint::parse("+18446744073709551615") == v);
		assert(bigint::parse("18446744073709551615", 10) == v);
		assert(bigint::parse("1777777777777777777777", 8) == v);
		assert(bigint::parse("0ffffffffffffffff", 16) == v);
		assert(bigint::parse("+18446744073709551615", 10) == v);
		assert(bigint::parse("1777777777777777777777", 8) == v);
		assert(bigint::parse("0ffffffffffffffff", 16) == v);
		assert(bigint::parse("18446744073709551615", 0) == v);
		assert(bigint::parse("0o1777777777777777777777", 0) == v);
		assert(bigint::parse("0x0ffffffffffffffff", 0) == v);
		assert(bigint::parse("+18446744073709551615", 0) == v);
		assert(bigint::parse("0o1777777777777777777777", 0) == v);
		assert(bigint::parse("0x0ffffffffffffffff", 0) == v);
		assert(bigint::parse("+18446744073709551615", 0) == v);
		assert(bigint::parse("0o1777777777777777777777", 0) == v);
		assert(bigint::parse("0x0ffffffffffffffff", 0) == v);
		assert(bigint::parse("+18446744073709551615", 0) == v);
		assert(bigint::parse("0o1777777777777777777777", 0) == v);
		assert(bigint::parse("0x0FFFFFFFFFFFFFFFF", 0) == v);
		assert(bigint::parse("0O1777777777777777777777", 0) == v);
		assert(bigint::parse("0X0FFFFFFFFFFFFFFFF", 0) == v);
	}
	{
		bigint v = 0777777777777777777777ull;
		assert(sstostr(v) == "9223372036854775807");
		assert(sstostr(std::setw(10), v) == "9223372036854775807");
		assert(sstostr(std::showpos, v) == "+9223372036854775807");
		assert(sstostr(std::showpos, std::setw(10), v) == "+9223372036854775807");
		assert(sstostr(std::dec, v) == "9223372036854775807");
		assert(sstostr(std::oct, v) == "0777777777777777777777");
		assert(sstostr(std::hex, v) == "7fffffffffffffff");
		assert(sstostr(std::dec, std::showpos, v) == "+9223372036854775807");
		assert(sstostr(std::oct, std::showpos, v) == "0777777777777777777777");
		assert(sstostr(std::hex, std::showpos, v) == "7fffffffffffffff");
		assert(sstostr(std::dec, std::showbase, v) == "9223372036854775807");
		assert(sstostr(std::oct, std::showbase, v) == "0o0777777777777777777777");
		assert(sstostr(std::hex, std::showbase, v) == "0x7fffffffffffffff");
		assert(sstostr(std::dec, std::showbase, std::showpos, v) == "+9223372036854775807");
		assert(sstostr(std::oct, std::showbase, std::showpos, v) == "0o0777777777777777777777");
		assert(sstostr(std::hex, std::showbase, std::showpos, v) == "0x7fffffffffffffff");
		assert(sstostr(std::dec, std::showbase, std::showpos, std::setw(10), v) == "+9223372036854775807");
		assert(sstostr(std::oct, std::showbase, std::showpos, std::setw(10), v) == "0o0777777777777777777777");
		assert(sstostr(std::hex, std::showbase, std::showpos, std::setw(10), v) == "0x7fffffffffffffff");

		assert(tobin(v) == "0111111111111111111111111111111111111111111111111111111111111111");
		assert(tobin_showbase(v) == "0b0111111111111111111111111111111111111111111111111111111111111111");
		assert(tobin_uppercase(v) == "0111111111111111111111111111111111111111111111111111111111111111");

		assert(frombin("0111111111111111111111111111111111111111111111111111111111111111") == v);
		assert(fromunk("0b0111111111111111111111111111111111111111111111111111111111111111") == v);

		assert(bigint::parse("9223372036854775807") == v);
		assert(bigint::parse("9223372036854775807") == v);
		assert(bigint::parse("+9223372036854775807") == v);
		assert(bigint::parse("+9223372036854775807") == v);
		assert(bigint::parse("9223372036854775807", 10) == v);
		assert(bigint::parse("0777777777777777777777", 8) == v);
		assert(bigint::parse("7fffffffffffffff", 16) == v);
		assert(bigint::parse("+9223372036854775807", 10) == v);
		assert(bigint::parse("0777777777777777777777", 8) == v);
		assert(bigint::parse("7fffffffffffffff", 16) == v);
		assert(bigint::parse("9223372036854775807", 0) == v);
		assert(bigint::parse("0o0777777777777777777777", 0) == v);
		assert(bigint::parse("0x7fffffffffffffff", 0) == v);
		assert(bigint::parse("+9223372036854775807", 0) == v);
		assert(bigint::parse("0o0777777777777777777777", 0) == v);
		assert(bigint::parse("0x7fffffffffffffff", 0) == v);
		assert(bigint::parse("+9223372036854775807", 0) == v);
		assert(bigint::parse("0o0777777777777777777777", 0) == v);
		assert(bigint::parse("0x7fffffffffffffff", 0) == v);
	}
	{
		bigint res;

		assert(!bigint::try_parse(res, "", 10));
		assert(!bigint::try_parse(res, "", 16));
		assert(!bigint::try_parse(res, "", 8));
		assert(!bigint::try_parse(res, "", 2));
		assert(!bigint::try_parse(res, "", 0));

		assert(!bigint::try_parse(res, " ", 10));
		assert(!bigint::try_parse(res, " ", 16));
		assert(!bigint::try_parse(res, " ", 8));
		assert(!bigint::try_parse(res, " ", 2));
		assert(!bigint::try_parse(res, " ", 0));

		assert(!bigint::try_parse(res, "0x", 0));
		assert(!bigint::try_parse(res, "0x ", 0));
		assert(!bigint::try_parse(res, "0x*", 0));
		assert(!bigint::try_parse(res, "0o", 0));
		assert(!bigint::try_parse(res, "0o ", 0));
		assert(!bigint::try_parse(res, "0o*", 0));
		assert(!bigint::try_parse(res, "0b", 0));
		assert(!bigint::try_parse(res, "0b ", 0));
		assert(!bigint::try_parse(res, "0b*", 0));
		assert(!bigint::try_parse(res, "+0x0", 0));
		assert(!bigint::try_parse(res, "+0o0", 0));
		assert(!bigint::try_parse(res, "+0b0", 0));
		assert(!bigint::try_parse(res, "-0x0", 0));
		assert(!bigint::try_parse(res, "-0o0", 0));
		assert(!bigint::try_parse(res, "-0b0", 0));
	}
	{
		assert(bigint::parse("0777", 0) == 777); // we don't use C-style prefixes (octal is "0o777")
	}
}
void block_parse_tests()
{
	{
		const char *strs[] =
		{
			"67039039649712985497870124991029230637396829102961966888617807218608820150367734884009371490834517138450159290932430254268769414059732849732168245030420",
			"6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042000",
			"6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042048"
		};
		const char *bigstr = "-58375738281737588279837147284739847598236414275903469456726347213103494759684960785976894576983754826435876982739387450956487904879065879047659834598264827454587609"
			"239859698657983479286431823749275039865490658759768596874957698345876438726347263598375469587698796587984769839827439028746598475684509845024376456890783984759834769846"
			"2938490345903845901798179038356859708045892374298734836938403943756845982379824983745984576509863984759826429374938474985675467348652984597659856798346598275876984576893";

		uint_t<512> u;
		bigint b;

		u = uint_t<512>::parse(strs[0]);
		assert(u.blocks[0] == 5165088340638674452ull);
		assert(u.blocks[1] == 1475739525896764129ull);
		assert(u.blocks[2] == 16233134784864405422ull);
		assert(u.blocks[3] == 12543785970122495098ull);
		assert(u.blocks[4] == 8854437155380584775ull);
		assert(u.blocks[5] == 5165088340638674452ull);
		assert(u.blocks[6] == 1475739525896764129ull);
		assert(u.blocks[7] == 92233720368547758ull);
		assert(sstostr(u) == strs[0]);
		u *= 100ull;
		assert(sstostr(u) == strs[1]);
		u += 48ull;
		assert(sstostr(u) == strs[2]);

		b = bigint::parse(strs[0]);
		assert(b.blocks.size() == 8);
		assert(b.blocks[0] == 5165088340638674452ull);
		assert(b.blocks[1] == 1475739525896764129ull);
		assert(b.blocks[2] == 16233134784864405422ull);
		assert(b.blocks[3] == 12543785970122495098ull);
		assert(b.blocks[4] == 8854437155380584775ull);
		assert(b.blocks[5] == 5165088340638674452ull);
		assert(b.blocks[6] == 1475739525896764129ull);
		assert(b.blocks[7] == 92233720368547758ull);
		assert(sstostr(b) == strs[0]);

		b = bigint::parse(bigstr);
		assert(sstostr(b) == bigstr);
	}
}
void overflow_parse_tests()
{
	{
		uint_t<512> u;
		int_t<512> s;

		assert(sstostr(uint_t<512>::parse("13407807929942597099574024998205846127479365820592393377723561443721764030073546976801874298166903427690031858186486050853753882811946569946433649006084095")) == "13407807929942597099574024998205846127479365820592393377723561443721764030073546976801874298166903427690031858186486050853753882811946569946433649006084095");
		assert(sstostr(int_t<512>::parse("-6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042048")) == "-6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042048");
		assert(sstostr(int_t<512>::parse("6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042047")) == "6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042047");

		assert_throws(std::invalid_argument, uint_t<512>::parse(u, "13407807929942597099574024998205846127479365820592393377723561443721764030073546976801874298166903427690031858186486050853753882811946569946433649006084096"));
		assert_throws(std::invalid_argument, int_t<512>::parse(s, "-6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042049"));
		assert_throws(std::invalid_argument, int_t<512>::parse(s, "6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042048"));

		assert(uint_t<512>::try_parse(u, "13407807929942597099574024998205846127479365820592393377723561443721764030073546976801874298166903427690031858186486050853753882811946569946433649006084095"));
		assert(int_t<512>::try_parse(s, "-6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042048"));
		assert(int_t<512>::try_parse(s, "6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042047"));

		assert(!uint_t<512>::try_parse(u, "13407807929942597099574024998205846127479365820592393377723561443721764030073546976801874298166903427690031858186486050853753882811946569946433649006084096"));
		assert(!int_t<512>::try_parse(s, "-6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042049"));
		assert(!int_t<512>::try_parse(s, "6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042048"));
	}
}
void promotion_demotion_tests()
{
	{
		auto f128 = [](auto v) {
			using namespace BiggerInts;
			typedef std::conditional_t<std::is_signed_v<decltype(v)>, int_t<128>, uint_t<128>> tt;

			tt val = v;
			tt val2(v);
			tt val3{ v };
			assert(val == v);
			assert(val2 == v);
			assert(val3 == v);
		};
		auto ff = [&](long long v)
		{
			f128((unsigned short)v);
			f128((unsigned int)v);
			f128((unsigned long)v);
			f128((unsigned long long)v);

			f128((short)v);
			f128((int)v);
			f128((long)v);
			f128((long long)v);
		};
		ff(46);
		ff(-46);
	}

	// -- demotion tests -- //

	{
		uint_t<256> u6 = 64;
		uint_t<128> u5 = (uint_t<128>)u6;
		uint_t<128> u5_1(u6);
		uint_t<128> u5_2{ u6 };
		u5 = (uint_t<128>)u6;
		unsigned long long u4 = (unsigned long long)u5;

		assert(u6 == 64u);
		assert(u5 == 64u);
		assert(u5_1 == 64u);
		assert(u5_2 == 64u);
		assert(u4 == 64u);

		int_t<256> s6 = -64;
		int_t<128> s5 = (int_t<128>)s6;
		int_t<128> s5_1(s6);
		int_t<128> s5_2{ s6 };
		long long s4 = (long long)s5;

		assert(s6 == -64);
		assert(s5 == -64);
		assert(s5_1 == -64);
		assert(s5_2 == -64);
		assert(s4 == -64);
	}

	// -- bigint tests -- //

	{
		uint_t<256> v = 0u;
		bigint big = v;
		assert(!detail::is_neg(big));
		assert(big.blocks.size() == 0);
		assert(big == 0);
		assert(detail::highest_set_bit(big) == 0);

		int_t<256> sv = -1;
		big = sv;
		assert(detail::is_neg(big));
		assert(big.blocks.size() == 1);
		assert(big == -1);
	}

	{
		bigint val(46);
		bigint val_1{ 46 };
		bigint val_2 = 46;

		assert((uint_t<256>)val == 46u);
		assert((unsigned long)val == 46);
		assert((unsigned int)val_1 == 46);
		assert((unsigned short)val_2 == 46);

		bigint val_3(-46);
		bigint val_4{ -46 };
		bigint val_5 = -46;

		assert((int_t<256>)val_3 == -46);
		assert((long)val_3 == -46);
		assert((int)val_4 == -46);
		assert((short)val_5 == -46);

		for (int i = 46; i >= -99; --i, --val)
		{
			assert(val == i);
		}
		for (int i = -46; i < 103; ++i, ++val_3)
		{
			assert(val_3 == i);
		}
	}

	{
		bigint v = int_t<256>::max();
		assert((int_t<128>)v == -1);
		assert((int_t<512>)v == int_t<256>::max());
	}

	{
		for (int i = 1; i < 102; ++i)
		{
			uint_t<256> v = i;
			assert(v == v);
			assert(!(v != v));
			assert(v <= v);
			assert(v >= v);
			assert(!(v < v));
			assert(!(v > v));

			uint_t<256> vm1 = v;
			--vm1;
			uint_t<256> vp1 = v;
			++vp1;

			assert(vm1 < v);
			assert(v < vp1);
			assert(vm1 < vp1);
			assert(vm1 <= v);
			assert(v <= vp1);
			assert(vm1 <= vp1);

			assert(vp1 > v);
			assert(v > vm1);
			assert(vp1 > vm1);
			assert(vp1 >= v);
			assert(v >= vm1);
			assert(vp1 >= vm1);
		}
		for (int i = -43; i < 56; ++i)
		{
			int_t<256> v = i;
			assert(v == v);
			assert(!(v != v));
			assert(v <= v);
			assert(v >= v);
			assert(!(v < v));
			assert(!(v > v));

			int_t<256> vm1 = v;
			--vm1;
			int_t<256> vp1 = v;
			++vp1;

			assert(vm1 < v);
			assert(v < vp1);
			assert(vm1 < vp1);
			assert(vm1 <= v);
			assert(v <= vp1);
			assert(vm1 <= vp1);

			assert(vp1 > v);
			assert(v > vm1);
			assert(vp1 > vm1);
			assert(vp1 >= v);
			assert(v >= vm1);
			assert(vp1 >= vm1);
		}
	}
}
void more_ops_tests()
{
	{
		long long vals[] = { 0, 2, -4, 1, 233, -1646, -233, 453567, -447453, 35, 2, -2, -1 };
		constexpr int count = sizeof(vals) / sizeof(*vals);
		for (int i = 0; i < count; ++i)
			for (int j = 0; j < count; ++j)
			{
				bigint sum = (bigint)vals[i] + (bigint)vals[j];
				assert(sum == vals[i] + vals[j]);
				assert(vals[i] + vals[j] == 0 ? sum.blocks.size() == 0 : sum.blocks.size() == 1);
				assert(sum <= sum);
				assert(sum >= sum);
				assert(sum == sum);
				assert(!(sum < sum));
				assert(!(sum > sum));
				assert(sum < sum + 1);
				assert(sum + 1 > sum);
				assert(sum > sum - 1);
				assert(sum - 1 < sum);
				assert(sstostr(sum) == sstostr(vals[i] + vals[j]));

				bigint notsum = ~sum;
				assert(notsum == -(vals[i] + vals[j]) - 1);
				assert(-(vals[i] + vals[j]) - 1 == 0 ? notsum.blocks.size() == 0 : notsum.blocks.size() == 1);
				assert(sstostr(notsum) == sstostr(-(vals[i] + vals[j]) - 1));

				bigint negsum = -sum;
				assert(negsum == -(vals[i] + vals[j]));
				assert(-(vals[i] + vals[j]) == 0 ? negsum.blocks.size() == 0 : negsum.blocks.size() == 1);
				assert(sstostr(negsum) == sstostr(-(vals[i] + vals[j])));

				assert(notsum + 1 == negsum);
				assert(sstostr(notsum + 1) == sstostr(negsum));

				bigint diff = (bigint)vals[i] - (bigint)vals[j];
				assert(diff == vals[i] - vals[j]);
				assert(vals[i] - vals[j] == 0 ? diff.blocks.size() == 0 : diff.blocks.size() == 1);
				assert(sstostr(diff) == sstostr(vals[i] - vals[j]));

				bigint prod = (bigint)vals[i] * (bigint)vals[j];
				assert(prod == vals[i] * vals[j]);
				assert(vals[i] * vals[j] == 0 ? prod.blocks.size() == 0 : prod.blocks.size() == 1);
				assert(sstostr(prod) == sstostr(vals[i] * vals[j]));

				if (vals[j] != 0)
				{
					bigint quo = (bigint)vals[i] / (bigint)vals[j];
					assert(quo == vals[i] / vals[j]);
					assert(vals[i] / vals[j] == 0 ? quo.blocks.size() == 0 : quo.blocks.size() == 1);
					assert(sstostr(quo) == sstostr(vals[i] / vals[j]));
				}

				bigint band = (bigint)vals[i] & (bigint)vals[j];
				assert(band == (vals[i] & vals[j]));
				assert((vals[i] & vals[j]) == 0 ? band.blocks.size() == 0 : band.blocks.size() == 1);
				assert(sstostr(band) == sstostr(vals[i] & vals[j]));

				bigint bor = (bigint)vals[i] | (bigint)vals[j];
				assert(bor == (vals[i] | vals[j]));
				assert((vals[i] | vals[j]) == 0 ? bor.blocks.size() == 0 : bor.blocks.size() == 1);
				assert(sstostr(bor) == sstostr(vals[i] | vals[j]));

				bigint bxor = (bigint)vals[i] ^ (bigint)vals[j];
				assert(bxor == (vals[i] ^ vals[j]));
				assert((vals[i] ^ vals[j]) == 0 ? bxor.blocks.size() == 0 : bxor.blocks.size() == 1);
				assert(sstostr(bxor) == sstostr(vals[i] ^ vals[j]));
			}

		{
			bigint bigpos = 0x8000000000000000ull;
			assert(bigpos.blocks.size() == 2);
			assert(bigpos.blocks[0] == 0x8000000000000000ull);
			assert(bigpos.blocks[1] == 0ull);
			assert(bigpos == bigpos);
		}

		{
			bigint bigpos = (bigint)std::numeric_limits<std::int64_t>::max() + (bigint)std::numeric_limits<std::int64_t>::max();
			assert(bigpos.blocks.size() == 2);
			assert(bigpos.blocks[0] == 0xfffffffffffffffeull);
			assert(bigpos.blocks[1] == 0);
		}

		{
			bigint min = (bigint)std::numeric_limits<std::int64_t>::min();
			assert(min.blocks.size() == 1);
			assert(min.blocks[0] == 0x8000000000000000ull);
			assert(min == std::numeric_limits<std::int64_t>::min());
		}

		{
			bigint bigneg = (bigint)0 - (bigint)std::numeric_limits<std::int64_t>::min();
			assert(bigneg.blocks.size() == 2);
			assert(bigneg.blocks[0] == 0x8000000000000000ull);
			assert(bigneg.blocks[1] == 0ull);
		}

		{
			bigint bigneg = (bigint)(-1) - (bigint)std::numeric_limits<std::int64_t>::min();
			assert(bigneg.blocks.size() == 1);
			assert(bigneg.blocks[0] == 0x7fffffffffffffffull);
		}

		{
			bigint bigneg = (bigint)(1) - (bigint)std::numeric_limits<std::int64_t>::min();
			assert(bigneg.blocks.size() == 2);
			assert(bigneg.blocks[0] == 0x8000000000000001ull);
			assert(bigneg.blocks[1] == 0ull);
		}

		{
			bigint bigneg = (bigint)(std::numeric_limits<std::int64_t>::min()) - (bigint)std::numeric_limits<std::int64_t>::min();
			assert(bigneg == 0);
			assert(bigneg.blocks.size() == 0);
		}

		{
			bigint bigneg = (bigint)(std::numeric_limits<std::int64_t>::min()) - (bigint)0x8000000000000000ull;
			assert(bigneg != 0);
			assert(!(bigneg == 0));
			assert(bigneg.blocks.size() == 2);
			assert(bigneg.blocks[0] == 0ull);
			assert(bigneg.blocks[1] == 0xffffffffffffffffull);
		}

		{
			bigint bigneg = (bigint)(std::numeric_limits<std::int64_t>::min()) - -(bigint)0x8000000000000000ull;
			assert(bigneg == 0);
			assert(!(bigneg != 0));
			assert(bigneg.blocks.size() == 0);
		}
	}
}
void big_promotion_tests()
{
	{
		bigint v = bigint::parse("0fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff", 16);
		bigint r = v * v;
		assert(sstostr(r) == "1026301351570233739499293428424157827113862320509474868877723027253195088579758740627615722253086934881590063514230345951830442386318246187269187562628311952630251852968692235997360173760639550065418089580517919902935726984364421732731351192069777176001110091235222448929542559681760328631249946802584338006151524151805763498018141131455523556616803713025");
		uint_t<1024> v1024;

		const char *bigint_strs[] =
		{
			"0fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
			"1fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
			"3fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
			"7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
			"0ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
			"1ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
			"3ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
			"7ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
			"0fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
		};
		const char *other_strs[] =
		{
			"fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
			"1fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
			"3fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
			"7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
			"ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
			"1ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
			"3ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
			"7ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
			"fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
		};
		const char *ans[] =
		{
			"52374249726338269920211035149241586435466272736689036631732661889538140742460318866977656952786520141288100957559302388064271034229234090852267300225025",
			"209496998905353079680844140596966345741865090946756146526930647558152562969870223490219956860001973311404575807200527048423494277926800759411047483310081",
			"837987995621412318723376562387865382967460363787024586107722590232610251879538790005498485537719678738122647182728743186026797393726931829648146498060289",
			"3351951982485649274893506249551461531869841455148098344430890360930441007518270952111231258346302285937499276638768242728772830138947184902600499121881089",
			"13407807929942597099574024998205846127479365820592393377723561443721764030073315392623399665776056285720014482370779510884422601683867654778417822746804225",
			"53631231719770388398296099992823384509917463282369573510894245774887056120293724738850547927885919426820092681114531123476352968991628449449702943505776641",
			"214524926879081553593184399971293538039669853129478294043576983099548224481175825292116090241107066275160440227720950653782737000478829458470875079060226049",
			"858099707516326214372737599885174152158679412517913176174307932398192897924705153841892158023555042236401899917409454934885598250939949155227626926315143169",
			"3432398830065304857490950399540696608634717650071652704697231729592771591698824320714424226212473723217127877682689124379051693501809059263598760925409050625",
		};
		static_assert(sizeof(bigint_strs) == sizeof(ans));
		static_assert(sizeof(other_strs) == sizeof(ans));
		constexpr std::size_t count = sizeof(bigint_strs) / sizeof(*bigint_strs);

		for (std::size_t i = 0; i < count; ++i)
		{
			bigint::parse(v, bigint_strs[i], 16);
			assert(sstostr(std::hex, v) == bigint_strs[i]);
			v *= v;
			assert(sstostr(v) == ans[i]);

			assert((v | 0) == v);
			assert((-v | 0) == -v);

			assert((0 | v) == v);
			assert((0 | -v) == -v);
		}

		for (std::size_t i = 0; i < count; ++i)
		{
			uint_t<1024>::parse(v1024, other_strs[i], 16);
			assert(sstostr(std::hex, v1024) == other_strs[i]);
			v1024 *= v1024;
			assert(sstostr(v1024) == ans[i]);
		}
	}
}
void misc_tests()
{
	{
		assert(uint_t<256>(5) == 5u);
		assert(int_t<256>(5) == 5);
		assert(bigint(5) == 5);

		assert(uint_t<256>(5) < uint_t<256>(6));
		assert(uint_t<256>(-5) > uint_t<256>(6));

		assert(int_t<256>(5) < int_t<256>(6));
		assert(int_t<256>(-5) < int_t<256>(6));
	}
	{
		assert(sstostr(bigint::parse("12") * bigint::parse("12")) == "144");
		assert(sstostr(bigint::parse("-12") * bigint::parse("12")) == "-144");
		assert(sstostr(bigint::parse("12") * bigint::parse("-12")) == "-144");
		assert(sstostr(bigint::parse("-12") * bigint::parse("-12")) == "144");

		bigint v = 12;
		v *= v;
		assert(v == 144);
		v = -12;
		v *= v;
		assert(v == 144);
	}
	{
		std::pair<bigint, bigint> v = detail::divmod((bigint)5, (bigint)3);
		assert(v.first == 1);
		assert(v.second == 2);
		v = detail::divmod((bigint)5, 3);
		assert(v.first == 1);
		assert(v.second == 2);
		v = detail::divmod(5, (bigint)3);
		assert(v.first == 1);
		assert(v.second == 2);
	}
	{
		uint_t<256> v;
		assert(v == 0u);
		assert(0u == v);
		assert(v == v);
		assert(cmp(v, 0u) == 0);
		assert(cmp(0u, v) == 0);
		assert(cmp(v, v) == 0);
	}
	{
		int_t<256> v;
		assert(v == 0);
		assert(0 == v);
		assert(v == v);
		assert(cmp(v, 0) == 0);
		assert(cmp(0, v) == 0);
		assert(cmp(v, v) == 0);
	}
	{
		bigint v;
		assert(v == 0);
		assert(0 == v);
		assert(v == v);
		assert(cmp(v, 0) == 0);
		assert(cmp(0, v) == 0);
		assert(cmp(v, v) == 0);
	}
}
void pow_tests()
{
	{
		assert(sstostr(bigint::pow(12, -1)) == "0");
		assert(sstostr(bigint::pow(12, 0)) == "1");
		assert(sstostr(bigint::pow(12, 1)) == "12");
		assert(sstostr(bigint::pow(12, 2)) == "144");
		assert(sstostr(bigint::pow(12, 3)) == "1728");
		assert(sstostr(bigint::pow(12, 4)) == "20736");

		assert(sstostr(bigint::pow(-12, -1)) == "0");
		assert(sstostr(bigint::pow(-12, 0)) == "1");
		assert(sstostr(bigint::pow(-12, 1)) == "-12");
		assert(sstostr(bigint::pow(-12, 2)) == "144");
		assert(sstostr(bigint::pow(-12, 3)) == "-1728");
		assert(sstostr(bigint::pow(-12, 4)) == "20736");

		assert(bigint::pow(122, 32).tostr() == "5801156497853265440881973185035019503869198362251315125887713673216");
		assert(bigint::pow(56, 73).tostr() == "41469229998734605907229666533990439080884762393847798172038084645287681671512979662987869210144861737361911662097025032274313216");
		assert(bigint::pow(75, 53).tostr(16) == "45e8826020c86b10db09a711f9702595ada5dd56f5544eeca9c33bdc2a8fd7c468f60c37edc9520351b");
		assert(bigint::pow(42, 42).tostr(2) == "010110010001100101000001011111001101011010100001000111011011110111110010111101000001001101100101011110111011110000000011101111001000010000101010110001110110111011110011100100110010011001000000000000000000000000000000000000000000");
		assert(bigint::pow(2, 400).tostr(2) == "010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
	}
}
void factorial_tests()
{
	{
		const char *factorial_tests[] = {
			"1",
			"1",
			"2",
			"6",
			"24",
			"120",
			"720",
			"5040",
			"40320",
			"362880",
			"3628800",
			"39916800",
			"479001600",
			"6227020800",
			"87178291200",
			"1307674368000",
			"20922789888000",
			"355687428096000",
			"6402373705728000",
			"121645100408832000",
			"2432902008176640000",
			"51090942171709440000",
			"1124000727777607680000",
			"25852016738884976640000",
			"620448401733239439360000",
			"15511210043330985984000000",
			"403291461126605635584000000",
			"10888869450418352160768000000",
			"304888344611713860501504000000",
			"8841761993739701954543616000000",
			"265252859812191058636308480000000",
			"8222838654177922817725562880000000",
			"263130836933693530167218012160000000",
			"8683317618811886495518194401280000000",
			"295232799039604140847618609643520000000",
			"10333147966386144929666651337523200000000",
			"371993326789901217467999448150835200000000",
			"13763753091226345046315979581580902400000000",
			"523022617466601111760007224100074291200000000",
			"20397882081197443358640281739902897356800000000",
			"815915283247897734345611269596115894272000000000",
			"33452526613163807108170062053440751665152000000000",
			"1405006117752879898543142606244511569936384000000000",
			"60415263063373835637355132068513997507264512000000000",
			"2658271574788448768043625811014615890319638528000000000",
			"119622220865480194561963161495657715064383733760000000000",
			"5502622159812088949850305428800254892961651752960000000000",
			"258623241511168180642964355153611979969197632389120000000000",
			"12413915592536072670862289047373375038521486354677760000000000",
			"608281864034267560872252163321295376887552831379210240000000000",
			"30414093201713378043612608166064768844377641568960512000000000000",
			"1551118753287382280224243016469303211063259720016986112000000000000",
			"80658175170943878571660636856403766975289505440883277824000000000000",
			"4274883284060025564298013753389399649690343788366813724672000000000000",
			"230843697339241380472092742683027581083278564571807941132288000000000000",
			"12696403353658275925965100847566516959580321051449436762275840000000000000",
			"710998587804863451854045647463724949736497978881168458687447040000000000000",
			"40526919504877216755680601905432322134980384796226602145184481280000000000000",
			"2350561331282878571829474910515074683828862318181142924420699914240000000000000",
			"138683118545689835737939019720389406345902876772687432540821294940160000000000000",
			"8320987112741390144276341183223364380754172606361245952449277696409600000000000000",
			"507580213877224798800856812176625227226004528988036003099405939480985600000000000000",
			"31469973260387937525653122354950764088012280797258232192163168247821107200000000000000",
			"1982608315404440064116146708361898137544773690227268628106279599612729753600000000000000",
			"126886932185884164103433389335161480802865516174545192198801894375214704230400000000000000",
			"8247650592082470666723170306785496252186258551345437492922123134388955774976000000000000000",
			"544344939077443064003729240247842752644293064388798874532860126869671081148416000000000000000",
			"36471110918188685288249859096605464427167635314049524593701628500267962436943872000000000000000",
			"2480035542436830599600990418569171581047399201355367672371710738018221445712183296000000000000000",
			"171122452428141311372468338881272839092270544893520369393648040923257279754140647424000000000000000",
			"11978571669969891796072783721689098736458938142546425857555362864628009582789845319680000000000000000",
			"850478588567862317521167644239926010288584608120796235886430763388588680378079017697280000000000000000",
			"61234458376886086861524070385274672740778091784697328983823014963978384987221689274204160000000000000000",
			"4470115461512684340891257138125051110076800700282905015819080092370422104067183317016903680000000000000000",
			"330788544151938641225953028221253782145683251820934971170611926835411235700971565459250872320000000000000000",
			"24809140811395398091946477116594033660926243886570122837795894512655842677572867409443815424000000000000000000",
			"1885494701666050254987932260861146558230394535379329335672487982961844043495537923117729972224000000000000000000",
			"145183092028285869634070784086308284983740379224208358846781574688061991349156420080065207861248000000000000000000",
			"11324281178206297831457521158732046228731749579488251990048962825668835325234200766245086213177344000000000000000000",
			"894618213078297528685144171539831652069808216779571907213868063227837990693501860533361810841010176000000000000000000",
			"71569457046263802294811533723186532165584657342365752577109445058227039255480148842668944867280814080000000000000000000",
			"5797126020747367985879734231578109105412357244731625958745865049716390179693892056256184534249745940480000000000000000000",
			"475364333701284174842138206989404946643813294067993328617160934076743994734899148613007131808479167119360000000000000000000",
			"39455239697206586511897471180120610571436503407643446275224357528369751562996629334879591940103770870906880000000000000000000",
			"3314240134565353266999387579130131288000666286242049487118846032383059131291716864129885722968716753156177920000000000000000000",
			"281710411438055027694947944226061159480056634330574206405101912752560026159795933451040286452340924018275123200000000000000000000",
			"24227095383672732381765523203441259715284870552429381750838764496720162249742450276789464634901319465571660595200000000000000000000",
			"2107757298379527717213600518699389595229783738061356212322972511214654115727593174080683423236414793504734471782400000000000000000000",
			"185482642257398439114796845645546284380220968949399346684421580986889562184028199319100141244804501828416633516851200000000000000000000",
			"16507955160908461081216919262453619309839666236496541854913520707833171034378509739399912570787600662729080382999756800000000000000000000",
			"1485715964481761497309522733620825737885569961284688766942216863704985393094065876545992131370884059645617234469978112000000000000000000000",
			"135200152767840296255166568759495142147586866476906677791741734597153670771559994765685283954750449427751168336768008192000000000000000000000",
			"12438414054641307255475324325873553077577991715875414356840239582938137710983519518443046123837041347353107486982656753664000000000000000000000",
			"1156772507081641574759205162306240436214753229576413535186142281213246807121467315215203289516844845303838996289387078090752000000000000000000000",
			"108736615665674308027365285256786601004186803580182872307497374434045199869417927630229109214583415458560865651202385340530688000000000000000000000",
			"10329978488239059262599702099394727095397746340117372869212250571234293987594703124871765375385424468563282236864226607350415360000000000000000000000",
			"991677934870949689209571401541893801158183648651267795444376054838492222809091499987689476037000748982075094738965754305639874560000000000000000000000",
			"96192759682482119853328425949563698712343813919172976158104477319333745612481875498805879175589072651261284189679678167647067832320000000000000000000000",
			"9426890448883247745626185743057242473809693764078951663494238777294707070023223798882976159207729119823605850588608460429412647567360000000000000000000000",
			"933262154439441526816992388562667004907159682643816214685929638952175999932299156089414639761565182862536979208272237582511852109168640000000000000000000000",
			"93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000",
		};
		constexpr std::size_t factorial_test_count = sizeof(factorial_tests) / sizeof(*factorial_tests);

		for (std::size_t i = 0; i < factorial_test_count; ++i)
		{
			bigint v = bigint::factorial(i);
			assert(sstostr(v) == factorial_tests[i]);
		}
	}
}
void divmod_sign_tests()
{
	static_assert(6 / 5 == 1); static_assert(6 % 5 == 1); // assert ISO2011 modulus requirements
	static_assert(6 / -5 == -1); static_assert(6 % -5 == 1);
	static_assert(-6 / 5 == -1); static_assert(-6 % 5 == -1);
	static_assert(-6 / -5 == 1); static_assert(-6 % -5 == -1);

	assert(detail::divmod(int_t<256>(6) COMMA int_t<256>(5)) == std::make_pair<int_t<256> COMMA int_t<256>>(1 COMMA 1));
	assert(detail::divmod(int_t<256>(6) COMMA int_t<256>(-5)) == std::make_pair<int_t<256> COMMA int_t<256>>(-1 COMMA 1));
	assert(detail::divmod(int_t<256>(-6) COMMA int_t<256>(5)) == std::make_pair<int_t<256> COMMA int_t<256>>(-1 COMMA -1));
	assert(detail::divmod(int_t<256>(-6) COMMA int_t<256>(-5)) == std::make_pair<int_t<256> COMMA int_t<256>>(1 COMMA -1));

	assert(detail::divmod(bigint(6) COMMA bigint(5)) == std::make_pair<bigint COMMA bigint>(1 COMMA 1));
	assert(detail::divmod(bigint(6) COMMA bigint(-5)) == std::make_pair<bigint COMMA bigint>(-1 COMMA 1));
	assert(detail::divmod(bigint(-6) COMMA bigint(5)) == std::make_pair<bigint COMMA bigint>(-1 COMMA - 1));
	assert(detail::divmod(bigint(-6) COMMA bigint(-5)) == std::make_pair<bigint COMMA bigint>(1 COMMA - 1));
}
void extract_termpos_tests()
{
	{
		std::istringstream ss("");
		uint_t<128> v;
		ss >> v;
		assert(!ss);
		assert(ss.get() == EOF);
	}
	{
		std::istringstream ss("");
		int_t<128> v;
		ss >> v;
		assert(!ss);
		assert(ss.get() == EOF);
	}
	{
		std::istringstream ss("\t");
		uint_t<128> v;
		ss >> v;
		assert(!ss);
		ss.clear();
		assert(ss.get() == EOF);
	}
	{
		std::istringstream ss("\t");
		int_t<128> v;
		ss >> v;
		assert(!ss);
		ss.clear();
		assert(ss.get() == EOF);
	}
	{
		std::istringstream ss("12");
		uint_t<128> v;
		ss >> v;
		assert((bool)ss);
		assert(ss.get() == EOF);
		assert(v == 12u);
	}
	{
		std::istringstream ss("12");
		int_t<128> v;
		ss >> v;
		assert((bool)ss);
		assert(ss.get() == EOF);
		assert(v == 12);
	}
	{
		std::istringstream ss("12 ");
		uint_t<128> v;
		ss >> v;
		assert((bool)ss);
		assert(ss.get() == ' ');
		assert(v == 12u);
	}
	{
		std::istringstream ss("12 ");
		int_t<128> v;
		ss >> v;
		assert((bool)ss);
		assert(ss.get() == ' ');
		assert(v == 12);
	}
	{
		std::istringstream ss("12j");
		uint_t<128> v;
		ss >> v;
		assert((bool)ss);
		assert(ss.get() == 'j');
		assert(v == 12u);
	}
	{
		std::istringstream ss("12j");
		int_t<128> v;
		ss >> v;
		assert((bool)ss);
		assert(ss.get() == 'j');
		assert(v == 12);
	}
	{
		std::istringstream ss("x12j");
		uint_t<128> v;
		ss >> v;
		assert(!ss);
		ss.clear();
		assert(ss.get() == 'x');
	}
	{
		std::istringstream ss(" x12j");
		int_t<128> v;
		ss >> v;
		assert(!ss);
		ss.clear();
		assert(ss.get() == 'x');
	}
	{
		std::istringstream ss("1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111010101");
		uint_t<128> v;
		ss >> std::oct >> v;
		assert(!ss);
		ss.clear();
		assert(ss.get() == EOF);
	}
	{
		std::istringstream ss("777777777777777777777777777777777777777777777777777777777777777777701234567");
		uint_t<128> v;
		ss >> std::oct >> v;
		assert(!ss);
		ss.clear();
		assert(ss.get() == EOF);
	}
	{
		std::istringstream ss("99999999999999999999999999999999999999999999999999999999999999999990123456789");
		uint_t<128> v;
		ss >> std::dec >> v;
		assert(!ss);
		ss.clear();
		assert(ss.get() == EOF);
	}
	{
		std::istringstream ss("fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff0123456789abcdefABCDEF");
		uint_t<128> v;
		ss >> std::hex >> v;
		assert(!ss);
		ss.clear();
		assert(ss.get() == EOF);
	}
	{
		std::istringstream ss("1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111010101 ");
		uint_t<128> v;
		ss >> std::oct >> v;
		assert(!ss);
		ss.clear();
		assert(ss.get() == ' ');
	}
	{
		std::istringstream ss("777777777777777777777777777777777777777777777777777777777777777777701234567 ");
		uint_t<128> v;
		ss >> std::oct >> v;
		assert(!ss);
		ss.clear();
		assert(ss.get() == ' ');
	}
	{
		std::istringstream ss("99999999999999999999999999999999999999999999999999999999999999999990123456789 ");
		uint_t<128> v;
		ss >> std::dec >> v;
		assert(!ss);
		ss.clear();
		assert(ss.get() == ' ');
	}
	{
		std::istringstream ss("fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff0123456789abcdefABCDEF ");
		uint_t<128> v;
		ss >> std::hex >> v;
		assert(!ss);
		ss.clear();
		assert(ss.get() == ' ');
	}
	{
		std::istringstream ss("11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111110101x");
		uint_t<128> v;
		ss >> std::oct >> v;
		assert(!ss);
		ss.clear();
		assert(ss.get() == 'x');
	}
	{
		std::istringstream ss("777777777777777777777777777777777777777777777777777777777777777777701234567x");
		uint_t<128> v;
		ss >> std::oct >> v;
		assert(!ss);
		ss.clear();
		assert(ss.get() == 'x');
	}
	{
		std::istringstream ss("99999999999999999999999999999999999999999999999999999999999999999990123456789x");
		uint_t<128> v;
		ss >> std::dec >> v;
		assert(!ss);
		ss.clear();
		assert(ss.get() == 'x');
	}
	{
		std::istringstream ss("fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff0123456789abcdefABCDEFx");
		uint_t<128> v;
		ss >> std::hex >> v;
		assert(!ss);
		ss.clear();
		assert(ss.get() == 'x');
	}
}
void compount_assignment()
{
	{
		uint_t<256> v = 0u;
		assert(v == 0u);

		v += 1; assert(v == 1u);
		v += 1u; assert(v == 2u);
		v += uint_t<128>{1}; assert(v == 3u);
		v += uint_t<256>{1}; assert(v == 4u);
		v += uint_t<512>{1}; assert(v == 5u);
		v += bigint{ 1 }; assert(v == 6u);

		v -= 1; assert(v == 5u);
		v -= 1u; assert(v == 4u);
		v -= uint_t<128>{1}; assert(v == 3u);
		v -= uint_t<256>{1}; assert(v == 2u);
		v -= uint_t<512>{1}; assert(v == 1u);
		v -= bigint{ 1 }; assert(v == 0u);
	}
	{
		int_t<256> v = 0;
		assert(v == 0);

		v += 1; assert(v == 1);
		v += 1u; assert(v == 2);
		v += uint_t<128>{1}; assert(v == 3);
		v += uint_t<256>{1}; assert(v == 4);
		v += uint_t<512>{1}; assert(v == 5);
		v += bigint{ 1 }; assert(v == 6);

		v -= 1; assert(v == 5);
		v -= 1u; assert(v == 4);
		v -= uint_t<128>{1}; assert(v == 3);
		v -= uint_t<256>{1}; assert(v == 2);
		v -= uint_t<512>{1}; assert(v == 1);
		v -= bigint{ 1 }; assert(v == 0);
	}
	{
		uint_t<256> v = 65u;
		assert(v == 65u);

		v &= 1; assert(v == 1u);
		v |= 2; assert(v == 3u);
		v &= 1u; assert(v == 1u);
		v |= 2u; assert(v == 3u);
		v ^= 1; assert(v == 2u);
		v ^= 1u; assert(v == 3u);

		v = 65u;
		v &= int_t<128>{1}; assert(v == 1u);
		v |= int_t<128>{2}; assert(v == 3u);
		v &= uint_t<128>{1}; assert(v == 1u);
		v |= uint_t<128>{2}; assert(v == 3u);
		v ^= int_t<128>{1}; assert(v == 2u);
		v ^= uint_t<128>{1}; assert(v == 3u);

		v = 65u;
		v &= int_t<256>{1}; assert(v == 1u);
		v |= int_t<256>{2}; assert(v == 3u);
		v &= uint_t<256>{1}; assert(v == 1u);
		v |= uint_t<256>{2}; assert(v == 3u);
		v ^= int_t<256>{1}; assert(v == 2u);
		v ^= uint_t<256>{1}; assert(v == 3u);

		v = 65u;
		v &= int_t<512>{1}; assert(v == 1u);
		v |= int_t<512>{2}; assert(v == 3u);
		v &= uint_t<512>{1}; assert(v == 1u);
		v |= uint_t<512>{2}; assert(v == 3u);
		v ^= int_t<512>{1}; assert(v == 2u);
		v ^= uint_t<512>{1}; assert(v == 3u);

		v = 65u;
		v &= bigint{1}; assert(v == 1u);
		v |= bigint{2}; assert(v == 3u);
		v &= bigint{1}; assert(v == 1u);
		v |= bigint{2}; assert(v == 3u);
		v ^= bigint{1}; assert(v == 2u);
		v ^= bigint{1}; assert(v == 3u);
	}
}

// --------------------------------------------

static_assert(std::is_same<uint_t<8>, std::uint8_t>::value, "type equivalence violation");
static_assert(std::is_same<uint_t<16>, std::uint16_t>::value, "type equivalence violation");
static_assert(std::is_same<uint_t<32>, std::uint32_t>::value, "type equivalence violation");
static_assert(std::is_same<uint_t<64>, std::uint64_t>::value, "type equivalence violation");

static_assert(std::is_same<int_t<8>, std::int8_t>::value, "type equivalence violation");
static_assert(std::is_same<int_t<16>, std::int16_t>::value, "type equivalence violation");
static_assert(std::is_same<int_t<32>, std::int32_t>::value, "type equivalence violation");
static_assert(std::is_same<int_t<64>, std::int64_t>::value, "type equivalence violation");

static_assert(std::is_same<std::common_type_t<int, long>, long>::value, "common type violation");
static_assert(std::is_same<std::common_type_t<int, int_t<256>>, int_t<256>>::value, "common type violation");
static_assert(std::is_same<std::common_type_t<int, bigint>, bigint>::value, "common type violation");
static_assert(std::is_same<std::common_type_t<int_t<256>, bigint>, bigint>::value, "common type violation");
static_assert(std::is_same<std::common_type_t<int_t<256>, int_t<512>>, int_t<512>>::value, "common type violation");

// --------------------------------------------

#define ASSERT_BIN_OP_RESULT_TYPE(OP, A, B, R) static_assert(std::is_same_v<decltype(std::declval<A>() OP std::declval<B>()), R>, "result type violation: " #A " " #OP " " #B " -> " #R)
#define ASSERT_BIN_OP_RESULT_TYPE_TESTS(OP) \
	ASSERT_BIN_OP_RESULT_TYPE(OP, uint_t<128>, uint_t<128>, uint_t<128>); \
	ASSERT_BIN_OP_RESULT_TYPE(OP, uint_t<128>, uint_t<256>, uint_t<256>); \
	ASSERT_BIN_OP_RESULT_TYPE(OP, uint_t<256>, uint_t<128>, uint_t<256>); \
	ASSERT_BIN_OP_RESULT_TYPE(OP, uint_t<256>, unsigned long long, uint_t<256>); \
	ASSERT_BIN_OP_RESULT_TYPE(OP, uint_t<256>, unsigned long, uint_t<256>); \
	ASSERT_BIN_OP_RESULT_TYPE(OP, uint_t<256>, unsigned int, uint_t<256>); \
	ASSERT_BIN_OP_RESULT_TYPE(OP, uint_t<256>, unsigned short, uint_t<256>); \
	ASSERT_BIN_OP_RESULT_TYPE(OP, int_t<128>, int_t<128>, int_t<128>); \
	ASSERT_BIN_OP_RESULT_TYPE(OP, int_t<128>, int_t<256>, int_t<256>); \
	ASSERT_BIN_OP_RESULT_TYPE(OP, int_t<256>, int_t<128>, int_t<256>); \
	ASSERT_BIN_OP_RESULT_TYPE(OP, int_t<256>, long long, int_t<256>); \
	ASSERT_BIN_OP_RESULT_TYPE(OP, int_t<256>, long, int_t<256>); \
	ASSERT_BIN_OP_RESULT_TYPE(OP, int_t<256>, int, int_t<256>); \
	ASSERT_BIN_OP_RESULT_TYPE(OP, int_t<256>, short, int_t<256>); \
	ASSERT_BIN_OP_RESULT_TYPE(OP, int_t<128>, bigint, bigint); \
	ASSERT_BIN_OP_RESULT_TYPE(OP, bigint, int_t<128>, bigint)
#define ASSERT_BIN_OP_RESULT_TYPE_TESTS_MIX(OP) \
	ASSERT_BIN_OP_RESULT_TYPE(OP, uint_t<256>, int_t<128>, uint_t<256>); \
	ASSERT_BIN_OP_RESULT_TYPE(OP, int_t<128>, uint_t<256>, uint_t<256>); \
	ASSERT_BIN_OP_RESULT_TYPE(OP, uint_t<128>, int_t<256>, int_t<256>); \
	ASSERT_BIN_OP_RESULT_TYPE(OP, int_t<256>, uint_t<128>, int_t<256>); \
	ASSERT_BIN_OP_RESULT_TYPE(OP, int_t<128>, uint_t<128>, uint_t<128>); \
	ASSERT_BIN_OP_RESULT_TYPE(OP, uint_t<128>, int_t<128>, uint_t<128>); \
	ASSERT_BIN_OP_RESULT_TYPE(OP, uint_t<128>, bigint, bigint); \
	ASSERT_BIN_OP_RESULT_TYPE(OP, bigint, uint_t<128>, bigint)
	
ASSERT_BIN_OP_RESULT_TYPE_TESTS(+);
ASSERT_BIN_OP_RESULT_TYPE_TESTS(-);

ASSERT_BIN_OP_RESULT_TYPE_TESTS(&);
ASSERT_BIN_OP_RESULT_TYPE_TESTS(|);
ASSERT_BIN_OP_RESULT_TYPE_TESTS(^);

ASSERT_BIN_OP_RESULT_TYPE_TESTS(*);
ASSERT_BIN_OP_RESULT_TYPE_TESTS(/);
ASSERT_BIN_OP_RESULT_TYPE_TESTS(%);

ASSERT_BIN_OP_RESULT_TYPE_TESTS_MIX(+);
ASSERT_BIN_OP_RESULT_TYPE_TESTS_MIX(-);

ASSERT_BIN_OP_RESULT_TYPE_TESTS_MIX(&);
ASSERT_BIN_OP_RESULT_TYPE_TESTS_MIX(|);
ASSERT_BIN_OP_RESULT_TYPE_TESTS_MIX(^);

// --------------------------------------------

int main()
{
	// -- core tests -- //

	detail_tests();
	smol_tests();
	bigboi_tests();
	numeric_limits_tests();
	simple_parse_tests();
	big_parse_tests();
	block_parse_tests();
	overflow_parse_tests();
	promotion_demotion_tests();
	more_ops_tests();
	big_promotion_tests();
	misc_tests();
	pow_tests();
	factorial_tests();
	divmod_sign_tests();
	extract_termpos_tests();
	compount_assignment();

	// -- benchmarks -- //

	benchmark_binary("multiply ", 200000, 20000, 1000, 5000, [](const auto &a, const auto &b) { return a * b; });
	benchmark_binary("divmod   ", 4000, 2000, 100, 500, [](const auto &a, const auto &b) { return detail::divmod(a, b); });
	benchmark_binary("tostr dec", 1000, 1000, 100, 100, [fmt = tostr_fmt{}](const auto &a, const auto &b) { return fmt(a) + fmt(b); });
	benchmark_binary("tostr hex", 1000, 1000, 1000, 1000, [fmt = tostr_fmt{ 16 }](const auto &a, const auto &b) { return fmt(a) + fmt(b); });
	benchmark_binary("tostr oct", 1000, 1000, 1000, 1000, [fmt = tostr_fmt{ 8 }](const auto &a, const auto &b) { return fmt(a) + fmt(b); });
	benchmark_binary("tostr bin", 1000, 1000, 1000, 1000, [fmt = tostr_fmt{ 2 }](const auto &a, const auto &b) { return fmt(a) + fmt(b); });
	benchmark_string(10000, 10000, 10000);

	// -- all tests completed -- //

	std::cerr << "all tests completed successfully\n";
	std::cin.get();
	return 0;
}