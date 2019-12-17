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
#define t_end(name, cnt) { auto _dt = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - _t).count(); std::cerr << name ": " << _dt << "ms (" << cnt << " iterations)\n"; }

#define t_test(init, expr, count) { init; t_start; for(int i = 0; i < count; ++i) (expr); t_end; t_print(expr, count); }

template<typename ...Args>
std::string tostr(const Args &...args)
{
	std::ostringstream ostr;
	(ostr << ... << args);
	return ostr.str();
}

template<typename BinaryFunction>
void benchmark_binary(const char *name, std::size_t count_1, std::size_t count_2, std::size_t count_3, BinaryFunction f)
{
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
				dd = f(v, val2);
				dd = f(val2, v);
				dd = f(v, v);
				dd = f(val2, val2);
			}
		}
		t_end(name << " unsigned 128", count_1);
	}

	{
		uint_t<256> vals[] =
		{
			uint_t<256>::parse("16495386453214655388564132657894532123546789897856532123548786534241323564132"),
			uint_t<256>::parse("68573958787958483948687979685743456879789253125456795321234537897865149088795"),
			uint_t<256>::parse("47785319425432789132452387662354568686531231234567898998686525252521750040594"),
			uint_t<256>::parse("62641345234867595653789352121111040400456054511231243235464748687576678624112"),
			uint_t<256>::parse("12642723264275606435798655323342606031302620024626783175264211264237286172641"),
			uint_t<256>::parse("41672765971459452836324591059406343537667837237611632347685660909991442640942"),
			uint_t<256>::parse("46176276219426832272141594256267006162343037603776858790046132642672704216412"),
			uint_t<256>::parse("16126745353512706456790768664624426526567679254020425698654214787864567825452"),
			uint_t<256>::parse("13786541237897643215480666435034357737389796460345034303430457675429731824613"),
			uint_t<256>::parse("45389783215348783514060649676056783783378666546506443031332264316267678795642"),
			uint_t<256>::parse("26461343345768768645949213246533037502642462567285494136406726534234357678665"),
			uint_t<256>::parse("45676597896543534326162612343506645726780246243623123126452738607954037367866"),
			uint_t<256>::parse("16435738689456523313264560646727382642602642564253728672940256457898883765642"),
			uint_t<256>::parse("46567383388656261334357573830646507628607249504146296456425675972608366565421"),
		};

		uint_t<256> val2 = uint_t<256>::parse("16579862451765276262462456641541236525214");
		decltype(f(vals[0], vals[0])) dd;

		t_start;
		for (std::size_t i = 0; i < count_2; ++i)
		{
			for (const auto &v : vals)
			{
				dd = f(v, val2);
				dd = f(val2, v);
				dd = f(v, v);
				dd = f(val2, val2);
			}
		}
		t_end(name << " unsigned 256", count_2);
	}

	{
		uint_t<8192> vals[] =
		{
			uint_t<8192>::parse("16495386453214655388564132657894532123546789897856532123548786534241323564132"),
			uint_t<8192>::parse("68573958787958483948687979685743456879789253125456795321234537897865149088795"),
			uint_t<8192>::parse("47785319425432789132452387662354568686531231234567898998686525252521750040594"),
			uint_t<8192>::parse("62641345234867595653789352121111040400456054511231243235464748687576678624112"),
			uint_t<8192>::parse("12642723264275606435798655323342606031302620024626783175264211264237286172641"),
			uint_t<8192>::parse("41672765971459452836324591059406343537667837237611632347685660909991442640942"),
			uint_t<8192>::parse("46176276219426832272141594256267006162343037603776858790046132642672704216412"),
			uint_t<8192>::parse("16126745353512706456790768664624426526567679254020425698654214787864567825452"),
			uint_t<8192>::parse("13786541237897643215480666435034357737389796460345034303430457675429731824613"),
			uint_t<8192>::parse("45389783215348783514060649676056783783378666546506443031332264316267678795642"),
			uint_t<8192>::parse("26461343345768768645949213246533037502642462567285494136406726534234357678665"),
			uint_t<8192>::parse("45676597896543534326162612343506645726780246243623123126452738607954037367866"),
			uint_t<8192>::parse("16435738689456523313264560646727382642602642564253728672940256457898883765642"),
			uint_t<8192>::parse("46567383388656261334357573830646507628607249504146296456425675972608366565421"),
		};

		uint_t<8192> val2 = uint_t<8192>::parse("16579862451765276262462456641541236525214");
		decltype(f(vals[0], vals[0])) dd;

		t_start;
		for (std::size_t i = 0; i < count_3; ++i)
		{
			for (const auto &v : vals)
			{
				dd = f(v, val2);
				dd = f(val2, v);
				dd = f(v, v);
				dd = f(val2, val2);
			}
		}
		t_end(name << " unsigned 8192 (small)", count_3);
	}

	std::cerr << '\n';
}

int main(int argc, const char *argv[])
{
	assert(detail::highest_set_bit(0) == 0);
	for (std::size_t i = 0; i < 64; ++i)
	{
		assert(detail::highest_set_bit((std::uint64_t)1 << i) == i);
	}
	assert(detail::highest_set_bit(uint_t<512>(0)) == 0);

	assert(tostr(uint_t<128>(2) * uint_t<128>(13)) == "26");
	assert(tostr(uint_t<128>(2) + uint_t<128>(13)) == "15");
	assert(tostr(uint_t<128>(20) - uint_t<128>(13)) == "7");
	assert(tostr(uint_t<128>(2) / uint_t<128>(13)) == "0");
	assert(tostr(uint_t<128>(2) % uint_t<128>(13)) == "2");

	assert(tostr(2u * uint_t<128>(13)) == "26");
	assert(tostr(2u + uint_t<128>(13)) == "15");
	assert(tostr(20u - uint_t<128>(13)) == "7");
	assert(tostr(uint_t<128>(2) / 13u) == "0");
	assert(tostr(uint_t<128>(2) % 13u) == "2");

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
	big_mul <<= 73;
	assert(tostr(big_mul) == "1701978554986102599652944396520696900597046164730926970313320209887640881966509866088480298291994"
		"7794655263581165454491067120838946815593625140970325732225722954705464939198103933056536439508186725987781864851905868153749504");

	big_mul += big_mul;
	assert(tostr(big_mul) == "34039571099722051993058887930413938011940923294618539406266404197752817639330197321769605965839895589"
		"310527162330908982134241677893631187250281940651464451445909410929878396207866113072879016373451975563729703811736307499008");

	big_mul <<= 1;
	assert(tostr(big_mul) == "68079142199444103986117775860827876023881846589237078812532808395505635278660394643539211931679791178"
		"621054324661817964268483355787262374500563881302928902891818821859756792415732226145758032746903951127459407623472614998016");

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

	// -- test string parsing -- //

	assert(tostr(uint_t<256>::parse("184738578493786576848368767654873647697975746763857664")) == "184738578493786576848368767654873647697975746763857664");
	assert(tostr(uint_t<256>::parse("18")) == "18");
	assert(uint_t<256>::parse("18") == 18u);

	assert(tostr(int_t<256>::parse("184738578493786576848368767654873647697975746763857664")) == "184738578493786576848368767654873647697975746763857664");
	assert(tostr(int_t<256>::parse("18")) == "18");
	assert(int_t<256>::parse("18") == 18);
	assert(int_t<256>::parse("+18") == 18);

	assert(tostr(int_t<256>::parse("-184738578493786576848368767654873647697975746763857664")) == "-184738578493786576848368767654873647697975746763857664");
	assert(tostr(int_t<256>::parse("-18")) == "-18");
	assert(int_t<256>::parse("-18") == -18);

	assert(tostr(int_t<256>::parse(" -18")) == "-18");
	assert(tostr(int_t<256>::parse("-18 ")) == "-18");
	assert(tostr(int_t<256>::parse(" -18 ")) == "-18");
	assert(int_t<256>::parse("-18") == -18);
	assert(int_t<256>::parse("-18 ") == -18);
	assert(int_t<256>::parse(" -18") == -18);

	assert_throws(std::invalid_argument, uint_t<512>::parse("- 18"));
	assert_throws(std::invalid_argument, uint_t<512>::parse("+ 18"));

	assert_throws(std::invalid_argument, uint_t<512>::parse("abc"));
	assert_throws(std::invalid_argument, uint_t<512>::parse(" abc"));
	assert_throws(std::invalid_argument, uint_t<512>::parse("abc "));
	assert_throws(std::invalid_argument, uint_t<512>::parse(" abc "));

	assert_throws(std::invalid_argument, uint_t<512>::parse("a18"));
	assert_throws(std::invalid_argument, uint_t<512>::parse("18a"));

	assert_throws(std::invalid_argument, uint_t<512>::parse("a 18"));
	assert_throws(std::invalid_argument, uint_t<512>::parse("18 a"));

	// ----------

	assert_throws(std::invalid_argument, int_t<512>::parse("- 18"));
	assert_throws(std::invalid_argument, int_t<512>::parse("+ 18"));

	assert_throws(std::invalid_argument, int_t<512>::parse("abc"));
	assert_throws(std::invalid_argument, int_t<512>::parse(" abc"));
	assert_throws(std::invalid_argument, int_t<512>::parse("abc "));
	assert_throws(std::invalid_argument, int_t<512>::parse(" abc "));

	assert_throws(std::invalid_argument, int_t<512>::parse("a18"));
	assert_throws(std::invalid_argument, int_t<512>::parse("18a"));

	assert_throws(std::invalid_argument, int_t<512>::parse("a 18"));
	assert_throws(std::invalid_argument, int_t<512>::parse("18 a"));

	// base parsing tests

	assert_throws(std::invalid_argument, uint_t<256>::parse("0x764"));
	assert_throws(std::invalid_argument, uint_t<256>::parse("0x764", 10));
	assert_throws(std::invalid_argument, uint_t<256>::parse("0x764", 8));
	assert_throws(std::invalid_argument, uint_t<256>::parse("0x764", 16));
	assert(tostr(uint_t<256>::parse("0x764", 0)) == "1892");
	assert(tostr(uint_t<256>::parse("0764", 0)) == "500");
	assert(tostr(uint_t<256>::parse("764", 0)) == "764");

	assert(tostr(uint_t<256>::parse("167345562314657356422211643535261643535621", 10)) == "167345562314657356422211643535261643535621");
	assert(tostr(uint_t<256>::parse("167345562314657356422211643535261643535621", 8)) == "19846815955032434075951735996930767761");
	assert(tostr(uint_t<256>::parse("167345562314657356422211643535261643535621", 16)) == "32811116214936888653588305588611479218962394469921");

	assert(tostr(int_t<256>::parse("167345562314657356422211643535261643535621", 10)) == "167345562314657356422211643535261643535621");
	assert(tostr(int_t<256>::parse("167345562314657356422211643535261643535621", 8)) == "19846815955032434075951735996930767761");
	assert(tostr(int_t<256>::parse("167345562314657356422211643535261643535621", 16)) == "32811116214936888653588305588611479218962394469921");

	assert(tostr(int_t<256>::parse("+167345562314657356422211643535261643535621", 10)) == "167345562314657356422211643535261643535621");
	assert(tostr(int_t<256>::parse("+167345562314657356422211643535261643535621", 8)) == "19846815955032434075951735996930767761");
	assert(tostr(int_t<256>::parse("+167345562314657356422211643535261643535621", 16)) == "32811116214936888653588305588611479218962394469921");

	assert(tostr(int_t<256>::parse("-167345562314657356422211643535261643535621", 10)) == "-167345562314657356422211643535261643535621");
	assert(tostr(int_t<256>::parse("-167345562314657356422211643535261643535621", 8)) == "-19846815955032434075951735996930767761");
	assert(tostr(int_t<256>::parse("-167345562314657356422211643535261643535621", 16)) == "-32811116214936888653588305588611479218962394469921");

	{
		uint_t<256> temp;

		assert(!uint_t<256>::try_parse(temp, "0x764"));
		assert(!uint_t<256>::try_parse(temp, "0x764", 10));
		assert(!uint_t<256>::try_parse(temp, "0x764", 8));
		assert(!uint_t<256>::try_parse(temp, "0x764", 16));
		assert(uint_t<256>::try_parse(temp, "0x764", 0));
	}

	{
		int_t<256> temp;
		assert(int_t<256>::try_parse(temp, "-167345562314657356422211643535261643535621", 10));
		assert(tostr(temp) == "-167345562314657356422211643535261643535621");

		assert(int_t<256>::try_parse(temp, "-167345562314657356422211643535261643535621", 8));
		assert(tostr(temp) == "-19846815955032434075951735996930767761");

		assert(int_t<256>::try_parse(temp, "-167345562314657356422211643535261643535621", 16));
		assert(tostr(temp) == "-32811116214936888653588305588611479218962394469921");
	}

	{
		uint_t<512> u;

		u = uint_t<512>::parse("67039039649712985497870124991029230637396829102961966888617807218608820150367734884009371490834517138450159290932430254268769414059732849732168245030420");
		assert(tostr(u) == "67039039649712985497870124991029230637396829102961966888617807218608820150367734884009371490834517138450159290932430254268769414059732849732168245030420");
		u *= 100ull;
		assert(tostr(u) == "6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042000");
		u += 48ull;
		assert(tostr(u) == "6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042048");
	}

	// overflow parsing tests

	{
		uint_t<512> u;
		int_t<512> s;

		assert(tostr(uint_t<512>::parse("13407807929942597099574024998205846127479365820592393377723561443721764030073546976801874298166903427690031858186486050853753882811946569946433649006084095")) == "13407807929942597099574024998205846127479365820592393377723561443721764030073546976801874298166903427690031858186486050853753882811946569946433649006084095");
		assert(tostr(int_t<512>::parse("-6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042048")) == "-6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042048");
		assert(tostr(int_t<512>::parse("6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042047")) == "6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042047");

		assert_throws(std::invalid_argument, uint_t<512>::parse("13407807929942597099574024998205846127479365820592393377723561443721764030073546976801874298166903427690031858186486050853753882811946569946433649006084096"));
		assert_throws(std::invalid_argument, int_t<512>::parse("-6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042049"));
		assert_throws(std::invalid_argument, int_t<512>::parse("6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042048"));

		assert(uint_t<512>::try_parse(u, "13407807929942597099574024998205846127479365820592393377723561443721764030073546976801874298166903427690031858186486050853753882811946569946433649006084095"));
		assert(int_t<512>::try_parse(s, "-6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042048"));
		assert(int_t<512>::try_parse(s, "6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042047"));

		assert(!uint_t<512>::try_parse(u, "13407807929942597099574024998205846127479365820592393377723561443721764030073546976801874298166903427690031858186486050853753882811946569946433649006084096"));
		assert(!int_t<512>::try_parse(s, "-6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042049"));
		assert(!int_t<512>::try_parse(s, "6703903964971298549787012499102923063739682910296196688861780721860882015036773488400937149083451713845015929093243025426876941405973284973216824503042048"));
	}

	// -- promotion tests -- //

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

		assert((unsigned long)val == 46);
		assert((unsigned int)val_1 == 46);
		assert((unsigned short)val_2 == 46);

		bigint val_3(-46);
		bigint val_4{ -46 };
		bigint val_5 = -46;

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

				bigint notsum = ~sum;
				assert(notsum == -(vals[i] + vals[j]) - 1);
				assert(-(vals[i] + vals[j]) - 1 == 0 ? notsum.blocks.size() == 0 : notsum.blocks.size() == 1);

				bigint negsum = -sum;
				assert(negsum == -(vals[i] + vals[j]));
				assert(-(vals[i] + vals[j]) == 0 ? negsum.blocks.size() == 0 : negsum.blocks.size() == 1);

				assert(notsum + 1 == negsum);

				bigint diff = (bigint)vals[i] - (bigint)vals[j];
				assert(diff == vals[i] - vals[j]);
				assert(vals[i] - vals[j] == 0 ? diff.blocks.size() == 0 : diff.blocks.size() == 1);

				bigint prod = (bigint)vals[i] * (bigint)vals[j];
				assert(prod == vals[i] * vals[j]);
				assert(vals[i] * vals[j] == 0 ? prod.blocks.size() == 0 : prod.blocks.size() == 1);

				if (vals[j] != 0)
				{
					bigint quo = (bigint)vals[i] / (bigint)vals[j];
					assert(quo == vals[i] / vals[j]);
					assert(vals[i] / vals[j] == 0 ? quo.blocks.size() == 0 : quo.blocks.size() == 1);
				}
			}

		{
			bigint bigpos = (bigint)0x8000000000000000ull;
			assert(bigpos.blocks.size() == 2);
			assert(bigpos.blocks[0] == 0x8000000000000000ull);
			assert(bigpos.blocks[1] == 0ull);
			assert(bigpos == 0x8000000000000000ull);
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
	
	// -- benchmarks -- //

	benchmark_binary("divmod  ", 40000, 20000, 1000, [](const auto &a, const auto &b) { return detail::divmod(a, b); });
	benchmark_binary("multiply", 100000, 20000, 1000, [](const auto &a, const auto &b) { return a * b; });

	// -- all tests completed -- //

	std::cerr << "all tests completed successfully\n";
	std::cin.get();
	return 0;
}
