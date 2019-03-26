#define BOOST_TEST_MODULE test_data

#include <boost/test/unit_test.hpp>

#include <calib/data.hpp>

BOOST_AUTO_TEST_SUITE(suite_of_data_tests)

template <typename T>
bool check_data_2d(const CALIB::t_data<T> &L, const CALIB::t_data<T> &R) {

	if ((L.nrow() != R.nrow()) || (L.ncol() != R.ncol())) return false;
	for (int i = 0; i < L.nrow(); ++ i)
	for (int j = 0; j < L.ncol(); ++ j) {
		if (L(i, j) == R(i, j)) continue;
		return false;
	}
	return true;
}

template <typename T>
bool check_data_1d(const CALIB::t_data<T> &L, const CALIB::t_data<T> &R) {

	if (L.size() != R.size()) return false;
	for (int i = 0; i < L.size(); ++ i) {
		if (L(i) == R(i)) continue;
		return false;
	}
	return true;
}

BOOST_AUTO_TEST_CASE(test_constructor) {

	BOOST_TEST_MESSAGE("Testing data constructors");

	CALIB::t_data<CALIB::t_real> DAT2(5, 2);
	BOOST_REQUIRE(DAT2.data() != nullptr);
	BOOST_REQUIRE(DAT2.size() == 10);
	BOOST_REQUIRE(DAT2.nrow() == 5);
	BOOST_REQUIRE(DAT2.ncol() == 2);

	CALIB::t_data<CALIB::t_real> DAT1(10);
	BOOST_REQUIRE(DAT1.data() != nullptr);
	BOOST_REQUIRE(DAT1.size() == 10);
	BOOST_REQUIRE(DAT1.nrow() == 10);
	BOOST_REQUIRE(DAT1.ncol() == 1);

	CALIB::t_data<CALIB::t_real> DAT0;
	BOOST_REQUIRE(DAT0.data() == nullptr);
	BOOST_REQUIRE(DAT0.size() == 0);
	BOOST_REQUIRE(DAT0.nrow() == 0);
	BOOST_REQUIRE(DAT0.ncol() == 0);
}

BOOST_AUTO_TEST_CASE(test_resize) {

	BOOST_TEST_MESSAGE("Testing data resize");

	CALIB::t_data<CALIB::t_real> DATA;

	DATA.resize(5, 2);
	BOOST_REQUIRE(DATA.data() != nullptr);
	BOOST_REQUIRE(DATA.size() == 10);
	BOOST_REQUIRE(DATA.nrow() == 5);
	BOOST_REQUIRE(DATA.ncol() == 2);

	DATA.resize(10);
	BOOST_REQUIRE(DATA.data() != nullptr);
	BOOST_REQUIRE(DATA.size() == 10);
	BOOST_REQUIRE(DATA.nrow() == 10);
	BOOST_REQUIRE(DATA.ncol() == 1);
}

BOOST_AUTO_TEST_CASE(test_format) {

	BOOST_TEST_MESSAGE("Testing data format");
	//Column-major order is used:
	CALIB::t_data<CALIB::t_long> DAT2(5, 2);
	CALIB::t_long k = 0;
	for (int j = 0; j < DAT2.ncol(); ++ j)
	for (int i = 0; i < DAT2.nrow(); ++ i)
		DAT2(i, j) = k ++;
	
	CALIB::t_data<CALIB::t_long> DAT1(5, 2);
	for (int i = 0; i < DAT1.size(); ++ i)
		DAT1(i) = i;

	BOOST_REQUIRE(
	check_data_2d(DAT2, DAT1)
	);
}

//...

BOOST_AUTO_TEST_SUITE_END()
