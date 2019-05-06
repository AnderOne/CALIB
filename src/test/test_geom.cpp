#define BOOST_TEST_MODULE test_geom

#include <boost/test/unit_test.hpp>

#include <calib/geom.hpp>

BOOST_AUTO_TEST_SUITE(suite_of_geom_tests)

BOOST_AUTO_TEST_CASE(test_geom_sub) {

	BOOST_TEST_MESSAGE("Testing for difference between points");

	CALIB::t_geom<CALIB::PERIOD0, int> g00(1, 6, 2, 7);
	CALIB::t_geom<CALIB::PERIODX, int> g10(1, 6, 2, 7);
	CALIB::t_geom<CALIB::PERIODY, int> g01(1, 6, 2, 7);
	CALIB::t_geom<CALIB::PERIOD2, int> g11(1, 6, 2, 7);
	
	BOOST_REQUIRE(g00.subx(2, 5) == -3);
	BOOST_REQUIRE(g00.subx(5, 2) == 3);
	BOOST_REQUIRE(g00.suby(3, 6) == -3);
	BOOST_REQUIRE(g00.suby(6, 3) == 3);

	BOOST_REQUIRE(g10.subx(2, 5) == 2);
	BOOST_REQUIRE(g10.subx(5, 2) == -2);
	BOOST_REQUIRE(g10.suby(3, 6) == -3);
	BOOST_REQUIRE(g10.suby(6, 3) == 3);

	BOOST_REQUIRE(g01.subx(2, 5) == -3);
	BOOST_REQUIRE(g01.subx(5, 2) == 3);
	BOOST_REQUIRE(g01.suby(3, 6) == 2);
	BOOST_REQUIRE(g01.suby(6, 3) == -2);

	BOOST_REQUIRE(g11.subx(2, 5) == 2);
	BOOST_REQUIRE(g11.subx(5, 2) == -2);
	BOOST_REQUIRE(g11.suby(3, 6) == 2);
	BOOST_REQUIRE(g11.suby(6, 3) == -2);

}

BOOST_AUTO_TEST_CASE(test_geom_mov) {

	BOOST_TEST_MESSAGE("Testing for shifting of points");

	CALIB::t_geom<CALIB::PERIOD0, int> g00(1, 6, 2, 7);
	CALIB::t_geom<CALIB::PERIODX, int> g10(1, 6, 2, 7);
	CALIB::t_geom<CALIB::PERIODY, int> g01(1, 6, 2, 7);
	CALIB::t_geom<CALIB::PERIOD2, int> g11(1, 6, 2, 7);
	
	BOOST_REQUIRE(g00.movx(0) == 1);
	BOOST_REQUIRE(g00.movx(3) == 3);
	BOOST_REQUIRE(g00.movx(7) == 6);
	BOOST_REQUIRE(g00.movy(1) == 2);
	BOOST_REQUIRE(g00.movy(4) == 4);
	BOOST_REQUIRE(g00.movy(8) == 7);

	BOOST_REQUIRE(g10.movx(0) == 5);
	BOOST_REQUIRE(g10.movx(3) == 3);
	BOOST_REQUIRE(g10.movx(7) == 2);
	BOOST_REQUIRE(g10.movy(1) == 2);
	BOOST_REQUIRE(g10.movy(4) == 4);
	BOOST_REQUIRE(g10.movy(8) == 7);

	BOOST_REQUIRE(g01.movx(0) == 1);
	BOOST_REQUIRE(g01.movx(3) == 3);
	BOOST_REQUIRE(g01.movx(7) == 6);
	BOOST_REQUIRE(g01.movy(1) == 6);
	BOOST_REQUIRE(g01.movy(4) == 4);
	BOOST_REQUIRE(g01.movy(8) == 3);

	BOOST_REQUIRE(g11.movx(0) == 5);
	BOOST_REQUIRE(g11.movx(3) == 3);
	BOOST_REQUIRE(g11.movx(7) == 2);
	BOOST_REQUIRE(g11.movy(1) == 6);
	BOOST_REQUIRE(g11.movy(4) == 4);
	BOOST_REQUIRE(g11.movy(8) == 3);

}

BOOST_AUTO_TEST_CASE(test_cond) {

	BOOST_TEST_MESSAGE("Testing of boundary conditions");

	CALIB::t_cond c00 = CALIB::PERIOD0;
	CALIB::t_cond c10 = CALIB::PERIODX;
	CALIB::t_cond c01 = CALIB::PERIODY;
	CALIB::t_cond c11 = CALIB::PERIOD2;

	BOOST_REQUIRE(!IS_PERIODX(c00));
	BOOST_REQUIRE(!IS_PERIODX(c01));
	BOOST_REQUIRE(IS_PERIODX(c10));
	BOOST_REQUIRE(IS_PERIODX(c11));

	BOOST_REQUIRE(!IS_PERIODY(c00));
	BOOST_REQUIRE(!IS_PERIODY(c10));
	BOOST_REQUIRE(IS_PERIODY(c01));
	BOOST_REQUIRE(IS_PERIODY(c11));

}

//...

BOOST_AUTO_TEST_SUITE_END()
