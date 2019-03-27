#define BOOST_TEST_MODULE test_name

#include <boost/test/unit_test.hpp>

#include <calib/data.hpp>

BOOST_AUTO_TEST_SUITE(suite_of_data_tests)

bool check_name(const CALIB::t_name &L, const CALIB::t_name &R) {

	if (L.size() != R.size()) return false;
	for (int i = 0; i < L.size(); ++ i) {
		if (L[i] == R[i]) continue;
		return false;
	}
	return true;
}

BOOST_AUTO_TEST_CASE(test_constructor) {

	BOOST_TEST_MESSAGE("Testing name constructors");

	CALIB::t_name A = "aBcXyZ";
	const char *P = "abcxyz";
	BOOST_REQUIRE(A.size() == 6);
	for (int i = 0; i < A.size(); ++ i) {
	BOOST_REQUIRE(A[i] == P[i]);
	}

	CALIB::t_name B("ABCXYZ");
	BOOST_REQUIRE(
	check_name(B, A)
	);

	std::string S = "abcxyz";
	CALIB::t_name C(S);
	BOOST_REQUIRE(
	check_name(C, B)
	);

	CALIB::t_name D(C);
	BOOST_REQUIRE(
	check_name(D, C)
	);
}

BOOST_AUTO_TEST_CASE(test_assignment) {

	BOOST_TEST_MESSAGE("Testing name assignment");

	CALIB::t_name A = "aBcXyZ";
	CALIB::t_name B;
	std::string S = "ABCXYZ";
	BOOST_REQUIRE(!check_name(B, A));
	B = "AbCxYz";
	BOOST_REQUIRE(check_name(B, A));
	B = S;
	BOOST_REQUIRE(check_name(B, A));
	A = "_";
	BOOST_REQUIRE(
	!check_name(B, A)
	);
	B = A;
	BOOST_REQUIRE(check_name(B, A));
	A = "";
	BOOST_REQUIRE(
	!check_name(B, A)
	);
}

BOOST_AUTO_TEST_CASE(test_compare) {

	BOOST_TEST_MESSAGE("Testing name compare");

	CALIB::t_name A = "ABCxyz";
	CALIB::t_name B = "abcXYZ";
	CALIB::t_name C = "abc";
	std::string S = "abcxyz";
	BOOST_REQUIRE(A == "abcxyz");
	BOOST_REQUIRE(A == S);
	BOOST_REQUIRE(A == B);
	BOOST_REQUIRE(A <= B);
	BOOST_REQUIRE(A >= B);
	BOOST_REQUIRE(C <= B);
	BOOST_REQUIRE(C < B);
	BOOST_REQUIRE(B >= B);
	BOOST_REQUIRE(B > C);
	BOOST_REQUIRE(A != C);
}

//...

BOOST_AUTO_TEST_SUITE_END()
