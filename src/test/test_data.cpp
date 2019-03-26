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

	CALIB::t_data<CALIB::t_real> A(10);
	for (int i = 0; i < A.size(); ++ i) { A(i) = i; }
	CALIB::t_data<CALIB::t_real> B(A);
	BOOST_REQUIRE(
	(B.data() != A.data()) && check_data_1d(B, A)
	);

	CALIB::t_data<CALIB::t_real> C(std::move(A));
	BOOST_REQUIRE(A.data() == 0);
	BOOST_REQUIRE(A.size() == 0);
	BOOST_REQUIRE(C.data() != 0);
	BOOST_REQUIRE(
	check_data_1d(C, B)
	);
}

BOOST_AUTO_TEST_CASE(test_assignment) {

	BOOST_TEST_MESSAGE("Testing data assignment");

	CALIB::t_data<CALIB::t_real> A(10);
	for (int i = 0; i < A.size(); ++ i) { A(i) = i; }

	CALIB::t_data<CALIB::t_real> B;
	B = A;
	BOOST_REQUIRE(
	(B.data() != A.data()) && check_data_1d(B, A)
	);

	CALIB::t_data<CALIB::t_real> C;
	C = std::move(A);
	BOOST_REQUIRE(A.data() == 0);
	BOOST_REQUIRE(A.size() == 0);
	BOOST_REQUIRE(C.data() != 0);
	BOOST_REQUIRE(
	check_data_1d(C, B)
	);
}

BOOST_AUTO_TEST_CASE(test_resize) {

	BOOST_TEST_MESSAGE("Testing data resize");

	CALIB::t_data<CALIB::t_real> A;

	A.resize(5, 2);
	BOOST_REQUIRE(A.data() != nullptr);
	BOOST_REQUIRE(A.size() == 10);
	BOOST_REQUIRE(A.nrow() == 5);
	BOOST_REQUIRE(A.ncol() == 2);

	A.resize(10);
	BOOST_REQUIRE(A.data() != nullptr);
	BOOST_REQUIRE(A.size() == 10);
	BOOST_REQUIRE(A.nrow() == 10);
	BOOST_REQUIRE(A.ncol() == 1);
}

BOOST_AUTO_TEST_CASE(test_format) {

	BOOST_TEST_MESSAGE("Testing data format");
	//Column-major order is used:
	CALIB::t_data<CALIB::t_long> A(5, 2); CALIB::t_long k = 0;
	for (int j = 0; j < A.ncol(); ++ j)
	for (int i = 0; i < A.nrow(); ++ i) { A(i, j) = k ++; }
	CALIB::t_data<CALIB::t_long> B(5, 2);
	for (int i = 0; i < B.size(); ++ i) { B(i) = i; }
	BOOST_REQUIRE(
	check_data_2d(A, B)
	);
}

BOOST_AUTO_TEST_CASE(test_slice) {

	BOOST_TEST_MESSAGE("Testing data slice");

	CALIB::t_data<CALIB::t_long> A(6, 8);
	for (int j = 0; j < A.ncol(); ++ j)
	for (int i = 0; i < A.nrow(); ++ i) { A(i, j) = -1; }
	CALIB::t_data<CALIB::t_long> B(6, 8);
	B() = -1;
	BOOST_REQUIRE(check_data_2d(A, B));

	for (int j = 0; j < A.ncol(); ++ j) { A(2, j) = 1; }
	B(2, CALIB::t_sect()) = 1;
	BOOST_REQUIRE(check_data_2d(A, B));

	for (int i = 0; i < A.nrow(); ++ i) { A(i, 3) = 2; }
	B(CALIB::t_sect(), 3) = 2;
	BOOST_REQUIRE(check_data_2d(A, B));

	for (int j = 1; j < A.ncol(); j += 3)
	for (int i = 1; i < A.nrow(); i += 2) {
		A(i, j) = 3;
	}
	B(CALIB::t_sect(1, B.nrow(), 2),
	  CALIB::t_sect(1, B.ncol(), 3)
	) = 3;
	BOOST_REQUIRE(
	check_data_2d(A, B)
	);
	//...
}

//...

BOOST_AUTO_TEST_SUITE_END()
