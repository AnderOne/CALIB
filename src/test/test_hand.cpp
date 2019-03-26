#define BOOST_TEST_MODULE test_hand

#include <boost/test/unit_test.hpp>

#include <calib/hand.hpp>

BOOST_AUTO_TEST_SUITE(suite_of_hand_and_weak_tests)

struct t_test_data_1 {

	virtual ~t_test_data_1() { end = val; }

	inline t_test_data_1(int a): val(a) {}

	virtual int type() const { return 1; }

	static int end;

	int val;
};

struct t_test_data_2: public t_test_data_1 {

	virtual ~t_test_data_2() { end = val; }

	inline t_test_data_2(int a):
	       t_test_data_1(a) {}

	virtual int type() const { return 2; }

	static int end;
};

int t_test_data_1::end = 0;
int t_test_data_2::end = 0;


BOOST_AUTO_TEST_CASE(test_polymorphism) {

	BOOST_TEST_MESSAGE("Testing polymorphic capabilities");
	CALIB::t_hand<t_test_data_1> A;
	CALIB::t_hand<t_test_data_2> B;
	//RTTI:
	A = new t_test_data_1(123);
	BOOST_REQUIRE(A->type() == 1);
	B = A.get<t_test_data_2>();
	BOOST_REQUIRE(A != 0);
	BOOST_REQUIRE(B == 0);
	//...
	A = new t_test_data_2(123);
	BOOST_REQUIRE(A->type() == 2);
	B = A.get<t_test_data_2>();
	BOOST_REQUIRE(A != 0);
	BOOST_REQUIRE(B == A);
	//...
	B.off();
	int &e1 = t_test_data_1::end;
	int &e2 = t_test_data_2::end;
	e1 = e2 = 0;
	A.off();
	BOOST_REQUIRE(e1 && e2);
}

BOOST_AUTO_TEST_CASE(test_qualifiers) {

	BOOST_TEST_MESSAGE("Testing qualifiers compatibility");
	CALIB::t_hand<const int> A = new int(1);
	CALIB::t_hand<int> B = new int(2);
	BOOST_REQUIRE(A != B);
	//B = A;	//TEST: Compilation error!
	A = B;
	BOOST_REQUIRE(A == B);
	//A.del();	//TEST: Compilation error!
	B.del();
	BOOST_REQUIRE(A == 0);
	BOOST_REQUIRE(B == 0);
}

BOOST_AUTO_TEST_CASE(test_constructor) {

	BOOST_TEST_MESSAGE("Testing constructors of t_hand");
	CALIB::t_hand<int> H1(new int(123));
	BOOST_REQUIRE(H1 != nullptr);
	BOOST_REQUIRE(*H1 == 123);
	CALIB::t_hand<int> H2(H1);
	BOOST_REQUIRE(H2 == H1);
	CALIB::t_hand<int> H3;
	BOOST_REQUIRE(H3 == nullptr);

	BOOST_TEST_MESSAGE("Testing constructors of t_weak");
	CALIB::t_weak<int> W1(H1);
	BOOST_REQUIRE(W1 == H1);
	CALIB::t_weak<int> W2(W1);
	BOOST_REQUIRE(W2 == W1);
	CALIB::t_weak<int> W3;
	BOOST_REQUIRE(W3 == nullptr);
}

BOOST_AUTO_TEST_CASE(test_assignment) {

	BOOST_TEST_MESSAGE("Testing assignment of t_hand");
	CALIB::t_hand<int> A, B;
	BOOST_REQUIRE(A == nullptr);
	BOOST_REQUIRE(B == nullptr);
	A = new int(1);
	BOOST_REQUIRE(A != nullptr);
	BOOST_REQUIRE(B == nullptr);
	const int *p = A;
	BOOST_REQUIRE(p == A);
	B = A;
	BOOST_REQUIRE(B == A);
	BOOST_REQUIRE(B == p);
	*A = 2;
	BOOST_REQUIRE(*A == 2);
	BOOST_REQUIRE(*B == 2);
	B = new int(3);
	BOOST_REQUIRE(*A == 2);
	BOOST_REQUIRE(*B == 3);
	A = nullptr;
	BOOST_REQUIRE(A == 0);
	BOOST_REQUIRE(B != 0);
	B = A;
	BOOST_REQUIRE(A == 0);
	BOOST_REQUIRE(B == 0);

	BOOST_TEST_MESSAGE("Testing assignment of t_weak");
	A = new int(4);
	CALIB::t_weak<int> W;
	W = A;
	BOOST_REQUIRE(W == A);
	B = W;
	BOOST_REQUIRE(B == A);
	W = B;
	BOOST_REQUIRE(W == A);
	B = nullptr;
	BOOST_REQUIRE(W != 0);
	BOOST_REQUIRE(W == A);
	A = nullptr;
	BOOST_REQUIRE(W == 0);
}

BOOST_AUTO_TEST_CASE(test_unlink) {

	BOOST_TEST_MESSAGE("Testing to unlink references");
	CALIB::t_hand<int> A, B;
	BOOST_REQUIRE(A == 0);
	BOOST_REQUIRE(B == 0);
	A = new int(1);
	BOOST_REQUIRE(A != 0);
	BOOST_REQUIRE(B == 0);
	B = A;
	BOOST_REQUIRE(B == A);
	B.off();
	BOOST_REQUIRE(A != 0);
	BOOST_REQUIRE(B == 0);
	B = A;
	BOOST_REQUIRE(B == A);
	A.del();
	BOOST_REQUIRE(A == 0);
	BOOST_REQUIRE(B == 0);

	A = new int(2);
	CALIB::t_weak<int> W;
	W = A;
	BOOST_REQUIRE(W == A);
	B = W;
	BOOST_REQUIRE(B == A);
	W = B;
	BOOST_REQUIRE(W == A);
	B.off();
	BOOST_REQUIRE(W != 0);
	BOOST_REQUIRE(W == A);
	A.off();
	BOOST_REQUIRE(W == 0);
}

BOOST_AUTO_TEST_CASE(test_move) {

	BOOST_TEST_MESSAGE("Testing to move resource");
	CALIB::t_hand<t_test_data_1> A, B;
	CALIB::t_weak<t_test_data_1> W;
	A = new t_test_data_1(1);
	B = A;
	W = B;
	BOOST_REQUIRE(W == B);
	BOOST_REQUIRE(B == A);
	A = new t_test_data_1(2);
	BOOST_REQUIRE(W == B);
	BOOST_REQUIRE(B != A);
	t_test_data_1::end = 0;
	B.mov(A);
	BOOST_REQUIRE(
	t_test_data_1::end == 1
	);
	BOOST_REQUIRE(W == B);
	BOOST_REQUIRE(B == A);
}

BOOST_AUTO_TEST_CASE(test_pick) {

	BOOST_TEST_MESSAGE("Testing to pick resource");
	CALIB::t_hand<t_test_data_1> A, B;
	CALIB::t_weak<t_test_data_1> W;
	A = new t_test_data_1(1);
	B = A;
	W = B;
	BOOST_REQUIRE(W == B);
	BOOST_REQUIRE(B == A);
	t_test_data_1::end = 0;
	A = A.out();
	BOOST_REQUIRE(
	t_test_data_1::end == 0
	);
	BOOST_REQUIRE(W == 0);
	BOOST_REQUIRE(B == 0);
	BOOST_REQUIRE(A != 0);
}

BOOST_AUTO_TEST_SUITE_END()
