#include <calib/calib.hpp>
#include "file.hpp"

#include <functional>
#include <exception>
#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>

using namespace CALIB;

t_real func1(t_real x, t_real y, t_real t) { t_real r = sqrt(x * x + y * y); return sin(y) * cos(x) + ((r < 1.0)? (1.0 - r): (0)); }
t_real func2(t_real x, t_real y, t_real t) { return sin(x - y) * cos(x + y); }
t_real func3(t_real x, t_real y, t_real t) { return x * y; }

t_flow::t_cond get_cond() {

	std::cout << "Select type of boundary condition:\n(0) -- PERIOD0;\n(1) -- PERIODX;\n(2) -- PERIODY;\n(3) -- PERIOD2;\n";
	static int t;
	if (!(std::cin >> t)) throw std::invalid_argument("Invalid input!");
	if (t == 0) return t_flow::t_cond::PERIOD0;
	if (t == 1) return t_flow::t_cond::PERIODX;
	if (t == 2) return t_flow::t_cond::PERIODY;
	if (t == 3) return t_flow::t_cond::PERIOD2;
	throw std::out_of_range("Unknown value!");
}

t_scal::t_func get_func() {

	std::cout << "Select test number [from 1 to 3]:\n";
	static int t;
	if (!(std::cin >> t)) throw std::invalid_argument("Invalid input!");
	if (t == 1) return func1;
	if (t == 2) return func2;
	if (t == 3) return func3;
	throw std::out_of_range("Unknown test!");
}

int main() {

	t_flow::t_cond cond = get_cond();
	t_scal::t_func func = get_func();

	//Формируем расчетную область и покрывающие ее регулярные сетки:
	c_flow_plane2d FLOW = t_meth::new_flow_plane2d(0, - M_PI, + M_PI, - M_PI, + M_PI, cond);

	//Формируем сеточное поле на основе аналитической функции:
	c_scal_regular SCAL = t_meth::run_meth_convert(
		t_meth::new_scal_analith(func), t_meth::new_grid_regular(FLOW, 1024, 1024)
	);

	//Формируем карту изолиний на основе сеточного поля:
	c_scal_isoline ISOL = t_meth::run_meth_convert(
		SCAL, t_scal::t_rect(min(SCAL->valz()), max(SCAL->valz()), 32)
	);

	//Задаем параметры коррекции:
	t_meth::set_data_correct(ISOL, 0.001, 0.5);

	//Выполняем перераспределение узлов на контуре:
	ISOL = t_meth::run_meth_correct(ISOL, RENODER);

	f_scal_isoline fisol(fopen("isol.dat", "wb"));
	fisol.write(ISOL);

	return 0;
}
