#include <calib/calib.hpp>
#include "file.hpp"

#include <functional>
#include <exception>
#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>

using namespace CALIB;

//Ellipse:
t_real func1(t_real x, t_real y, t_real t) {

	t_real r = 2.0 * sqrt(0.25 * x * x + y * y), f = (2 * M_PI) * ((r < 1.0)? sqrt(1 - r): (0)); return f;
}

//Диполь:
t_real func2(t_real x, t_real y, t_real t) {

	t_real x1 = x - 0.5, y1 = y - 0.5, r1 = sqrt(x1 * x1 + y1 * y1);
	t_real x2 = x + 0.5, y2 = y + 0.5, r2 = sqrt(x2 * x2 + y2 * y2);
	return - (((r1 < 1.5)? (1.0 - exp(1.5 - r1)): 0) + ((r2 < 1.5)? (1.0 - exp(1.5 - r2)): 0));
}

//Модон:
t_real func3(t_real x, t_real y, t_real t) {

	t_real x1 = x - 0.5, y1 = y - 0.5, r1 = sqrt(x1 * x1 + y1 * y1);
	t_real x2 = x + 0.5, y2 = y + 0.5, r2 = sqrt(x2 * x2 + y2 * y2);
	return - (((r1 < 1.5)? (1.0 - exp(1.5 - r1)): 0) + ((r2 < 1.5)? (exp(1.5 - r2) - 1.0): 0));
}

//...
t_real func4(t_real x, t_real y, t_real t) {

	return sin(x) * cos(y) + sin(x + y) * cos(x - y);
}

//Jet:
t_real func5(t_real x, t_real y, t_real t) {

	t_real y1 = y - 0.1 * (sin(2 * x) - sin(3 * x));
	if (fabs(y1) < 1.0) {
		t_long s0 = ((y1 < 0)? (-1): (+1));
		y1 = 10 * (0.5 - fabs(fabs(y1) - 0.5));
		y1 = s0 * y1;
	}
	else {
		y1 = 0;
	}
	return y1;
}

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

	std::cout << "Select test number [from 1 to 5]:\n";
	static int t;
	if (!(std::cin >> t)) throw std::invalid_argument("Invalid input!");
	if (t == 1) return func1;
	if (t == 2) return func2;
	if (t == 3) return func3;
	if (t == 4) return func4;
	if (t == 5) return func5;
	throw std::out_of_range("Unknown test!");
}

int main() {

	t_flow::t_cond cond = get_cond();
	t_scal::t_func func = get_func();

	//Формируем расчетную область и покрывающие ее регулярные сетки:
	c_flow_plane2d FLOW = t_meth::new_flow_plane2d(0, - M_PI, + M_PI, - M_PI, + M_PI, cond);
	c_grid_regular FINE = t_meth::new_grid_regular(FLOW, 1024, 1024);
	c_grid_regular GRID = t_meth::new_grid_regular(FLOW, 256, 256);

	//Формируем сеточное поле на основе аналитической функции:
	c_scal_regular SCAL = t_meth::run_meth_convert(t_meth::new_scal_analith(func), FINE);
	t_real minz = min(SCAL->valz());
	t_real maxz = max(SCAL->valz());

	//Формируем карту изолиний на основе сеточного поля:
	c_scal_isoline ISOL = t_meth::run_meth_convert(SCAL, t_scal::t_rect(minz, maxz, 32));

	//Устанавливаем параметры конверсии:
	t_meth::set_data_convert(ISOL, mean(SCAL->valz()), 0);
	//Устанавливаем параметры коррекции:
	t_meth::set_data_correct(ISOL, 0.0001, 0.01);

	//Устанавливаем параметры инверсии:
	t_meth::set_data_inverse(FLOW, ISOL, GRID);
	//Устанавливаем параметры адвекции:
	t_meth::set_data_dynamic(FLOW, 0.01);

	//Выполняем перераспределение узлов:
	ISOL = t_meth::run_meth_correct(ISOL, t_meth::RENODER);
	//Отображаем поле контуров на сетку:
	SCAL = t_meth::run_meth_convert(ISOL, FINE);

	f_scal_isoline fisol(fopen("isol.dat", "wb")); fisol.write(ISOL);
	f_scal_regular fscal(fopen("scal.dat", "wb")); fscal.write(SCAL);

	//Запускаем динамику:
	t_meth::new_flux_isoline(ISOL, "ISOL");
	for (int i = 1; i <= 2000; ++ i) {
		std::cout << "time-step: " << i << std::endl;
		t_meth::run_meth_dynamic(FLOW);
		if (i % 20 == 0) {
			t_meth::run_meth_correct(FLOW, t_meth::SURGERY);
		}
		t_meth::run_meth_correct(FLOW, t_meth::RENODER);
		if (i % 20 != 0) continue;
		ISOL = t_meth::get_flux_isoline(
			FLOW, "ISOL"
		)->scal();
		SCAL = t_meth::run_meth_convert(
			ISOL, FINE
		);
		fisol.write(ISOL);
		fscal.write(SCAL);
	}

	return 0;
}
