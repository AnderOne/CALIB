#include <calib/calib.hpp>
#include "file.hpp"

#include <iostream>

#define _USE_MATH_DEFINES
#include <tr1/cmath>
#include <cmath>

using std::tr1::cyl_bessel_j;

using namespace CALIB;

//[J.H.G.M. van Geffen, P.A. Davies, 1999-2000]:
static const t_real f0 = + 5.0, beta = 0.3;	//Северное полушарие;
//static const t_real f0 = - 5.0, beta = 0.3;	//Южное полушарие;
//static const t_real f0 = 0.0, beta = 0.3;	//Экватор;

t_real depth(t_real x, t_real y, t_real t) {
	t_real h = 0.0, r = sqrt(x * x + y * y), h0 = 0.4, R = 1.0; if (r < R) h = h0 * cos(r / R * (M_PI / 2.0));
	return 1.0 - h;
}

t_real pvort(t_real x, t_real y, t_real t) {

	const t_real a = 0.5, k = 2.4048 / a, G = 4, j1 = cyl_bessel_j(1, k * a);
	t_real dx = x - 3.0, dy = y + 3.0, r = sqrt(dx * dx + dy * dy);
	t_real f = (r < a)? (
	(((k * a) * G) / (2 * M_PI * (a * a) * j1)) * cyl_bessel_j(0, k * r)
	): (0);
	return (f + (f0 + beta * y)) / depth(x, y, t);
}

int main() {

	//Формируем расчетную область и покрывающие ее регулярные сетки:
	c_flow_plane2d FLOW = t_meth::new_flow_plane2d(0, - 5 * M_PI, + 5 * M_PI, - 5 * M_PI, + 5 * M_PI, t_flow::t_cond::PERIODX);
	c_grid_regular FINE = t_meth::new_grid_regular(FLOW, 1024, 1024);
	c_grid_regular GRID = t_meth::new_grid_regular(FLOW, 256, 256);

	//Формируем сеточное поле на основе аналитической функции:
	c_scal_regular SCAL = t_meth::run_meth_convert(t_meth::new_scal_analith(pvort), FINE);
	t_real minz = min(SCAL->valz()), maxz = max(SCAL->valz());

	//Формируем карту изолиний на основе сеточного поля:
	c_scal_isoline ISOL = t_meth::run_meth_convert(SCAL, t_scal::t_rect(minz, maxz, 256));

	//Устанавливаем параметры конверсии:
	t_meth::set_data_convert(ISOL, mean(SCAL->valz()), 0);
	//Устанавливаем параметры коррекции:
	t_meth::set_data_correct(ISOL, 0.0025,  0.07);

	//Устанавливаем параметры инверсии:
	t_meth::set_data_inverse(
	FLOW, ISOL, t_meth::run_meth_convert(t_meth::new_scal_analith(depth), GRID), f0, beta
	);
	//Устанавливаем параметры адвекции:
	t_meth::set_data_dynamic(FLOW, 0.025);

	//Выполняем перераспределение узлов:
	ISOL = t_meth::run_meth_correct(ISOL, RENODER);
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
			t_meth::run_meth_correct(FLOW, SURGERY);
		}
		t_meth::run_meth_correct(FLOW, RENODER);
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
