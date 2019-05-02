#include <calib/calib.hpp>
#include "file.hpp"

#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>

using namespace CALIB;

int main() {

	//Формируем расчетную область:
	//c_flow_plane2d FLOW = t_meth::new_flow_plane2d(0, - M_PI, + M_PI, - INF, + INF, t_flow::t_cond::PERIODX);
	//c_flow_plane2d FLOW = t_meth::new_flow_plane2d(0, - INF, + INF, - M_PI, + M_PI, t_flow::t_cond::PERIODY);
	c_flow_plane2d FLOW = t_meth::new_flow_plane2d(0, - INF, + INF, - INF, + INF, t_flow::t_cond::PERIOD0);

	//Формируем поле трассера:
	t_size numk = 1000, numc = 1, nums = 1;
	t_grid::t_node NODE(numk * numc); t_grid::t_cont CONT(numc); t_grid::t_step STEP(nums);
	t_real xx[] = {0};//-1.0, -1.0,  1.0, 1.0};
	t_real yy[] = {0};//-1.0,  1.0, -1.0, 1.0};
	for (t_size k = 0; k < numc; ++ k) {
		t_real h = 2 * M_PI / numk;
		t_real p = 0;
		for (t_size i = 0; i < numk; ++ i) {
			NODE(k * numk + i) = t_flow::t_item{xx[k] + cos(p), yy[k] + sin(p) / 2};
			p += h;
		}
	}
	for (t_size i = 0; i < numc; ++ i) {
		CONT.head(i) = i * numk; CONT.tail(i) = numk; CONT.stat(i) = 1;
	}
	STEP.head(0) = 0; STEP.tail(0) = numc;
	c_grid_isoline GRID = t_meth::new_grid_isoline(
		FLOW, std::move(NODE), std::move(CONT), std::move(STEP)
	);
	c_scal_isoline ISOL = t_meth::new_scal_isoline(
		GRID, -1.0, +1.0
	);
	//Устанавливаем параметры коррекции:
	t_meth::set_data_correct(ISOL, 0.0001, 0.01);

	//Устанавливаем параметры инверсии:
	t_meth::set_data_inverse(FLOW, ISOL);
	//Устанавливаем параметры адвекции:
	t_meth::set_data_dynamic(FLOW, 0.1);

	//Выполняем перераспределение узлов:
	ISOL = t_meth::run_meth_correct(ISOL, t_meth::RENODER);

	f_scal_isoline fisol(fopen("isol.dat", "wb"));
	fisol.write(ISOL);

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
		fisol.write(ISOL);
	}

	return 0;
}
