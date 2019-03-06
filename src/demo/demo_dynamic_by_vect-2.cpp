#include <calib/calib.hpp>
#include "file.hpp"

#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>

using namespace CALIB;

//[Doswell, 1984]
t_vect::t_item veloc(t_real x, t_real y, t_real t) {

	t_real p = atan2(y, x), r = sqrt(x * x + y * y);
#if false
	t_real vt = 0, d = cosh(r); if (fabs(d) > 1.e-14) { vt = tanh(r) / (d * d); vt = vt * (1 - sin(t)); }	//Нестационарное поле скоростей!
#else
	t_real vt = 0, d = cosh(r); if (fabs(d) > 1.e-14) { vt = tanh(r) / (d * d); }	//Стационарное поле скоростей!
#endif
	t_real u = - sin(p) * vt;
	t_real v = cos(p) * vt;

	return t_vect::t_item(u, v);
}

t_scal::t_item field(t_real x, t_real y, t_real t) {
	return tanh(y);
}

int main() {

	//Формируем расчетную область:
	c_flow_plane2d FLOW = t_meth::new_flow_plane2d(0, - M_PI, + M_PI, - M_PI, + M_PI, t_flow::t_cond::PERIOD0);

	//Формируем поле трассера:
	c_grid_regular GRID = t_meth::new_grid_regular(FLOW, 256, 256);
	c_scal_regular SCAL = t_meth::run_meth_convert(t_meth::new_scal_analith(field), GRID);

	//Устанавливаем параметры инверсии:
	t_meth::set_data_inverse(FLOW, t_meth::new_vect_analith(veloc));
	//Устанавливаем параметры адвекции:
	t_meth::set_data_dynamic(FLOW, 0.01);

	f_scal_regular fscal(fopen("scal.dat", "wb"));
	fscal.write(SCAL);

	//Запускаем адвекцию:
	t_meth::new_flux_regular(SCAL, "SCAL");
	for (int i = 1; i <= 2000; ++ i) {
		std::cout << "time-step: " << i << std::endl;
		t_meth::run_meth_dynamic(FLOW);
		if (i % 20 != 0) continue;
		SCAL = t_meth::get_flux_regular(
			FLOW, "SCAL"
		)->scal();
		fscal.write(SCAL);
	}

	return 0;
}
