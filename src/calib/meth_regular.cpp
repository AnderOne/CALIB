/**
 * Copyright (c) 2015-2017, 2019 Andrey Baranov <armath123@gmail.com>
 *
 * This file is part of CALIB (Contour Advection Library).
 *
 * CALIB is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 2 of the License, or (at your
 * option) any later version.
 *
 * CALIB is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
 * License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with CALIB;
 * if not, see <http://www.gnu.org/licenses/>
**/

#include <calib/calib.hpp>

namespace CALIB {

//Обобщенный метод взятия производной:
//по X:
static t_scal::t_data get_data_regular_differ_x(const t_scal::t_data &DATA, const t_grid::t_rect &RECT, t_flow::t_cond cond) {
	//...
	t_size numx = RECT.numx(), numy = RECT.numy(); t_real lenx = RECT.lenx(), leny = RECT.leny(); t_scal::t_data DIFF(numx, numy);
	//Вычисляем значения во внутренних узлах:
	DIFF(t_sect(1, numx - 2), t_sect::all()) = (DATA(t_sect(2, numx - 1), t_sect::all()) -
	                                            DATA(t_sect(0, numx - 3), t_sect::all())
	                                           ) / (2.0 * lenx);
	//Вычисляем значения в граничных узлах:
	if (IS_PERIODX(cond)) {
		DIFF(numx - 1, t_sect::all()) = (DATA(1, t_sect::all()) - DATA(numx - 2, t_sect::all())) / (2.0 * lenx);
		DIFF(0, t_sect::all()) = DIFF(numx - 1, t_sect::all());
	}
	else {
		DIFF(numx - 1, t_sect::all()) = (DATA(numx - 1, t_sect::all()) -
		                                 DATA(numx - 2, t_sect::all())
		                                ) / (lenx);
		DIFF(0, t_sect::all()) = (DATA(1, t_sect::all()) -
		                          DATA(0, t_sect::all())
		                         ) / (lenx);
	}
	//...
	return std::move(DIFF);
}
//по Y:
static t_scal::t_data get_data_regular_differ_y(const t_scal::t_data &DATA, const t_grid::t_rect &RECT, t_flow::t_cond cond) {
	//...
	t_size numx = RECT.numx(), numy = RECT.numy(); t_real lenx = RECT.lenx(), leny = RECT.leny(); t_scal::t_data DIFF(numx, numy);
	//Вычисляем значения во внутренних узлах:
	DIFF(t_sect::all(), t_sect(1, numy - 2)) = (DATA(t_sect::all(), t_sect(2, numy - 1)) -
	                                            DATA(t_sect::all(), t_sect(0, numy - 3))
	                                           ) / (2.0 * leny);
	//Вычисляем значения в граничных узлах:
	if (IS_PERIODY(cond)) {
		DIFF(t_sect::all(), numy - 1) = (DATA(t_sect::all(), 1) - DATA(t_sect::all(), numy - 2)) / (2.0 * leny);
		DIFF(t_sect::all(), 0) = DIFF(t_sect::all(), numy - 1);
	}
	else {
		DIFF(t_sect::all(), numy - 1) = (DATA(t_sect::all(), numy - 1) -
		                                 DATA(t_sect::all(), numy - 2)
		                                ) / (leny);
		DIFF(t_sect::all(), 0) = (DATA(t_sect::all(), 1) -
		                          DATA(t_sect::all(), 0)
		                         ) / (leny);
	}
	//...
	return std::move(DIFF);
}

static t_vect::t_data get_data_regular_gradient(const t_scal::t_data &DATA, const t_grid::t_rect &RECT, t_flow::t_cond cond) {
	//...
	t_size numx = RECT.numx(), numy = RECT.numy(); t_vect::t_data VECT(numx, numy);
	VECT()[0] = get_data_regular_differ_x(DATA, RECT, cond)();
	VECT()[1] = get_data_regular_differ_y(DATA, RECT, cond)();
	//...
	return std::move(VECT);
}

static t_vect::t_data get_data_regular_velocity(const t_scal::t_data &DATA, const t_grid::t_rect &RECT, t_flow::t_cond cond) {
	//...
	t_size numx = RECT.numx(), numy = RECT.numy(); t_vect::t_data VECT(numx, numy);
	VECT()[0] = - get_data_regular_differ_y(DATA, RECT, cond)();
	VECT()[1] = get_data_regular_differ_x(DATA, RECT, cond)();
	//...
	return std::move(VECT);
}

static t_scal::t_data get_data_regular_jacobian(const t_scal::t_data &LHSF, const t_scal::t_data &RHSF,
                                                const t_grid::t_rect &RECT, t_flow::t_cond cond) {

	t_size numx = RECT.numx(), numy = RECT.numy();
	t_real lenx = RECT.lenx(), leny = RECT.leny();

	t_scal::t_data JAC(numx, numy), JAC1(numx, numy), JAC2(numx, numy);
	t_scal::t_data DPX, DPY, DZX, DZY;

	//Вычисляем производные сеточных полей (с учетом краевых условий):
	DPX = get_data_regular_differ_x(RHSF, RECT, cond);
	DPY = get_data_regular_differ_y(RHSF, RECT, cond);
	DZX = get_data_regular_differ_x(LHSF, RECT, cond);
	DZY = get_data_regular_differ_y(LHSF, RECT, cond);
	//Вычисляем якобиан по схеме Аракавы:
	JAC1() = DPX() * DZY() - DPY() * DZX();
	DPX() *= LHSF();
	DPY() *= LHSF();
	//...
	DPX = get_data_regular_differ_y(DPX, RECT, cond);
	DPY = get_data_regular_differ_x(DPY, RECT, cond);
	//...
	JAC2() = DPX() - DPY();
	DZX() *= RHSF();
	DZY() *= RHSF();
	//...
	DZX = get_data_regular_differ_y(DZX, RECT, cond);
	DZY = get_data_regular_differ_x(DZY, RECT, cond);
	//...
	//JAC3() = DZX() - DZY(); JAC() = (JAC1() + JAC2() + JAC3()) / 3.0;
	JAC() = (
		(JAC1() + JAC2()) + (DZX() - DZY())
	) / 3.0;
	//...
	return std::move(JAC);
}

//Вычисляет градиент скалярного поля, заданного на регулярной сетке:
c_vect_regular t_meth::run_meth_regular_gradient(const c_scal_regular &SCAL) {

	//Определяем параметры текущего потока:
	t_flow::t_cond cond = SCAL->grid()->flow()->cond();
	//...
	t_vect::t_data VECT = get_data_regular_gradient(
	SCAL->data(), SCAL->grid()->rect(), cond
	);
	//...
	return t_meth::new_vect_regular(
	SCAL->grid(), std::move(VECT)
	);
}

//Вычисляет поле скорости по полю функции тока, заданного на регулярной сетке:
c_vect_regular t_meth::run_meth_regular_velocity(const c_scal_regular &SCAL) {

	//Определяем параметры текущего потока:
	t_flow::t_cond cond = SCAL->grid()->flow()->cond();
	//...
	t_vect::t_data VECT = get_data_regular_velocity(
	SCAL->data(), SCAL->grid()->rect(), cond
	);
	//...
	return t_meth::new_vect_regular(
	SCAL->grid(), std::move(VECT)
	);
}

//Вычисляет якобиан для двух скалярных полей, заданных на регулярной сетке:
c_scal_regular t_meth::run_meth_regular_jacobian(const c_scal_regular &LHSF,
                                                 const c_scal_regular &RHSF) {

	//Определяем параметры текущего потока:
	if (LHSF->grid() != RHSF->grid()) {
	__ERR_METH("Attempt to call jacobian between different grids!");
	}
	t_flow::t_cond cond = LHSF->grid()->flow()->cond();
	c_grid_regular GRID = LHSF->grid();
	//...
	return t_meth::new_scal_regular(
	GRID, get_data_regular_jacobian(
	LHSF->data(), RHSF->data(),
	GRID->rect(), cond
	)
	);
}

//...

}
