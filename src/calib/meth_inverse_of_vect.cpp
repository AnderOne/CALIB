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

/** Инвертор на основе аналитического поля скорости **/

struct t_data_inverse_of_vect_analith:
public t_data_inverse,
public t_meth {

	inline explicit t_data_inverse_of_vect_analith(const c_vect_analith &_vect): VECT(_vect) {}

	c_scal_regular run(const c_scal_regular &_scal) {

		c_vect_regular DIFF = t_meth::run_meth_regular_gradient(_scal);
		c_grid_regular GRID = _scal->grid();
		c_vect_regular TEMP = t_meth::run_meth_convert(VECT, GRID);
		t_scal::t_data DATA(
		         GRID->rect().numx(), GRID->rect().numy()
		);
		DATA() = TEMP->valx() * DIFF->valx() +
		         TEMP->valy() * DIFF->valy();
		return t_meth::new_scal_regular(
		GRID, std::move(DATA)
		);
	}

	c_vect_scatter run(const c_grid_scatter &_grid) {
		return
		t_meth::run_meth_convert(VECT, _grid);
	}

	t_bool run(t_real dt) {
		return true;
	}
	t_bool run() {
		return true;
	}
private:
	c_vect_analith VECT;
};

t_bool t_meth::set_data_inverse(const c_flow_plane2d &FLOW, const c_vect_analith &VECT) {

	if (get_data_inverse(FLOW) != nullptr) {
		__ERR_METH("Can't overwrite existing invertor!");
		return false;
	}
	set_data_inverse(
	FLOW, new t_data_inverse_of_vect_analith(VECT)
	);
	return true;
}

/** Инвертор на основе сеточного поля скорости **/

struct t_data_inverse_of_vect_regular:
public t_data_inverse,
public t_meth {

	inline explicit t_data_inverse_of_vect_regular(const c_vect_regular &_vect):
	                                               LAST(_vect), NEXT(_vect), VECT(_vect) {}

	c_scal_regular run(const c_scal_regular &_scal) {

		c_vect_regular DIFF = t_meth::run_meth_regular_gradient(_scal);
		c_grid_regular GRID = _scal->grid();
		t_scal::t_data DATA(GRID->rect().numx(), GRID->rect().numy());
		c_vect_regular TEMP;
		if (GRID != VECT->grid())
		         TEMP = t_meth::run_meth_convert(VECT, GRID);
		else {
		         TEMP = VECT;
		}
		DATA() = TEMP->valx() * DIFF->valx() +
		         TEMP->valy() * DIFF->valy();
		return t_meth::new_scal_regular(
		GRID, std::move(DATA)
		);
	}

	c_vect_scatter run(const c_grid_scatter &_grid) {
		return
		t_meth::run_meth_convert(VECT, _grid);
	}

	t_bool run(t_real dt) {
		if (dt == 1) { VECT = NEXT; return true; }
		if (dt == 0) { VECT = LAST; return true; }
		if (LAST != NEXT) {
			t_vect::t_data DATA(LAST->data().nrow(), LAST->data().ncol());
			DATA() = LAST->valz() + (NEXT->valz() - LAST->valz()) * dt;
			VECT = t_meth::new_vect_regular(
			LAST->grid(), std::move(DATA)
			);
			return true;
		}
		VECT = LAST;
		return true;
	}

	t_bool run() {
		LAST = NEXT;
		return true;
	}
private:
	friend struct t_meth;
	c_vect_regular LAST;
	c_vect_regular NEXT;
	c_vect_regular VECT;
};

t_bool t_meth::set_data_inverse(const c_flow_plane2d &FLOW, const c_vect_regular &VECT) {

	auto TEMP = get_data_inverse(FLOW);
	auto DATA = TEMP.get<t_data_inverse_of_vect_regular> ();
	if (TEMP == nullptr) {
		set_data_inverse(FLOW, new t_data_inverse_of_vect_regular(VECT));
		return true;
	}
	if (DATA == nullptr) {
		__ERR_METH("Can't overwrite existing invertor!");
		return false;
	}
	if (DATA->LAST->grid() != VECT->grid()) {
		__ERR_METH("The fields belong different grids!");
		return false;
	}
	DATA->NEXT = VECT;
	return true;
}

//...

}
