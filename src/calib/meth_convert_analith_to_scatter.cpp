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

inline static c_vect_scatter _new(const c_grid_scatter &GRID, t_vect::t_data &&DATA) { return t_meth::new_vect_scatter(GRID, std::move(DATA)); }

inline static c_scal_scatter _new(const c_grid_scatter &GRID, t_scal::t_data &&DATA) { return t_meth::new_scal_scatter(GRID, std::move(DATA)); }

template <typename T_INP, typename T_OUT>
struct t_meth_convert_analith_to_scatter:
public t_meth {

	inline static t_hand<const T_OUT> run(const t_hand<const T_INP> &FROM, const c_grid_scatter &GRID) {

		t_real time = GRID->flow()->time();
		t_size numk = GRID->node().size(); typename T_OUT::t_data DATA(numk);
		//...
		for (int ik = 0; ik < numk; ++ ik) {
			DATA(ik) = FROM->valz(
				GRID->node().valx(ik), GRID->node().valy(ik), time
			);
		}
		return _new(
			GRID, std::move(DATA)
		);
	}
};

c_vect_scatter t_meth::run_meth_convert(const c_vect_analith &FROM, const c_grid_scatter &GRID) {
	c_vect_scatter VECT =
	t_meth_convert_analith_to_scatter<
	t_vect_analith, t_vect_scatter
	>::run(FROM, GRID);
	t_meth::cpy_data(FROM, VECT);
	return VECT;
}

c_scal_scatter t_meth::run_meth_convert(const c_scal_analith &FROM, const c_grid_scatter &GRID) {
	c_scal_scatter SCAL =
	t_meth_convert_analith_to_scatter<
	t_scal_analith, t_scal_scatter
	>::run(FROM, GRID);
	t_meth::cpy_data(FROM, SCAL);
	return SCAL;
}

}
