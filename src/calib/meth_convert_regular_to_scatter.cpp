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
struct t_meth_convert_regular_to_scatter:
public t_meth {

	//Выполняет билинейную интерполяцию от регулярной сетки к точкам:
	//NOTE: Bi-linear interpolation is used.
	template <t_flow::t_cond cond>
	inline static t_hand<const T_OUT> run_convert(const t_hand<const T_INP> &FROM, const c_grid_scatter &GRID) {

		t_real minx = FROM->grid()->rect().minx(), miny = FROM->grid()->rect().miny();
		t_real maxx = FROM->grid()->rect().maxx(), maxy = FROM->grid()->rect().maxy();
		t_real lenx = FROM->grid()->rect().lenx(), leny = FROM->grid()->rect().leny();
		t_size numx = FROM->grid()->rect().numx(), numy = FROM->grid()->rect().numy();
		t_size numk = GRID->node().size();

		const typename T_INP::t_data &BUFF = FROM->data(); const t_grid::t_node &NODE = GRID->node();
		typename T_OUT::t_data DATA(numk);

		for (t_long ix1, iy1, ix2, iy2, ik = 0; ik < numk; ++ ik) {
			t_real cx1, cy1, wx1, wy1, wx2, wy2;
			//Определяем индексы текущей ячейки:
			cx1 = (NODE.valx(ik) - minx) / lenx; cy1 = (NODE.valy(ik) - miny) / leny;
			ix2 = ((ix1 = floor(cx1)) + 1);
			iy2 = ((iy1 = floor(cy1)) + 1);
			//Учитываем граничные условия:
			if (!IS_PERIODX(cond)) {
				if (ix1 >= numx) ix1 = numx - 1; else if (ix1 < 0) ix1 = 0;
				if (ix2 >= numx) ix2 = numx - 1; else if (ix2 < 0) ix2 = 0;
			}
			else {
				ix1 = (ix1 + numx) % numx; ix2 = (ix2 + numx) % numx;
			}
			if (!IS_PERIODY(cond)) {
				if (iy1 >= numy) iy1 = numy - 1; else if (iy1 < 0) iy1 = 0;
				if (iy2 >= numy) iy2 = numy - 1; else if (iy2 < 0) iy2 = 0;
			}
			else {
				iy1 = (iy1 + numy) % numy; iy2 = (iy2 + numy) % numy;
			}
			//...
			wx1 = 1.0 - (wx2 = cx1 - ix1); wy1 = 1.0 - (wy2 = cy1 - iy1);
			DATA(ik) =
			wy1 * (wx1 * BUFF(ix1, iy1) +
			       wx2 * BUFF(ix2, iy1)
			) +
			wy2 * (wx1 * BUFF(ix1, iy2) +
			       wx2 * BUFF(ix2, iy2)
			);
		}
		//...
		return _new(
			GRID, std::move(DATA)
		);
	}

	inline static t_hand<const T_OUT> run(const t_hand<const T_INP> &FROM, const c_grid_scatter &GRID) {

		t_flow::t_cond cond = GRID->flow()->cond();
		if (FROM->grid()->flow() != GRID->flow()) {
			__ERR_METH("Attempt to call conversion between different flows!");
		}
		if (cond == t_flow::t_cond::PERIOD0)
			return run_convert<t_flow::t_cond::PERIOD0> (FROM, GRID);
		if (cond == t_flow::t_cond::PERIODX)
			return run_convert<t_flow::t_cond::PERIODX> (FROM, GRID);
		if (cond == t_flow::t_cond::PERIODY)
			return run_convert<t_flow::t_cond::PERIODY> (FROM, GRID);
		if (cond == t_flow::t_cond::PERIOD2)
			return run_convert<t_flow::t_cond::PERIOD2> (FROM, GRID);
		//...
		return nullptr;
	}
};

c_vect_scatter t_meth::run_meth_convert(const c_vect_regular &FROM, const c_grid_scatter &GRID) {
	c_vect_scatter VECT =
	t_meth_convert_regular_to_scatter<
	t_vect_regular, t_vect_scatter
	>::run(FROM, GRID);
	t_meth::cpy_data(FROM, VECT);
	return VECT;
}

c_scal_scatter t_meth::run_meth_convert(const c_scal_regular &FROM, const c_grid_scatter &GRID) {
	c_scal_scatter SCAL =
	t_meth_convert_regular_to_scatter<
	t_scal_regular, t_scal_scatter
	>::run(FROM, GRID);
	t_meth::cpy_data(FROM, SCAL);
	return SCAL;
}

}
