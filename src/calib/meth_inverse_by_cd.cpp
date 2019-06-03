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

#define PI (3.14159265358979323846)

namespace CALIB {

//Абстрактный класс инвертора по замкнутым интегралам:
struct t_meth_integrc  {

	virtual t_vect::t_data get(const c_scal_isoline &VORT, const t_grid::t_node &MARK) = 0;
	virtual ~t_meth_integrc() {}
};

template <t_flow::t_cond cond>
struct t_meth_integrc_temp:
public t_meth_integrc {

	virtual t_vect::t_data get(const c_scal_isoline &VORT, const t_grid::t_node &MARK) {

		const t_grid::t_cont &CONT = VORT->grid()->cont();
		const t_grid::t_node &NODE = VORT->grid()->node();
		t_real lenx = VORT->grid()->flow()->rect().lenx();
		t_real leny = VORT->grid()->flow()->rect().leny();

		t_vect::t_data DATA(MARK.size());

		for (t_size k = 0; k < MARK.size(); ++ k) {

			t_real x  = MARK.valx(k), y  = MARK.valy(k), dz = VORT->rect().lenz();
			t_real fx = 0, fy = 0;

			for (t_size c = 0; c < CONT.size(); ++ c) {

				t_size i0  = CONT.head(c), nn = CONT.tail(c), in = i0 + nn;
				t_real x1  = NODE.valx(in - 2), y1  = NODE.valy(in - 2);
				t_real x2  = NODE.valx(in - 1), y2  = NODE.valy(in - 1);
				t_real dx1 = GEOM.subx(x2, x1), dy1 = GEOM.suby(y2, y1);
				t_real sx = 0, sy = 0;

				for (t_size i = i0; i < in; ++ i) {

					t_real x3 = NODE.valx(i), y3 = NODE.valy(i);
					t_real dx2 = GEOM.subx(x3, x2);
					t_real dy2 = GEOM.suby(y3, y2);
					t_real dx = GEOM.subx(x, x2);
					t_real dy = GEOM.suby(y, y2);
					t_real fz;
					if (std::abs(dx) + std::abs(dy) > 1.e-14) {
						//NOTE: Green's function!
						if (cond == t_flow::t_cond::PERIODX) {
							fz = log(std::abs(cosh((2 * PI) * dy / lenx) - cos((2 * PI) * dx / lenx)));
						}
						else
						if (cond == t_flow::t_cond::PERIODY) {
							fz = log(std::abs(cosh((2 * PI) * dx / leny) - cos((2 * PI) * dy / leny)));
						}
						else {
							fz = log(dx * dx + dy * dy);
						}
						sx += fz * (dx1 + dx2);
						sy += fz * (dy1 + dy2);
					}

					dx1 = dx2; dy1 = dy2;
					x1 = x2; y1 = y2;
					x2 = x3; y2 = y3;
				}

				fx += sx / (2 * nn);
				fy += sy / (2 * nn);
				
			}

			DATA.valx(k) = - dz * fx;
			DATA.valy(k) = - dz * fy;
		}
		//...
		return DATA;
	}

	inline t_meth_integrc_temp(const t_flow::t_rect &RECT): GEOM(
	RECT.minx(), RECT.maxx(), RECT.miny(), RECT.maxy()
	) {}

	virtual ~t_meth_integrc_temp() {}
private:
	t_flow::t_geom<cond, t_real> GEOM;
};

struct t_meth_inverse_of_cd:
public t_data_inverse,
public t_meth {

	inline explicit t_meth_inverse_of_cd(const c_flow_plane2d &_flow): FLOW(_flow) {

		t_flow::t_cond cond = FLOW->cond();

		if (cond == t_flow::t_cond::PERIOD0) METH = new t_meth_integrc_temp<t_flow::t_cond::PERIOD0> (FLOW->rect());
		if (cond == t_flow::t_cond::PERIODX) METH = new t_meth_integrc_temp<t_flow::t_cond::PERIODX> (FLOW->rect());
		if (cond == t_flow::t_cond::PERIODY) METH = new t_meth_integrc_temp<t_flow::t_cond::PERIODY> (FLOW->rect());
		//...
	}

	//NOTE: Для бесконечных потоков трассеры на регулярных сетках не предусмотрены!
	c_scal_regular run(const c_scal_regular &_scal) { return nullptr; }

	c_vect_scatter run(const c_grid_scatter &_grid) {
		c_scal_isoline VORT = t_meth::get_flux_isoline(FLOW, "RVORT")->scal();
		return t_meth::new_vect_scatter(
			_grid, METH->get(VORT, _grid->node())
		);
	}

	t_bool run(t_real dt) { return true; }

	t_bool run() { return true; }
protected:
	t_weak<const t_flow_plane2d> FLOW;
	t_hand<t_meth_integrc> METH;
};

t_bool t_meth::set_data_inverse(const c_flow_plane2d &FLOW, const c_scal_isoline &RVORT) {

	if ((RVORT->grid()->flow() != FLOW)) {
		__ERR_METH("Incompatible parameters of inversion! (objects belonging different flows)");
	}
	//Определяем параметры потока:
	t_flow::t_cond cond = FLOW->cond();
	//...
	if (FLOW->rect().fin2()) {
		__ERR_METH("This version not support Contour Dynamics for a finite flow!");
	}
	//Проверяем идентификаторы:
	if (get_flux(FLOW, "RVORT") != nullptr) {
		__ERR_METH("Special name of the flux \"RVORT\" is busy!");
	}
	t_meth::new_flux_isoline(
		RVORT, "RVORT"
	);
	//...
	t_hand<t_data_inverse> METH = new t_meth_inverse_of_cd(FLOW);
	set_data_inverse(
		FLOW, METH
	);
	//...
	return true;
}

//...

}
