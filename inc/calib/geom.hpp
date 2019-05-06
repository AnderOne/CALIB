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

#ifndef __INCLUDE_GEOM_H
#define __INCLUDE_GEOM_H

namespace CALIB {

//Перечислимый тип возможных краевых условий:
enum t_cond { PERIOD0 = 0, PERIODX = 1, PERIODY = 2, PERIOD2 = 3 };

#define IS_PERIODX(cond) (cond & 0x1)
#define IS_PERIODY(cond) (cond & 0x2)

//Шаблонный класс операций над координатами:
template <t_cond cond, typename TYPE>
struct t_geom {

	t_geom(TYPE _minx, TYPE _maxx, TYPE _miny, TYPE _maxy):
	min{_minx, _miny}, max{_maxx, _maxy} {
		hlf[0] = (len[0] = max[0] - min[0]) / 2;
		hlf[1] = (len[1] = max[1] - min[1]) / 2;
	}
	t_geom() {}

	//Разность между координатами:
	TYPE subx(TYPE lhsx, TYPE rhsx) const { return (IS_PERIODX(cond))? (sub1<0> (lhsx, rhsx)): (sub0<0> (lhsx, rhsx)); }
	TYPE suby(TYPE lhsy, TYPE rhsy) const { return (IS_PERIODY(cond))? (sub1<1> (lhsy, rhsy)): (sub0<1> (lhsy, rhsy)); }
	//Контроль координат:
	TYPE movx(TYPE valx) const { return (IS_PERIODX(cond))? (mov1<0> (valx)): (mov0<0> (valx)); }
	TYPE movy(TYPE valy) const { return (IS_PERIODY(cond))? (mov1<1> (valy)): (mov0<1> (valy)); }

private:
	template <int ind>
	TYPE sub1(TYPE lhs, TYPE rhs) const {
		return ((lhs = lhs - rhs) < 0)? ((- lhs > hlf[ind])? (lhs + len[ind]): (lhs)):
		                                ((lhs > hlf[ind])? (lhs - len[ind]): (lhs));
	}
	template <int ind>
	TYPE sub0(TYPE lhs, TYPE rhs) const {
		return lhs - rhs;
	}
	template <int ind>
	TYPE mov1(TYPE val) const {
		if (val < min[ind]) return val + len[ind];
		if (val > max[ind]) return val - len[ind];
		return val;
	}
	template <int ind>
	TYPE mov0(TYPE val) const {
		if (val < min[ind]) return min[ind];
		if (val > max[ind]) return max[ind];
		return val;
	}

	TYPE min[2], max[2];
	TYPE len[2], hlf[2];
};

}

#endif //__INCLUDE_GEOM_H
