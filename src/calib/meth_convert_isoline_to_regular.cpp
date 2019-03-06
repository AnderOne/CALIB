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

//Отображает поле изолиний на прямоугольную сетку:
//NOTE: The main idea of this procedure was described in [Dritschel, 1997].
//Scalar field is assumed as piece-constant between any two contours.
//Also for its correction average value is considered,
//because it isn't dependent on the time parameter
//into incompressible flows.
c_scal_regular t_meth::run_meth_convert(const c_scal_isoline &FROM, const c_grid_regular &GRID) {

	//Определяем параметры текущего потока:
	t_flow::t_cond cond = GRID->flow()->cond();
	if (FROM->grid()->flow() != GRID->flow()) __ERR_METH("Attempt to call conversion between different flows!");

	//Определяем параметры отображения:
	auto dats = get_data_convert(FROM);
	if (!dats) __ERR_METH("You must first call t_meth::set_data_convert(...)!");
	t_real avgz = dats->avgz;
	t_size fact = dats->fact;

	//Определяем параметры сетки:
	t_real minx = GRID->rect().minx(), miny = GRID->rect().miny();
	t_real maxx = GRID->rect().maxx(), maxy = GRID->rect().maxy();
	t_real lenx = GRID->rect().lenx(), leny = GRID->rect().leny();
	t_long numx = GRID->rect().numx(), numy = GRID->rect().numy();
	const t_grid::t_node &NODE = FROM->grid()->node();
	const t_grid::t_cont &CONT = FROM->grid()->cont();
	const t_grid::t_step &STEP = FROM->grid()->step();
	t_size numk = NODE.size(), numc = CONT.size();
	t_size numz = STEP.size();

	//Генерируем крупную сетку:
	numx = (1 << fact) * (numx - 1) + 1;
	numy = (1 << fact) * (numy - 1) + 1;
	c_grid_regular FINE = t_meth::new_grid_regular(GRID->flow(), numx, numy);
	lenx = FINE->rect().lenx();
	leny = FINE->rect().leny();
	//...
	t_flow::t_geom<t_flow::t_cond::PERIOD2, t_real> POINT(0, numx, 0, numy);
	t_flow::t_geom<t_flow::t_cond::PERIOD2, t_long> INDEX(0, 0, 0, numy);
	CALIB::t_data<t_long> JUMPX(numx, numy);
	CALIB::t_data<t_long> JUMPY(numy);
	JUMPX() = 0;
	JUMPY() = 0;
	//Заполняем прыжки скалярного поля в узлах регулярной сетки:
	for (t_long ix, iy, iz = 0; iz < numz; ++ iz) {
		t_long c0, c1, i0, i1; c1 = (c0 = STEP.head(iz)) + STEP.tail(iz);
		for (t_long ic = c0; ic < c1; ++ ic) {
			i1 = (i0 = CONT.head(ic)) + CONT.tail(ic);
			t_real x2, y2;
			if (CONT.stat(ic)) {
				x2 = NODE.valx(i1 - 1);
				y2 = NODE.valy(i1 - 1);
			}
			//Обрабатываем концы незамкнутых контуров:
			else {
				x2 = NODE.valx(i0);
				y2 = NODE.valy(i0);
				if (NODE.valy(i1 - 1) == miny) {
					JUMPX(((t_long) ((NODE.valx(i1 - 1) - minx) / lenx) + numx) % numx, 0) -= 1;
				}
				if (y2 == miny) {
					JUMPX(((t_long) ((x2 - minx) / lenx) + numx) % numx, 0) += 1;
				}
				if (NODE.valx(i1 - 1) == minx) {
					JUMPY(((t_long) ((NODE.valy(i1 - 1) - miny) / leny) + numy) % numy) += 1;
				}
				if (x2 == minx) {
					JUMPY(((t_long) ((y2 - miny) / leny) + numy) % numy) -= 1;
				}
			}
			//Обходим узлы отдельно взятого контура:
			for (t_long ik = i0; ik < i1; ++ ik) {
				//Переходим к следующему узлу:
				t_real x1 = x2, y1 = y2; x2 = NODE.valx(ik); y2 = NODE.valy(ik);
				//Выполняем нормировку координат узлов:
				t_real cx1 = (x1 - minx) / lenx, cy1 = (y1 - miny) / leny;
				t_real cx2 = (x2 - minx) / lenx, cy2 = (y2 - miny) / leny;
				//Определяем допустимые диапазоны:
				t_long ix1 = (t_long) (cx1), iy1 = (t_long) (cy1);
				t_long ix2 = (t_long) (cx2), iy2 = (t_long) (cy2);
				//Определяем направление обхода:
				t_real dx = POINT.subx(cx2, cx1), dy = POINT.suby(cy2, cy1);
				t_long tx = (dx < 0)? (-1): (+1);
				t_long ty = (dy < 0)? (-1): (+1);
				//Определяем наличие пересечений с левой границей области:
				t_bool cross = (dx * (x2 - x1) < 0);
				//Если пересечения с горизонтальными ребрами отсутствуют:
				if (iy1 == iy2) {
					if (cross) JUMPY((iy1 + numy) % numy) -= tx;
					continue;
				}
				//Определяем число ребер, с которыми пересекается текущий отрезок:
				if (ty < 0) ++ iy2; else ++ iy1;
				t_long my = INDEX.suby(iy2, iy1);
				//assert(my * ty >= 0);
				//Выполняем обход всех таких ребер:
				t_bool event = false;
				for (t_long jy = iy1; jy != (iy1 + my) + ty; jy += ty) {
					//Вычисляем координаты пересечения:
					ix = (t_long) (cx1 + ((jy - cy1) / dy) * dx);
					//Выполняем коррекцию координат (в соответствии с допустимым диапазоном):
					if (cross) {
						if (tx > 0) {
							if ((!event) && (ix >= numx)) {
								//Заполняем прыжки скалярного поля вдоль левой границы:
								JUMPY((jy - (ty > 0) + numy) % numy) -= tx;
								event = true;
								ix -= numx;
							}
						}
						else {
							if ((!event) && (ix < 0)) {
								//Заполняем прыжки скалярного поля вдоль левой границы:
								JUMPY((jy - (ty > 0) + numy) % numy) -= tx;
								event = true;
								ix += numx;
							}
						}
					}
					else {
						if (tx > 0) {
							if (ix < ix1) ix = ix1; if (ix > ix2) ix = ix2;
						}
						else {
							if (ix < ix2) ix = ix2; if (ix > ix1) ix = ix1;
						}
					}
					//Заполняем прыжки скалярного поля вдоль левой границы:
					ix = (ix + numx) % numx;
					iy = (jy + numy) % numy;
					JUMPX(ix, iy) += ty;
				}
				if (!event && cross) {
					//Заполняем прыжки скалярного поля вдоль левой границы:
					iy = (iy2 - (ty <= 0) + numy) % numy;
					JUMPY(iy) -= tx;
				}
			}
		}
	}

	//Суммируем прыжки скалярного поля на сетке:
	t_scal::t_data DATAZ(numx, numy);
	t_real minz = FROM->rect().minz();
	t_real lenz = FROM->rect().lenz();
	DATAZ(0, 0) = minz;	//NOTE: Здесь может стоять любое значение, т.к. все равно используется корректировка!
	for (t_size ix1 = 0, ix = 1; ix < numx; ix1 = ix ++) DATAZ(ix, 0) = DATAZ(ix1, 0) + lenz * JUMPX(ix1, 0);
	for (t_size iy1 = 0, iy = 1; iy < numy; iy1 = iy ++) {
		DATAZ(0, iy) = DATAZ(0, iy1) + lenz * JUMPY(iy1);
		for (t_size ix1 = 0, ix = 1; ix < numx;
		     ix1 = ix ++) {
			DATAZ(ix, iy) = DATAZ(ix1, iy) + lenz * JUMPX(ix1, iy);
		}
	}
	//Проводим корректировку по среднему значению:
	t_real tmpz = blitz::mean(DATAZ());
	DATAZ() += avgz - tmpz;

	//Осредняем значения с мелкой сетки на более грубую (в 2^{fact} раз), используя 9-ти точечный шаблон:
	t_real rx0, rx1, rx2, rx3, ry0, ry1, ry2, ry3;
	t_size tmpx = numx, tmpy = numy; numx = GRID->rect().numx(); numy = GRID->rect().numy();
	while ((fact --) != 0) {
		t_long n, ix, iy, halx = (tmpx + 1) / 2, haly = (tmpy + 1) / 2;
		for (iy = 0; iy < tmpy; ++ iy) {
			if (IS_PERIODX(cond)) {
				rx1 = DATAZ(2 * halx - 3, iy); rx0 = DATAZ(1, iy);
			}
			else {
				rx1 = DATAZ(2, iy); rx0 = DATAZ(2 * halx - 4, iy);
			}
			for (ix = 0; ix < halx - 1; ++ ix) {
				rx3 = DATAZ(2 * ix + 1, iy);
				rx2 = DATAZ(2 * ix, iy);
				DATAZ(ix, iy) = (rx1 + rx2 + rx3) / 3.0;
				rx1 = rx3;
			}
			rx2 = DATAZ(2 * halx - 2, iy);
			DATAZ(halx - 1, iy) = (rx1 + rx2 + rx0) / 3.0;
		}
		for (ix = 0; ix < halx; ++ ix) {
			if (IS_PERIODY(cond)) {
				ry1 = DATAZ(ix, 2 * haly - 3); ry0 = DATAZ(ix, 1);
			}
			else {
				ry1 = DATAZ(ix, 2); ry0 = DATAZ(ix, 2 * haly - 4);
			}
			for (iy = 0; iy < haly - 1; ++ iy) {
				ry3 = DATAZ(ix, 2 * iy + 1);
				ry2 = DATAZ(ix, 2 * iy);
				DATAZ(ix, iy) = (ry1 + ry2 + ry3) / 3.0;
				ry1 = ry3;
			}
			ry2 = DATAZ(ix, 2 * haly - 2);
			DATAZ(ix, haly - 1) = (ry1 + ry2 + ry0) / 3.0;
		}
		tmpx = halx;
		tmpy = haly;
	}
	//Формируем новое поле:
	t_scal::t_data BUFF(numx, numy);
	BUFF(t_sect::all(), t_sect::all()) = DATAZ(
		t_sect(0, numx - 1),
		t_sect(0, numy - 1)
	);
	c_scal_regular SCAL;
	SCAL = t_meth::new_scal_regular(
	GRID, std::move(BUFF)
	);
	//Копируем параметры:
	cpy_data(FROM, SCAL);
	//...
	return SCAL;
}

//Инициализирует конвертор:
t_bool t_meth::set_data_convert(const c_scal &SCAL, t_real avgz, t_size fact) {
	t_hand<t_data_convert> DATA(
		new t_data_convert(avgz, fact)
	);
	set_data_convert(SCAL, DATA);
	return true;
}

}
