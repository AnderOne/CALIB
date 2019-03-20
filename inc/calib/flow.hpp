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

#ifndef __INCLUDE_FLOW_H
#define __INCLUDE_FLOW_H

#include "hand.hpp"
#include "data.hpp"

#include <functional>

namespace CALIB {

/** Дескрипторы потоков жидкости **/

typedef t_weak<const struct t_flow> w_flow;

typedef t_hand<const struct t_flow_plane2d> c_flow_plane2d;
typedef t_hand<const struct t_flow_polygon> c_flow_polygon;
typedef t_hand<const struct t_flow_polar2d> c_flow_polar2d;
typedef t_hand<const struct t_flow> c_flow;

typedef t_hand<struct t_flow_plane2d> h_flow_plane2d;
typedef t_hand<struct t_flow_polygon> h_flow_polygon;
typedef t_hand<struct t_flow_polar2d> h_flow_polar2d;
typedef t_hand<struct t_flow> h_flow;

/** Дескрипторы расчетных сеток **/

typedef t_weak<const struct t_grid> w_grid;

typedef t_hand<const struct t_grid_scatter> c_grid_scatter;
typedef t_hand<const struct t_grid_isoline> c_grid_isoline;
typedef t_hand<const struct t_grid_regular> c_grid_regular;
typedef t_hand<const struct t_grid_triangl> c_grid_triangl;
typedef t_hand<const struct t_grid> c_grid;

typedef t_hand<struct t_grid_scatter> h_grid_scatter;
typedef t_hand<struct t_grid_isoline> h_grid_isoline;
typedef t_hand<struct t_grid_regular> h_grid_regular;
typedef t_hand<struct t_grid_triangl> h_grid_triangl;
typedef t_hand<struct t_grid> h_grid;

/** Дескрипторы векторных полей **/

typedef t_weak<const struct t_vect> w_vect;

typedef t_hand<const struct t_vect_dataset> c_vect_dataset;
typedef t_hand<const struct t_vect_scatter> c_vect_scatter;
typedef t_hand<const struct t_vect_regular> c_vect_regular;
typedef t_hand<const struct t_vect_analith> c_vect_analith;
typedef t_hand<const struct t_vect> c_vect;

typedef t_hand<struct t_vect_dataset> h_vect_dataset;
typedef t_hand<struct t_vect_scatter> h_vect_scatter;
typedef t_hand<struct t_vect_regular> h_vect_regular;
typedef t_hand<struct t_vect_analith> h_vect_analith;
typedef t_hand<struct t_vect> h_vect;

/** Дескрипторы скалярных полей **/

typedef t_weak<const struct t_scal> w_scal;

typedef t_hand<const struct t_scal_dataset> c_scal_dataset;
typedef t_hand<const struct t_scal_scatter> c_scal_scatter;
typedef t_hand<const struct t_scal_isoline> c_scal_isoline;
typedef t_hand<const struct t_scal_regular> c_scal_regular;
typedef t_hand<const struct t_scal_analith> c_scal_analith;
typedef t_hand<const struct t_scal> c_scal;

typedef t_hand<struct t_scal_dataset> h_scal_dataset;
typedef t_hand<struct t_scal_scatter> h_scal_scatter;
typedef t_hand<struct t_scal_isoline> h_scal_isoline;
typedef t_hand<struct t_scal_regular> h_scal_regular;
typedef t_hand<struct t_scal_analith> h_scal_analith;
typedef t_hand<struct t_scal> h_scal;

/** Дескрипторы трассеров **/

typedef t_weak<const struct t_flux> w_flux;

typedef t_hand<const struct t_flux_scatter> c_flux_scatter;
typedef t_hand<const struct t_flux_isoline> c_flux_isoline;
typedef t_hand<const struct t_flux_regular> c_flux_regular;
typedef t_hand<const struct t_flux> c_flux;

typedef t_hand<struct t_flux_scatter> h_flux_scatter;
typedef t_hand<struct t_flux_isoline> h_flux_isoline;
typedef t_hand<struct t_flux_regular> h_flux_regular;
typedef t_hand<struct t_flux> h_flux;

//...

/** Структуры потоков жидкости **/

#define IS_PERIODX(cond) (cond & 0x1)
#define IS_PERIODY(cond) (cond & 0x2)

struct t_flow {

	//Перечислимый тип возможных краевых условий:
	typedef enum { PERIOD0 = 0, PERIODX = 1, PERIODY = 2, PERIOD2 = 3 } t_cond;

	//Тип элемента области (двумерной точки):
	typedef CALIB::t_item<t_real, 2> t_item;

	//Шаблонный класс операций над точками:
	template <t_cond cond, typename TYPE>
	struct t_geom {
		//Разность между точками:
		inline TYPE subx(TYPE lhsx, TYPE rhsx) const { if (IS_PERIODX(cond)) return sub1<0> (lhsx, rhsx); else return sub0<0> (lhsx, rhsx); }
		inline TYPE suby(TYPE lhsy, TYPE rhsy) const { if (IS_PERIODY(cond)) return sub1<1> (lhsy, rhsy); else return sub0<1> (lhsy, rhsy); }
		template <int ind>
		inline TYPE sub1(TYPE lhs, TYPE rhs) const {
			lhs -= rhs; return (lhs < 0)? ((- lhs > hlf[ind])? (lhs + len[ind]): (lhs)): ((lhs > hlf[ind])? (lhs - len[ind]): (lhs));
		}
		template <int ind>
		inline TYPE sub0(TYPE lhs, TYPE rhs) const {
			return lhs - rhs;
		}
		//Контроль точек:
		inline TYPE movx(TYPE valx) const { if (IS_PERIODX(cond)) return mov1<0> (valx); else return mov0<0> (valx); }
		inline TYPE movy(TYPE valy) const { if (IS_PERIODY(cond)) return mov1<1> (valy); else return mov0<1> (valy); }
		template <int ind>
		inline TYPE mov1(TYPE val) const {
			if (val < min[ind]) val += len[ind]; else if (val > max[ind]) val -= len[ind];
			return val;
		}
		template <int ind>
		inline TYPE mov0(TYPE val) const {
			if (val < min[ind]) val = min[ind]; else if (val > max[ind]) val = max[ind];
			return val;
		}
		inline t_geom(TYPE _minx, TYPE _maxx, TYPE _miny, TYPE _maxy):
		              min{_minx, _miny}, max{_maxx, _maxy} {
			len[0] = max[0] - min[0];
			len[1] = max[1] - min[1];
			hlf[0] = len[0] / 2; hlf[1] = len[1] / 2;
		}
		inline t_geom() {}
	private:
		TYPE min[2], max[2], len[2], hlf[2];
	};

	//Прямоугольная рамка:
	struct t_rect {
		inline t_rect(t_real minx, t_real maxx, t_real miny, t_real maxy):
		              MIN{minx, miny}, MAX{maxx, maxy}, LEN{maxx - minx, maxy - miny} {}
		inline t_bool finx() const { return (MIN[0] != -INF) && (MAX[0] != +INF); }
		inline t_bool finy() const { return (MIN[1] != -INF) && (MAX[1] != +INF); }
		inline t_bool fin2() const { return finx() && finy(); }
		inline t_real minx() const { return MIN[0]; }
		inline t_real miny() const { return MIN[1]; }
		inline t_item min2() const { return MIN; }
		inline t_real maxx() const { return MAX[0]; }
		inline t_real maxy() const { return MAX[1]; }
		inline t_item max2() const { return MAX; }
		inline t_real lenx() const { return LEN[0]; }
		inline t_real leny() const { return LEN[1]; }
		inline t_item len2() const { return LEN; }
	private:
		t_item MIN{-INF, -INF};
		t_item MAX{+INF, +INF};
		t_item LEN{+INF, +INF};
	};
	//...
	inline t_real time(t_real _time) { return TIME = _time; }
	inline t_real time() const { return TIME; }
	//...
	friend struct t_meth;
	//...
	virtual ~t_flow();
protected:
	h_flow hand() const;
	t_real TIME = 0;
};

struct t_flow_plane2d:
public t_flow {
	inline const t_flow::t_rect &rect() const { return RECT; }
	inline t_cond cond(t_cond _cond) { return COND = _cond; }
	inline t_cond cond() const { return COND; }
protected:
	inline h_flow_plane2d hand() const {
		return t_flow::hand().get<t_flow_plane2d> ();
	}
	inline t_flow_plane2d(const t_flow::t_rect &_rect,
	                      t_flow::t_cond _cond,
	                      t_real _time):
	                      RECT(_rect), COND(_cond) {
		TIME = _time;
	}
	//...
	friend struct t_meth;
	//...
	t_flow::t_rect RECT;
	t_flow::t_cond COND;
};

/** Структуры расчетных сеток **/

struct t_grid {
	//Массив узлов:
	struct t_node: public t_data<t_real, 2> {
		template <typename ... T> inline auto valx(T ... args) const -> decltype(auto) { return (*this)(args ...)[0]; }
		template <typename ... T> inline auto valy(T ... args) const -> decltype(auto) { return (*this)(args ...)[1]; }
		template <typename ... T> inline auto valx(T ... args) -> decltype(auto) { return (*this)(args ...)[0]; }
		template <typename ... T> inline auto valy(T ... args) -> decltype(auto) { return (*this)(args ...)[1]; }
		inline explicit t_node(t_size size):
		                t_data<t_real, 2> (size) {}
		inline explicit t_node():
		                t_data<t_real, 2> () {}
	};
	//Массив шагов:
	struct t_step: public t_data<t_size, 2> {
		template <typename ... T> inline auto head(T ... args) const -> decltype(auto) { return (*this)(args ...)[0]; }
		template <typename ... T> inline auto tail(T ... args) const -> decltype(auto) { return (*this)(args ...)[1]; }
		template <typename ... T> inline auto head(T ... args) -> decltype(auto) { return (*this)(args ...)[0]; }
		template <typename ... T> inline auto tail(T ... args) -> decltype(auto) { return (*this)(args ...)[1]; }
		inline explicit t_step(t_size size):
		                t_data<t_size, 2> (size) {}
		inline explicit t_step():
		                t_data<t_size, 2> () {}
	};
	//Массив контуров:
	struct t_cont: public t_step {
	private:
		t_data<t_byte, 1> STAT;
	public:
		template <typename ... T> inline auto stat(T ... args) const -> decltype(auto) { return STAT(args ...); }
		template <typename ... T> inline auto stat(T ... args) -> decltype(auto) { return STAT(args ...); }
		inline explicit t_cont(t_size size):
		                t_step(size), STAT(size) {}
		inline explicit t_cont():
		                t_step() {}
	};
	//Рамка:
	struct t_rect {
		inline t_rect(t_real minx, t_real maxx, t_size numx, t_real miny, t_real maxy, t_size numy):
		              LEN{(maxx - minx) / (numx - 1), (maxy - miny) / (numy - 1)},
		              MIN{minx, miny}, MAX{maxx, maxy}, NUM{numx, numy} {}
		inline t_real valx(t_size i) const { return minx() + i * lenx(); }
		inline t_real valy(t_size i) const { return miny() + i * leny(); }
		inline t_real minx() const { return MIN[0]; }
		inline t_real miny() const { return MIN[1]; }
		inline t_real maxx() const { return MAX[0]; }
		inline t_real maxy() const { return MAX[1]; }
		inline t_size numx() const { return NUM[0]; }
		inline t_size numy() const { return NUM[1]; }
		inline t_real lenx() const { return LEN[0]; }
		inline t_real leny() const { return LEN[1]; }
	private:
		CALIB::t_item<t_real, 2> MIN, MAX, LEN;
		CALIB::t_item<t_size, 2> NUM;
	};
	//...
	inline c_flow flow() const { return FLOW; }
	//...
	friend struct t_meth;
	//...
	virtual ~t_grid();
	//...
protected:
	inline t_grid(c_flow _flow): FLOW(_flow) {}
	h_grid hand() const;
	//...
	w_flow FLOW;
};

struct t_grid_regular:
public t_grid {
	inline const t_grid::t_rect &rect() const { return RECT; }
	inline c_flow_plane2d flow() const { return t_grid::flow().get<const t_flow_plane2d> (); }
protected:
	inline t_grid_regular(c_flow_plane2d _flow, const t_grid::t_rect &_rect):
	       t_grid(_flow), RECT(_rect) {}
	//...
	inline h_grid_regular hand() const {
	return t_grid::hand().
	get<t_grid_regular>();
	}
	//...
	friend struct t_meth;
	//...
	t_grid::t_rect RECT;
};

struct t_grid_scatter:
public t_grid {
	inline const t_grid::t_node &node() const { return NODE; }
	inline c_flow_plane2d flow() const { return t_grid::flow().get<const t_flow_plane2d> (); }
protected:
	inline t_grid_scatter(c_flow_plane2d _flow, const t_grid::t_node &_node):
	       t_grid(_flow), NODE(_node) {}
	inline t_grid_scatter(c_flow_plane2d _flow, t_grid::t_node &&_node):
	       t_grid(_flow), NODE(std::move(_node)) {}
	//...
	inline h_grid_scatter hand() const {
	return t_grid::hand().
	get<t_grid_scatter>();
	}
	//...
	friend struct t_meth;
	//...
	t_grid::t_node NODE;
};

struct t_grid_isoline:
public t_grid_scatter {
	inline const t_grid::t_cont &cont() const { return CONT; }
	inline const t_grid::t_step &step() const { return STEP; }
	inline c_flow_plane2d flow() const { return t_grid::flow().get<const t_flow_plane2d> (); }
protected:
	inline t_grid_isoline(c_flow_plane2d _flow, const t_grid::t_node &_node,
	                                            const t_grid::t_cont &_cont,
	                                            const t_grid::t_step &_step):
	       t_grid_scatter(_flow, _node), CONT(_cont), STEP(_step) {}
	inline t_grid_isoline(c_flow_plane2d _flow, t_grid::t_node &&_node,
	                                            t_grid::t_cont &&_cont,
	                                            t_grid::t_step &&_step):
	       t_grid_scatter(_flow, std::move(_node)),
	       CONT(std::move(_cont)),
	       STEP(std::move(_step)) {}
	//...
	inline h_grid_isoline hand() const {
	return t_grid::hand().
	get<t_grid_isoline>();
	}
	//...
	friend struct t_meth;
	//...
	t_grid::t_cont CONT;
	t_grid::t_step STEP;
};

/** Структуры векторных полей **/

struct t_vect {
	//Массив данных:
	struct t_data: public CALIB::t_data<t_real, 2> {
		template <typename ... T>
		inline auto valx(T ... args) const -> decltype(auto) { return (*this)(args ...)[0]; }
		template <typename ... T>
		inline auto valy(T ... args) const -> decltype(auto) { return (*this)(args ...)[1]; }
		template <typename ... T>
		inline auto valx(T ... args) -> decltype(auto) { return (*this)(args ...)[0]; }
		template <typename ... T>
		inline auto valy(T ... args) -> decltype(auto) { return (*this)(args ...)[1]; }
		template <typename ... T>
		inline t_data(T ... args): CALIB::t_data<t_real, 2> (args ...) {}
	};
	//Элемент поля:
	typedef t_data::t_item t_item;
	//Функция поля:
	typedef std::function<
	t_item(t_real, t_real, t_real)
	> t_func;
	//...
	friend struct t_meth;
	//...
	virtual ~t_vect();
protected:
	h_vect hand() const;
};

struct t_vect_dataset:
public t_vect {
	inline c_grid grid() const { return GRID; }
protected:
	inline t_vect_dataset(c_grid _grid):
	                      GRID(_grid) {}
	//...
	friend struct t_meth;
	//...
	c_grid GRID;
};

struct t_vect_scatter:
public t_vect_dataset {
	template <typename ... T> inline auto valx(T ...  args) const { return DATA(args ...)[0]; }
	template <typename ... T> inline auto valy(T ...  args) const { return DATA(args ...)[1]; }
	template <typename ... T> inline auto valz(T ...  args) const { return DATA(args ...); }
	inline const t_vect::t_data &data() const { return DATA; }
	inline c_grid_scatter grid() const { return GRID.get<const t_grid_scatter> (); }
protected:
	inline t_vect_scatter(c_grid _grid, const t_vect::t_data &_data):
	       t_vect_dataset(_grid), DATA(_data) {}
	inline t_vect_scatter(c_grid _grid, t_vect::t_data &&_data):
	       t_vect_dataset(_grid),
	       DATA(std::move(_data)) {}
	//...
	inline h_vect_scatter hand() const {
	return t_vect::hand().
	get<t_vect_scatter>();
	}
	//...
	friend struct t_meth;
	//...
	t_vect::t_data DATA;
};

struct t_vect_regular:
public t_vect_dataset {
	template <typename ... T> inline auto valx(T ...  args) const { return DATA(args ...)[0]; }
	template <typename ... T> inline auto valy(T ...  args) const { return DATA(args ...)[1]; }
	template <typename ... T> inline auto valz(T ...  args) const { return DATA(args ...); }
	inline const t_vect::t_data &data() const { return DATA; }
	inline c_grid_regular grid() const { return GRID.get<const t_grid_regular> (); }
protected:
	inline t_vect_regular(c_grid _grid, const t_vect::t_data &_data):
	       t_vect_dataset(_grid), DATA(_data) {}
	inline t_vect_regular(c_grid _grid, t_vect::t_data &&_data):
	       t_vect_dataset(_grid), DATA(std::move(_data)) {}
	//...
	inline h_vect_regular hand() const {
	return t_vect::hand().
	get<t_vect_regular>();
	}
	//...
	friend struct t_meth;
	//...
	t_vect::t_data DATA;
};

struct t_vect_analith:
public t_vect {
	template <typename ... T> inline auto valx(T ...  args) const { return FUNC(args ...)[0]; }
	template <typename ... T> inline auto valy(T ...  args) const { return FUNC(args ...)[1]; }
	template <typename ... T> inline auto valz(T ...  args) const { return FUNC(args ...); }
	inline const t_vect::t_func &func() const { return FUNC; }
protected:
	inline t_vect_analith(const t_vect::t_func &_func): FUNC(_func) {}
	//...
	inline h_vect_analith hand() const {
	return t_vect::hand().
	get<t_vect_analith>();
	}
	//...
	friend struct t_meth;
	//...
	t_vect::t_func FUNC;
};

/** Структуры скалярных полей **/

struct t_scal {
	//Массив данных:
	struct t_data: public CALIB::t_data<t_real, 1> {
		template <typename ... T>
		inline auto valz(T ... args) const -> decltype(auto) { return (*this)(args ...); }
		template <typename ... T>
		inline auto valz(T ... args) -> decltype(auto) { return (*this)(args ...); }
		template <typename ... T>
		inline t_data(T ... args): CALIB::t_data<t_item, 1> (args ...) {}
	};
	//Уровни поля:
	struct t_rect {
		inline t_rect(t_real minz, t_real maxz, t_size numz): LEN((maxz - minz) / numz),
		       MIN(minz), MAX(maxz), NUM(numz) {}
		inline t_real valz(t_size i) const { return MIN + (i + 0.5) * LEN; }
		inline t_real minz() const { return MIN; }
		inline t_real maxz() const { return MAX; }
		inline t_real lenz() const { return LEN; }
		inline t_size numz() const { return NUM; }
	private:
		t_real MIN, MAX, LEN; t_size NUM;
	};
	//Элемент поля:
	typedef t_data::t_item t_item;
	//Функция поля:
	typedef std::function<
	t_item(t_real, t_real, t_real)
	> t_func;
	//...
	friend struct t_meth;
	//...
	virtual ~t_scal();
protected:
	h_scal hand() const;
};

struct t_scal_dataset:
public t_scal {
	inline c_grid grid() const { return GRID; }
protected:
	inline t_scal_dataset(c_grid _grid):
	                      GRID(_grid) {}
	//...
	friend struct t_meth;
	//...
	c_grid GRID;
};

struct t_scal_isoline:
public t_scal_dataset {
	template <typename ... T> inline auto valz(T ...  args) const { return RECT(args ...); }
	inline const t_scal::t_rect &rect() const { return RECT; }
	inline c_grid_isoline grid() const { return GRID.get<const t_grid_isoline> (); }
protected:
	inline t_scal_isoline(c_grid _grid, const t_scal::t_rect &_rect):
	       t_scal_dataset(_grid), RECT(_rect) {}
	//...
	inline h_scal_isoline hand() const {
	return t_scal::hand().
	get<t_scal_isoline>();
	}
	//...
	friend struct t_meth;
	//...
	t_scal::t_rect RECT;
};

struct t_scal_scatter:
public t_scal_dataset {
	template <typename ... T> inline auto valz(T ...  args) const { return DATA(args ...); }
	inline const t_scal::t_data &data() const { return DATA; }
	inline c_grid_scatter grid() const { return GRID.get<const t_grid_scatter> (); }
protected:
	inline t_scal_scatter(c_grid _grid, const t_scal::t_data &_data):
	       t_scal_dataset(_grid), DATA(_data) {}
	inline t_scal_scatter(c_grid _grid, t_scal::t_data &&_data):
	       t_scal_dataset(_grid),
	       DATA(std::move(_data)) {}
	//...
	inline h_scal_scatter hand() const {
	return t_scal::hand().
	get<t_scal_scatter>();
	}
	//...
	friend struct t_meth;
	//...
	t_scal::t_data DATA;
};

struct t_scal_regular:
public t_scal_dataset {
	template <typename ... T> inline auto valz(T ...  args) const { return DATA(args ...); }
	inline const t_scal::t_data &data() const { return DATA; }
	inline c_grid_regular grid() const { return GRID.get<const t_grid_regular> (); }
protected:
	inline t_scal_regular(c_grid _grid, const t_scal::t_data &_data):
	       t_scal_dataset(_grid), DATA(_data) {}
	inline t_scal_regular(c_grid _grid, t_scal::t_data &&_data):
	       t_scal_dataset(_grid),
	       DATA(std::move(_data)) {}
	//...
	inline h_scal_regular hand() const {
	return t_scal::hand().
	get<t_scal_regular>();
	}
	//...
	friend struct t_meth;
	//...
	t_scal::t_data DATA;
};

struct t_scal_analith:
public t_scal {
	template <typename ... T> inline auto valz(T ...  args) const { return FUNC(args ...); }
	inline const t_scal::t_func &func() const { return FUNC; }
protected:
	inline t_scal_analith(const t_scal::t_func &_func):
	                      FUNC(_func) {}
	//...
	inline h_scal_analith hand() const {
	return t_scal::hand().
	get<t_scal_analith>();
	}
	//...
	friend struct t_meth;
	//...
	t_scal::t_func FUNC;
};

/** Структуры трассеров **/

struct t_flux {
	inline c_scal_dataset scal() const { return SCAL; }
	inline c_grid grid() const { return GRID; }
	inline c_flow flow() const { return GRID != nullptr? GRID->flow(): nullptr; }
	virtual ~t_flux();
protected:
	inline t_flux(c_scal_dataset _scal):
	              GRID(_scal->grid()),
	              SCAL(_scal) {}
	inline t_flux(c_grid_scatter _grid):
	              GRID(_grid) {}
	//...
	friend struct t_meth;
	//...
	c_scal_dataset SCAL;
	c_grid GRID;
	//...
	h_flux hand() const;
};

struct t_flux_scatter:
public t_flux {
	inline c_scal_scatter scal() const { return SCAL.get<const t_scal_scatter> (); }
	inline c_grid_scatter grid() const { return GRID.get<const t_grid_scatter> (); }
protected:
	inline t_flux_scatter(c_scal_scatter _scal):
	       t_flux(_scal) {}
	inline t_flux_scatter(c_grid_scatter _grid):
	       t_flux(_grid) {}
	//...
	inline h_flux_scatter hand() const {
	return t_flux::hand().
	get<t_flux_scatter>();
	}
	//...
	friend struct t_meth;
};

struct t_flux_isoline:
public t_flux {
	inline c_scal_isoline scal() const { return SCAL.get<const t_scal_isoline> (); }
	inline c_grid_isoline grid() const { return GRID.get<const t_grid_isoline> (); }
protected:
	inline t_flux_isoline(c_scal_isoline _scal):
	       t_flux(_scal) {}
	//...
	inline h_flux_isoline hand() const {
	return t_flux::hand().
	get<t_flux_isoline>();
	}
	//...
	friend struct t_meth;
};

struct t_flux_regular:
public t_flux {
	inline c_scal_regular scal() const { return SCAL.get<const t_scal_regular> (); }
	inline c_grid_regular grid() const { return GRID.get<const t_grid_regular> (); }
protected:
	inline t_flux_regular(c_scal_regular _scal):
	       t_flux(_scal) {}
	//...
	inline h_flux_regular hand() const {
	return t_flux::hand().
	get<t_flux_regular>();
	}
	//...
	friend struct t_meth;
};

//...

}

#endif //__INCLUDE_FLOW_H
