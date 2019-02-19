#ifndef __INCLUDE_METH_H
#define __INCLUDE_METH_H

#include "data.hpp"
#include "hand.hpp"
#include "flow.hpp"

#include <exception>
#include <iostream>

#define __ERR_METH(text) { throw std::logic_error(text); }

namespace CALIB {

typedef t_hand<struct t_data_of_scal> h_data_of_scal;
typedef t_hand<struct t_data_of_flow> h_data_of_flow;
typedef t_hand<struct t_data_convert> h_data_convert;
typedef t_hand<struct t_data_correct> h_data_correct;
typedef t_hand<struct t_data_dynamic> h_data_dynamic;
typedef t_hand<struct t_data_inverse> h_data_inverse;

enum t_type_correct { RENODER = 0, SURGERY = 1 };

struct t_meth {

	inline static c_flow_plane2d get_flow_plane2d(const c_flow &_flow) { return _flow.get<const t_flow_plane2d>(); }
	inline static c_grid_scatter get_grid_scatter(const c_grid &_grid) { return _grid.get<const t_grid_scatter>(); }
	inline static c_grid_isoline get_grid_isoline(const c_grid &_grid) { return _grid.get<const t_grid_isoline>(); }
	inline static c_grid_regular get_grid_regular(const c_grid &_grid) { return _grid.get<const t_grid_regular>(); }
	inline static c_vect_scatter get_vect_scatter(const c_vect &_vect) { return _vect.get<const t_vect_scatter>(); }
	inline static c_vect_regular get_vect_regular(const c_vect &_vect) { return _vect.get<const t_vect_regular>(); }
	inline static c_vect_analith get_vect_analith(const c_vect &_vect) { return _vect.get<const t_vect_analith>(); }
	inline static c_scal_scatter get_scal_scatter(const c_scal &_scal) { return _scal.get<const t_scal_scatter>(); }
	inline static c_scal_isoline get_scal_isoline(const c_scal &_scal) { return _scal.get<const t_scal_isoline>(); }
	inline static c_scal_regular get_scal_regular(const c_scal &_scal) { return _scal.get<const t_scal_regular>(); }
	inline static c_scal_analith get_scal_analith(const c_scal &_scal) { return _scal.get<const t_scal_analith>(); }
	inline static c_flux_scatter get_flux_scatter(const c_flux &_flux) { return _flux.get<const t_flux_scatter>(); }
	inline static c_flux_isoline get_flux_isoline(const c_flux &_flux) { return _flux.get<const t_flux_isoline>(); }
	inline static c_flux_regular get_flux_regular(const c_flux &_flux) { return _flux.get<const t_flux_regular>(); }
	inline static h_flow_plane2d get_flow_plane2d(const h_flow &_flow) { return _flow.get<t_flow_plane2d>(); }
	inline static h_grid_scatter get_grid_scatter(const h_grid &_grid) { return _grid.get<t_grid_scatter>(); }
	inline static h_grid_isoline get_grid_isoline(const h_grid &_grid) { return _grid.get<t_grid_isoline>(); }
	inline static h_grid_regular get_grid_regular(const h_grid &_grid) { return _grid.get<t_grid_regular>(); }
	inline static h_vect_scatter get_vect_scatter(const h_vect &_vect) { return _vect.get<t_vect_scatter>(); }
	inline static h_vect_regular get_vect_regular(const h_vect &_vect) { return _vect.get<t_vect_regular>(); }
	inline static h_vect_analith get_vect_analith(const h_vect &_vect) { return _vect.get<t_vect_analith>(); }
	inline static h_scal_scatter get_scal_scatter(const h_scal &_scal) { return _scal.get<t_scal_scatter>(); }
	inline static h_scal_isoline get_scal_isoline(const h_scal &_scal) { return _scal.get<t_scal_isoline>(); }
	inline static h_scal_regular get_scal_regular(const h_scal &_scal) { return _scal.get<t_scal_regular>(); }
	inline static h_scal_analith get_scal_analith(const h_scal &_scal) { return _scal.get<t_scal_analith>(); }
	inline static h_flux_scatter get_flux_scatter(const h_flux &_flux) { return _flux.get<t_flux_scatter>(); }
	inline static h_flux_isoline get_flux_isoline(const h_flux &_flux) { return _flux.get<t_flux_isoline>(); }
	inline static h_flux_regular get_flux_regular(const h_flux &_flux) { return _flux.get<t_flux_regular>(); }
	//Генераторы:
	static c_flow_plane2d new_flow_plane2d(t_real mint, t_real minx, t_real maxx, t_real miny, t_real maxy, t_flow::t_cond cond);
	static c_grid_regular new_grid_regular(const c_flow_plane2d &FLOW, t_size numx, t_size numy);
	static c_grid_isoline new_grid_isoline(const c_flow_plane2d &FLOW, const t_grid::t_node &NODE,
	                                                                   const t_grid::t_cont &CONT, const t_grid::t_step &STEP);
	static c_grid_scatter new_grid_scatter(const c_flow_plane2d &FLOW, const t_grid::t_node &NODE);
	static c_grid_isoline new_grid_isoline(const c_flow_plane2d &FLOW, t_grid::t_node &&NODE,
	                                                                   t_grid::t_cont &&CONT, t_grid::t_step &&STEP);
	static c_grid_scatter new_grid_scatter(const c_flow_plane2d &FLOW, t_grid::t_node &&NODE);
	static c_vect_regular new_vect_regular(const c_grid_regular &GRID, const t_vect::t_data &DATA);
	static c_vect_scatter new_vect_scatter(const c_grid_scatter &GRID, const t_vect::t_data &DATA);
	static c_vect_regular new_vect_regular(const c_grid_regular &GRID, t_vect::t_data &&DATA);
	static c_vect_scatter new_vect_scatter(const c_grid_scatter &GRID, t_vect::t_data &&DATA);
	static c_vect_analith new_vect_analith(const t_vect::t_func &FUNC);
	static c_scal_regular new_scal_regular(const c_grid_regular &GRID, const t_scal::t_data &DATA);
	static c_scal_scatter new_scal_scatter(const c_grid_scatter &GRID, const t_scal::t_data &DATA);
	static c_scal_regular new_scal_regular(const c_grid_regular &GRID, t_scal::t_data &&DATA);
	static c_scal_scatter new_scal_scatter(const c_grid_scatter &GRID, t_scal::t_data &&DATA);
	static c_scal_isoline new_scal_isoline(const c_grid_isoline &GRID, t_real minz, t_real maxz);
	static c_scal_analith new_scal_analith(const t_scal::t_func &FUNC);
	static c_flux_regular new_flux_regular(const c_scal_regular &SCAL, t_name name);
	static c_flux_isoline new_flux_isoline(const c_scal_isoline &SCAL, t_name name);
	static c_flux_scatter new_flux_scatter(const c_scal_scatter &SCAL, t_name name);
	static c_flux_scatter new_flux_scatter(const c_grid_scatter &GRID, t_name name);
	static c_flux_regular new_flux_regular(const c_scal_regular &SCAL);
	static c_flux_isoline new_flux_isoline(const c_scal_isoline &SCAL);
	static c_flux_scatter new_flux_scatter(const c_scal_scatter &SCAL);
	static c_flux_scatter new_flux_scatter(const c_grid_scatter &GRID);
	//Конверторы скалярных полей:
	static c_scal_isoline run_meth_convert(const c_scal_regular &FROM, const t_scal::t_rect &RECT);
	static c_scal_regular run_meth_convert(const c_scal_isoline &FROM, const c_grid_regular &GRID);
	static c_scal_regular run_meth_convert(const c_scal_analith &FROM, const c_grid_regular &GRID);
	static c_scal_scatter run_meth_convert(const c_scal_analith &FROM, const c_grid_scatter &GRID);
	static c_scal_scatter run_meth_convert(const c_scal_regular &FROM, const c_grid_scatter &GRID);
	//Конверторы векторных полей:
	static c_vect_regular run_meth_convert(const c_vect_analith &FROM, const c_grid_regular &GRID);
	static c_vect_scatter run_meth_convert(const c_vect_analith &FROM, const c_grid_scatter &GRID);
	static c_vect_scatter run_meth_convert(const c_vect_regular &FROM, const c_grid_scatter &GRID);
	//Корректоры скалярных полей:
	static c_scal_isoline run_meth_correct(const c_scal_isoline &FROM, t_type_correct type);
	//Обработчики потоков:
	static t_bool run_meth_correct(const c_flow_plane2d &FLOW, t_type_correct type);
	static t_bool run_meth_dynamic(const c_flow_plane2d &FLOW);
	static t_bool run_meth_inverse(const c_flow_plane2d &FLOW);
	//Инициализирует инвертор:
	static t_bool set_data_inverse(const c_flow_plane2d &FLOW, const c_scal_isoline &PVORT,
	                                                           const c_scal_regular &DEPTH, t_real f0 = 0, t_real bt = 0);
	static t_bool set_data_inverse(const c_flow_plane2d &FLOW, const c_scal_isoline &PVORT,
	                                                           const c_grid_regular &GRID, t_real f0 = 0, t_real bt = 0);
	static t_bool set_data_inverse(const c_flow_plane2d &FLOW, const c_vect_regular &VECT);
	static t_bool set_data_inverse(const c_flow_plane2d &FLOW, const c_vect_analith &VECT);
	//Инициализирует адвектор:
	static t_bool set_data_dynamic(const c_flow &FLOW, t_real lent);
	//Инициализирует конвертор:
	static t_bool set_data_convert(const c_scal &SCAL, t_real avgz, t_size fact);
	//Инициализирует корректор:
	static t_bool set_data_correct(const c_scal &SCAL, t_real minl, t_real maxl);

	//Вычисляет якобиан для двух скалярных полей, заданных на регулярной сетке:
	static c_scal_regular run_meth_regular_jacobian(const c_scal_regular &LHSF, const c_scal_regular &RHSF);
	//Вычисляет поле скорости по полю функции тока, заданного на регулярной сетке:
	static c_vect_regular run_meth_regular_velocity(const c_scal_regular &SCAL);
	//Вычисляет градиент скалярного поля, заданного на регулярной сетке:
	static c_vect_regular run_meth_regular_gradient(const c_scal_regular &SCAL);

	//Геттеры:
	static c_flux_regular get_flux_regular(const c_flow &FLOW, t_name name) {
		return get_flux_regular(get_flux(FLOW, name.tolower()));
	}
	static c_flux_isoline get_flux_isoline(const c_flow &FLOW, t_name name) {
		return get_flux_isoline(get_flux(FLOW, name.tolower()));
	}
	static c_flux_scatter get_flux_scatter(const c_flow &FLOW, t_name name) {
		return get_flux_scatter(get_flux(FLOW, name.tolower()));
	}
	static c_flux get_flux(const c_flow &FLOW, t_name name);
	static t_real get_test(const c_flow &FLOW, t_name name);

	virtual ~t_meth() {}
protected:
	template <typename T> inline static auto get_hand(const T &_hand) {
		return _hand->hand();
	}
	static const t_dict<t_name, c_flux> &get_flux(const c_flow &FLOW);
	//...
	static void cpy_data(const t_flow *FROM, const t_flow *FLOW);
	static void cpy_data(const t_grid *FROM, const t_grid *GRID);
	static void cpy_data(const t_vect *FROM, const t_vect *VECT);
	static void cpy_data(const t_scal *FROM, const t_scal *SCAL);
	//...
	static void del_data(const t_flow *FLOW);
	static void del_data(const t_grid *GRID);
	static void del_data(const t_vect *VECT);
	static void del_data(const t_scal *SCAL);
	static void del_data(const t_flux *FLUX);
	//...
	static h_data_convert
	set_data_convert(const t_scal *SCAL, h_data_convert DATA);
	static h_data_correct
	set_data_correct(const t_scal *SCAL, h_data_correct DATA);
	static h_data_dynamic
	set_data_dynamic(const t_flow *FLOW, h_data_dynamic DATA);
	static h_data_inverse
	set_data_inverse(const t_flow *FLOW, h_data_inverse DATA);
	//...
	static h_data_of_scal
	new_data_of_scal(const t_scal *SCAL);
	static h_data_of_flow
	new_data_of_flow(const t_flow *FLOW);
	//...
	static h_data_convert
	get_data_convert(const t_scal *SCAL);
	static h_data_correct
	get_data_correct(const t_scal *SCAL);
	static h_data_of_scal
	get_data_of_scal(const t_scal *SCAL);
	static h_data_dynamic
	get_data_dynamic(const t_flow *FLOW);
	static h_data_inverse
	get_data_inverse(const t_flow *FLOW);
	static h_data_of_flow
	get_data_of_flow(const t_flow *FLOW);
	//...
	friend struct t_flow;
	friend struct t_grid;
	friend struct t_vect;
	friend struct t_scal;
	friend struct t_flux;
	//...
	inline t_meth() {}
};

//...
struct t_data_convert {
	inline t_data_convert(t_real _avgz, t_size _fact): avgz(_avgz), fact(_fact) {}
	t_real avgz; t_size fact;
};
//...
struct t_data_correct {
	inline t_data_correct(t_real _minl, t_real _maxl): minl(_minl), maxl(_maxl) {}
	/*virtual c_scal_isoline run_renoder(const c_scal_isoline &_scal) = 0;
	virtual c_scal_isoline run_surgery(const c_scal_isoline &_scal) = 0;
	virtual ~t_data_correct();*/
	t_real minl; t_real maxl;
};
//...
struct t_data_dynamic {
	inline t_data_dynamic(t_real _lent): lent(_lent) {}
	t_real lent;
};
//...
struct t_data_inverse {
	virtual c_scal_regular run(const c_scal_regular &_scal) = 0;
	virtual c_vect_scatter run(const c_grid_scatter &_grid) = 0;
	virtual t_bool run(t_real dt) = 0;
	virtual t_bool run() = 0;
	virtual ~t_data_inverse() {}
};

struct t_data_of_scal {
	t_hand<t_data_convert> CONVERT;
	t_hand<t_data_correct> CORRECT;
};

struct t_data_of_flow {
	t_hand<t_data_dynamic> DYNAMIC;
	t_hand<t_data_inverse> INVERSE;
};

//...

}

#endif //__INCLUDE_METH_H
