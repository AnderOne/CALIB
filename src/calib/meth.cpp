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

#include <algorithm>
#include <map>
#include <set>

namespace CALIB {

namespace HIDE {

	static std::map<const t_flow *, std::set<const t_grid *> > GRID_OF_FLOW;
	static std::map<const t_grid *, std::set<const t_vect *> > VECT_OF_GRID;
	static std::map<const t_grid *, std::set<const t_scal *> > SCAL_OF_GRID;
	//...
	static std::map<const t_flow *, std::map<t_name, c_flux> > FLUX_OF_FLOW;
	//...
	static std::map<const t_scal *, h_data_of_scal> DATA_OF_SCAL;
	static std::map<const t_flow *, h_data_of_flow> DATA_OF_FLOW;
	//...
	static std::map<const t_flow *, t_weak<t_flow> > FLOW;
	static std::map<const t_grid *, t_weak<t_grid> > GRID;
	static std::map<const t_vect *, t_weak<t_vect> > VECT;
	static std::map<const t_scal *, t_weak<t_scal> > SCAL;
	static std::map<const t_flux *, t_weak<t_flux> > FLUX;
	//...
}

/** Генераторы потоков жидкости **/

c_flow_plane2d t_meth::new_flow_plane2d(t_real mint, t_real minx, t_real maxx, t_real miny, t_real maxy, t_flow::t_cond cond) {

	t_flow::t_rect RECT = t_flow::t_rect{minx, maxx, miny, maxy};
	if ((!RECT.finx() && IS_PERIODX(cond)) || (!RECT.finy() && IS_PERIODY(cond))) {
		__ERR_METH("Incorrect the boundary conditions to infinity flow!");
	}
	h_flow_plane2d F = new t_flow_plane2d(RECT, cond, mint);
	HIDE::FLOW[F] = F;
	//...
	return F;
}

/** Генераторы расчетных сеток **/

c_grid_regular t_meth::new_grid_regular(const c_flow_plane2d &FLOW, t_size numx, t_size numy) {

	const t_flow::t_rect &R = FLOW->rect();
	if (!R.fin2()) {
		__ERR_METH("Can't create a regular grid based on infinity flow!");
	}
	h_grid_regular G = new t_grid_regular(
	FLOW, t_grid::t_rect(R.minx(), R.maxx(), numx, R.miny(), R.maxy(), numy)
	);
	HIDE::GRID_OF_FLOW[FLOW].insert(G);
	HIDE::GRID[G] = G;
	//...
	return G;
}

c_grid_isoline t_meth::new_grid_isoline(const c_flow_plane2d &FLOW, const t_grid::t_node &NODE,
                                                                    const t_grid::t_cont &CONT, const t_grid::t_step &STEP) {
	h_grid_isoline G = new t_grid_isoline(FLOW, NODE, CONT, STEP);
	//...
	HIDE::GRID_OF_FLOW[FLOW].insert(G);
	HIDE::GRID[G] = G;
	//...
	return G;
}

c_grid_isoline t_meth::new_grid_isoline(const c_flow_plane2d &FLOW, t_grid::t_node &&NODE,
                                                                    t_grid::t_cont &&CONT, t_grid::t_step &&STEP) {
	h_grid_isoline G = new t_grid_isoline(
	FLOW, std::move(NODE), std::move(CONT), std::move(STEP)
	);
	//...
	HIDE::GRID_OF_FLOW[FLOW].insert(G);
	HIDE::GRID[G] = G;
	//...
	return G;
}

c_grid_scatter t_meth::new_grid_scatter(const c_flow_plane2d &FLOW, const t_grid::t_node &NODE) {
	//...
	h_grid_scatter G = new t_grid_scatter(FLOW, NODE);
	//...
	HIDE::GRID_OF_FLOW[FLOW].insert(G);
	HIDE::GRID[G] = G;
	//...
	return G;
}

c_grid_scatter t_meth::new_grid_scatter(const c_flow_plane2d &FLOW, t_grid::t_node &&NODE) {
	//...
	h_grid_scatter G = new t_grid_scatter(
	FLOW, std::move(NODE)
	);
	//...
	HIDE::GRID_OF_FLOW[FLOW].insert(G);
	HIDE::GRID[G] = G;
	//...
	return G;
}

/** Генераторы векторных полей **/

c_vect_regular t_meth::new_vect_regular(const c_grid_regular &GRID, const t_vect::t_data &DATA) {
	//...
	h_vect_regular V = new t_vect_regular(GRID, DATA);
	//...
	HIDE::VECT_OF_GRID[GRID].insert(V);
	HIDE::VECT[V] = V;
	//...
	return V;
}

c_vect_scatter t_meth::new_vect_scatter(const c_grid_scatter &GRID, const t_vect::t_data &DATA) {
	//...
	h_vect_scatter V = new t_vect_scatter(GRID, DATA);
	//...
	HIDE::VECT_OF_GRID[GRID].insert(V);
	HIDE::VECT[V] = V;
	//...
	return V;
}

c_vect_regular t_meth::new_vect_regular(const c_grid_regular &GRID, t_vect::t_data &&DATA) {
	//...
	h_vect_regular V = new t_vect_regular(
	GRID, std::move(DATA)
	);
	//...
	HIDE::VECT_OF_GRID[GRID].insert(V);
	HIDE::VECT[V] = V;
	//...
	return V;
}

c_vect_scatter t_meth::new_vect_scatter(const c_grid_scatter &GRID, t_vect::t_data &&DATA) {
	//...
	h_vect_scatter V = new t_vect_scatter(
	GRID, std::move(DATA)
	);
	//...
	HIDE::VECT_OF_GRID[GRID].insert(V);
	HIDE::VECT[V] = V;
	//...
	return V;
}

c_vect_analith t_meth::new_vect_analith(const t_vect::t_func &FUNC) {
	//...
	h_vect_analith V =
	new t_vect_analith(FUNC);
	//...
	HIDE::VECT[V] = V;
	//...
	return V;
}

/** Генераторы скалярных полей **/

c_scal_regular t_meth::new_scal_regular(const c_grid_regular &GRID, const t_scal::t_data &DATA) {
	//...
	h_scal_regular S = new t_scal_regular(GRID, DATA);
	//...
	HIDE::SCAL_OF_GRID[GRID].insert(S);
	HIDE::SCAL[S] = S;
	//...
	return S;
}

c_scal_scatter t_meth::new_scal_scatter(const c_grid_scatter &GRID, const t_scal::t_data &DATA) {
	//...
	h_scal_scatter S = new t_scal_scatter(GRID, DATA);
	//...
	HIDE::SCAL_OF_GRID[GRID].insert(S);
	HIDE::SCAL[S] = S;
	//...
	return S;
}

c_scal_regular t_meth::new_scal_regular(const c_grid_regular &GRID, t_scal::t_data &&DATA) {
	//...
	h_scal_regular S = new t_scal_regular(
	GRID, std::move(DATA)
	);
	//...
	HIDE::SCAL_OF_GRID[GRID].insert(S);
	HIDE::SCAL[S] = S;
	//...
	return S;
}

c_scal_scatter t_meth::new_scal_scatter(const c_grid_scatter &GRID, t_scal::t_data &&DATA) {
	//...
	h_scal_scatter S = new t_scal_scatter(
	GRID, std::move(DATA)
	);
	//...
	HIDE::SCAL_OF_GRID[GRID].insert(S);
	HIDE::SCAL[S] = S;
	//...
	return S;
}

c_scal_isoline t_meth::new_scal_isoline(const c_grid_isoline &GRID, t_real minz, t_real maxz) {
	//...
	h_scal_isoline S = new t_scal_isoline(
	GRID, t_scal::t_rect(minz, maxz,
	GRID->step().size())
	);
	//...
	HIDE::SCAL_OF_GRID[GRID].insert(S);
	HIDE::SCAL[S] = S;
	//...
	return S;
}

c_scal_analith t_meth::new_scal_analith(const t_scal::t_func &FUNC) {
	//...
	h_scal_analith S =
	new t_scal_analith(FUNC);
	//...
	HIDE::SCAL[S] = S;
	//...
	return S;
}

/** Генераторы трассеров **/

c_flux_regular t_meth::new_flux_regular(const c_scal_regular &SCAL, t_name name) {
	name = name.tolower(); c_flow FLOW = SCAL->grid()->flow();
	if (HIDE::FLUX_OF_FLOW[FLOW][name] != nullptr) {
		__ERR_METH("Can't overwrite existing flux object by same name!");
	}
	c_flux_regular F = new_flux_regular(SCAL);
	HIDE::FLUX_OF_FLOW[FLOW][name] = F;
	return F;
}

c_flux_isoline t_meth::new_flux_isoline(const c_scal_isoline &SCAL, t_name name) {
	name = name.tolower(); c_flow FLOW = SCAL->grid()->flow();
	if (HIDE::FLUX_OF_FLOW[FLOW][name] != nullptr) {
		__ERR_METH("Can't overwrite existing flux object by same name!");
	}
	c_flux_isoline F = new_flux_isoline(SCAL);
	HIDE::FLUX_OF_FLOW[FLOW][name] = F;
	return F;
}

c_flux_scatter t_meth::new_flux_scatter(const c_scal_scatter &SCAL, t_name name) {
	name = name.tolower(); c_flow FLOW = SCAL->grid()->flow();
	if (HIDE::FLUX_OF_FLOW[FLOW][name] != nullptr) {
		__ERR_METH("Can't overwrite existing flux object by same name!");
	}
	c_flux_scatter F = new_flux_scatter(SCAL);
	HIDE::FLUX_OF_FLOW[FLOW][name] = F;
	return F;
}

c_flux_scatter t_meth::new_flux_scatter(const c_grid_scatter &GRID, t_name name) {
	name = name.tolower(); c_flow FLOW = GRID->flow();
	if (HIDE::FLUX_OF_FLOW[FLOW][name] != nullptr) {
		__ERR_METH("Can't overwrite existing flux object by same name!");
	}
	c_flux_scatter F = new_flux_scatter(GRID);
	HIDE::FLUX_OF_FLOW[FLOW][name] = F;
	return F;
}

c_flux_regular t_meth::new_flux_regular(const c_scal_regular &SCAL) {
	h_flux_regular F =
	new t_flux_regular(SCAL);
	HIDE::FLUX[F] = F;
	return F;
}

c_flux_isoline t_meth::new_flux_isoline(const c_scal_isoline &SCAL) {
	h_flux_isoline F =
	new t_flux_isoline(SCAL);
	HIDE::FLUX[F] = F;
	return F;
}

c_flux_scatter t_meth::new_flux_scatter(const c_scal_scatter &SCAL) {
	h_flux_scatter F =
	new t_flux_scatter(SCAL);
	HIDE::FLUX[F] = F;
	return F;
}

c_flux_scatter t_meth::new_flux_scatter(const c_grid_scatter &GRID) {
	h_flux_scatter F =
	new t_flux_scatter(GRID);
	HIDE::FLUX[F] = F;
	return F;
}

//...

/** Обработчики служебных данных **/

h_data_convert t_meth::set_data_convert(const t_scal *SCAL, h_data_convert DATA) {
	auto dat = new_data_of_scal(SCAL);
	if (dat == nullptr) return nullptr;
	dat->CONVERT = DATA;
	return DATA;
}

h_data_correct t_meth::set_data_correct(const t_scal *SCAL, h_data_correct DATA) {
	auto dat = new_data_of_scal(SCAL);
	if (dat == nullptr) return nullptr;
	dat->CORRECT = DATA;
	return DATA;
}

h_data_dynamic t_meth::set_data_dynamic(const t_flow *FLOW, h_data_dynamic DATA) {
	auto dat = new_data_of_flow(FLOW);
	if (dat == nullptr) return nullptr;
	dat->DYNAMIC = DATA;
	return DATA;
}

h_data_inverse t_meth::set_data_inverse(const t_flow *FLOW, h_data_inverse DATA) {
	auto dat = new_data_of_flow(FLOW);
	if (dat == nullptr) return nullptr;
	dat->INVERSE = DATA;
	return DATA;
}

//...

h_data_of_scal t_meth::new_data_of_scal(const t_scal *SCAL) {
	auto it = HIDE::DATA_OF_SCAL.find(SCAL);
	if (it == HIDE::DATA_OF_SCAL.end()) {
		HIDE::DATA_OF_SCAL[SCAL] =
		new t_data_of_scal();
	}
	return HIDE::DATA_OF_SCAL[SCAL];
}

h_data_of_scal t_meth::get_data_of_scal(const t_scal *SCAL) {
	auto it = HIDE::DATA_OF_SCAL.find(SCAL);
	return
	(it != HIDE::DATA_OF_SCAL.end())?
	it->second:
	nullptr;
}

h_data_of_flow t_meth::new_data_of_flow(const t_flow *FLOW) {
	auto it = HIDE::DATA_OF_FLOW.find(FLOW);
	if (it == HIDE::DATA_OF_FLOW.end()) {
		HIDE::DATA_OF_FLOW[FLOW] =
		new t_data_of_flow();
	}
	return HIDE::DATA_OF_FLOW[FLOW];
}

h_data_of_flow t_meth::get_data_of_flow(const t_flow *FLOW) {
	auto it = HIDE::DATA_OF_FLOW.find(FLOW);
	return
	(it != HIDE::DATA_OF_FLOW.end())?
	it->second:
	nullptr;
}

//...

h_data_convert t_meth::get_data_convert(const t_scal *SCAL) {
	auto dat = get_data_of_scal(SCAL);
	return (dat != nullptr)?
	dat->CONVERT:
	nullptr;
}

h_data_correct t_meth::get_data_correct(const t_scal *SCAL) {
	auto dat = get_data_of_scal(SCAL);
	return (dat != nullptr)?
	dat->CORRECT:
	nullptr;
}

h_data_dynamic t_meth::get_data_dynamic(const t_flow *FLOW) {
	auto dat = get_data_of_flow(FLOW);
	return (dat != nullptr)?
	dat->DYNAMIC:
	nullptr;
}

h_data_inverse t_meth::get_data_inverse(const t_flow *FLOW) {
	auto dat = get_data_of_flow(FLOW);
	return (dat != nullptr)?
	dat->INVERSE:
	nullptr;
}

//...

const t_dict<t_name, c_flux> &t_meth::get_flux(const c_flow &FLOW) { 
	return HIDE::FLUX_OF_FLOW[FLOW];
}

c_flux t_meth::get_flux(const c_flow &FLOW, t_name name) {
	name = name.tolower();
	c_flux F =
	HIDE::FLUX_OF_FLOW
	[FLOW][name];
	return F;
}

//...

void t_meth::cpy_data(const t_grid *_from, const t_grid *_grid) {
	//...
}

void t_meth::cpy_data(const t_vect *_from, const t_vect *_vect) {
	//...
}

void t_meth::cpy_data(const t_scal *_from, const t_scal *_scal) {
	auto data = get_data_of_scal(_from); if (!data) return;
	auto h = new t_data_of_scal();
	h->CONVERT = data->CONVERT;
	h->CORRECT = data->CORRECT;
	HIDE::DATA_OF_SCAL[_scal] = std::move(h);
}

void t_meth::cpy_data(const t_flow *_from, const t_flow *_flow) {
	auto data = get_data_of_flow(_from); if (!data) return;
	auto h = new t_data_of_flow();
	h->DYNAMIC = data->DYNAMIC;
	h->INVERSE = data->INVERSE;
	HIDE::DATA_OF_FLOW[_flow] = std::move(h);
}

//...

/** Сборщики мусора **/

void t_meth::del_data(const t_flow *_flow) {
	//NOTE: Следует иметь в виду, что на данный момент все ссылки на _flow уже были освобождены!
	auto  &list = HIDE::GRID_OF_FLOW[_flow];
	while (list.size()) {
		auto it = list.begin();
		(*it)->hand().del();	//Удаляем привязанную сетку (ссылка также обнуляется)!
		list.erase(it);	//Удаляем ссылку из set'а!
	}
	HIDE::GRID_OF_FLOW.
	erase(_flow);
	HIDE::FLUX_OF_FLOW.
	erase(_flow);
	HIDE::DATA_OF_FLOW.
	erase(_flow);
	HIDE::FLOW.
	erase(_flow);
}

void t_meth::del_data(const t_grid *_grid) {
	//NOTE: Следует иметь в виду, что на данный момент все ссылки на _grid уже были освобождены!
	c_flow FLOW = _grid->flow();
	if (FLOW != nullptr) HIDE::GRID_OF_FLOW[FLOW].erase(_grid);
	//...
	auto  &vect = HIDE::VECT_OF_GRID[_grid];
	while (vect.size()) {
		auto it = vect.begin();
		(*it)->hand().del();	//Удаляем привязанное поле!
		vect.erase(it);	//Удаляем ссылку из set'а!
	}
	auto  &scal = HIDE::SCAL_OF_GRID[_grid];
	while (scal.size()) {
		auto it = scal.begin();
		(*it)->hand().del();	//Удаляем привязанное поле!
		scal.erase(it);	//Удаляем ссылку из set'а!
	}
	HIDE::VECT_OF_GRID.
	erase(_grid);
	HIDE::SCAL_OF_GRID.
	erase(_grid);
	HIDE::GRID.
	erase(_grid);
}

void t_meth::del_data(const t_vect *_vect) {
	//NOTE: Следует иметь в виду, что на данный момент все ссылки на _vect уже были освобождены!
	HIDE::VECT.erase(_vect);
}

void t_meth::del_data(const t_scal *_scal) {
	//NOTE: Следует иметь в виду, что на данный момент все ссылки на _scal уже были освобождены!
	HIDE::DATA_OF_SCAL.erase(_scal);
	HIDE::SCAL.erase(_scal);
}

void t_meth::del_data(const t_flux *_flux) {
	//NOTE: Следует иметь в виду, что на данный момент все ссылки на _flux уже были освобождены!
	HIDE::FLUX.erase(_flux);
}

//...

/** Возвращают обратные ссылки **/

h_flow t_flow::hand() const { return (HIDE::FLOW.count(this))? HIDE::FLOW[this]: nullptr; }

h_grid t_grid::hand() const { return (HIDE::GRID.count(this))? HIDE::GRID[this]: nullptr; }

h_vect t_vect::hand() const { return (HIDE::VECT.count(this))? HIDE::VECT[this]: nullptr; }

h_scal t_scal::hand() const { return (HIDE::SCAL.count(this))? HIDE::SCAL[this]: nullptr; }

h_flux t_flux::hand() const { return (HIDE::FLUX.count(this))? HIDE::FLUX[this]: nullptr; }

/** Деструкторы **/

t_flow::~t_flow() { t_meth::del_data(this); }

t_grid::~t_grid() { t_meth::del_data(this); }

t_vect::~t_vect() { t_meth::del_data(this); }

t_scal::~t_scal() { t_meth::del_data(this); }

t_flux::~t_flux() { t_meth::del_data(this); }

//...

}
