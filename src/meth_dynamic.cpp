#include <calib/calib.hpp>

namespace CALIB {

template <t_flow::t_cond cond> struct t_meth_dynamic: public t_meth {

	t_meth_dynamic(const c_flow_plane2d &_flow, const t_hand<t_data_dynamic> &_meth_dyn, const t_hand<t_data_inverse> &_meth_inv):
	               FLOW(_flow), METH_DYN(_meth_dyn), METH_INV(_meth_inv) {

		const t_flow::t_rect &R = FLOW->rect(); GEOM = t_flow::t_geom<cond, t_real> (R.minx(), R.maxx(), R.miny(), R.maxy());
		//...
		dt1 = METH_DYN->lent; dt2 = dt1 / 2.0; dt3 = dt1 / 3.0; dt6 = dt1 / 6.0;
		//Запоминаем поля трассеров:
		for (auto it: get_flux(FLOW)) {
			c_flux_scatter FS = get_flux_scatter(it.second);
			if (FS != nullptr) { GRID_OF_FLUX.insert(get_hand(FS->grid())); continue; }
			c_flux_isoline FI = get_flux_isoline(it.second);
			if (FI != nullptr) { GRID_OF_FLUX.insert(get_hand(FI->grid())); continue; }
			c_flux_regular FR = get_flux_regular(it.second);
			if (FR != nullptr) { SCAL_OF_FLUX.insert(get_hand(FR->scal())); continue; }
		}
		for (auto it: GRID_OF_FLUX) {
			TEMP_OF_GRID[it] = t_grid::t_node(it->node().size());
			INIT_OF_GRID[it] = it->node();
		}
		for (auto it: SCAL_OF_FLUX) {
			TEMP_OF_SCAL[it] = t_scal::t_data(it->data().size());
			INIT_OF_SCAL[it] = it->data();
		}
	}

	//Выполняет один шаг расчетной схемы Рунге-Кутты:
	virtual t_bool run(int step) {

		std::map<const t_scal *, c_scal_regular> JACB_OF_SCAL;
		std::map<const t_grid *, c_vect_scatter> VECT_OF_GRID;

		//NOTE: Важно заполнить эти массивы до того, как будут расчитываться перемещения трассеров,
		//т.к. в некоторых случаях положение трассера влияет на поле скорости!
		for (h_scal_regular it: SCAL_OF_FLUX) {
			JACB_OF_SCAL[it] = METH_INV->run(it);
		}
		for (h_grid_scatter it: GRID_OF_FLUX) {
			VECT_OF_GRID[it] = METH_INV->run(it);
		}

		//Перенос эйлеровых полей:
		for (h_scal_regular it: SCAL_OF_FLUX) {
			//...
			t_scal::t_data &DATA = const_cast<t_scal::t_data &> (it->data());
			t_scal::t_data &TEMP = TEMP_OF_SCAL[it];
			t_scal::t_data &INIT = INIT_OF_SCAL[it];
			t_scal::t_data JACB = JACB_OF_SCAL[it]->data();
			//...
			if (step == 1) {
				TEMP() = INIT() - dt6 * JACB();
				DATA() = INIT() - dt2 * JACB();
			}
			if (step == 2) {
				TEMP() = TEMP() - dt3 * JACB();
				DATA() = INIT() - dt2 * JACB();
			}
			if (step == 3) {
				TEMP() = TEMP() - dt3 * JACB();
				DATA() = INIT() - dt1 * JACB();
			}
			if (step == 4) {
				DATA() = TEMP() - dt6 * JACB();
			}
		}
		//Перенос лагранжевых полей:
		for (h_grid_scatter it: GRID_OF_FLUX) {
			//...
			t_grid::t_node &NODE = const_cast<t_grid::t_node &> (it->node());
			t_grid::t_node &TEMP = TEMP_OF_GRID[it];
			t_grid::t_node &INIT = INIT_OF_GRID[it];
			t_vect::t_data VECT = VECT_OF_GRID[it]->data();
			//...
			if (step == 1) {
				TEMP() = INIT() + dt6 * VECT();
				NODE() = INIT() + dt2 * VECT();
			}
			if (step == 2) {
				TEMP() = TEMP() + dt3 * VECT();
				NODE() = INIT() + dt2 * VECT();
			}
			if (step == 3) {
				TEMP() = TEMP() + dt3 * VECT();
				NODE() = INIT() + dt1 * VECT();
			}
			if (step == 4) {
				NODE() = TEMP() + dt6 * VECT();
			}
			for (t_size i = 0; i < NODE.size(); ++ i) {
				NODE.valx(i) =
				GEOM.movx(NODE.valx(i));
				NODE.valy(i) =
				GEOM.movy(NODE.valy(i));
			}
		}
		//...
		return true;
	}
	//Выполняет интегрирование методом Рунге-Кутты:
	virtual t_bool run() {
		const t_real ST[] = {0, 0.5, 0.5, 1.0}, mint = FLOW->time();
		for (int t = 0; t < 4; ++ t) {
			t_meth::get_hand(FLOW)->time(mint + ST[t] * dt1);
			METH_INV->run(ST[t]);
			this->run(t);
		}
		METH_INV->run();
		return true;
	}
private:
	std::map<const t_grid *, t_grid::t_node> TEMP_OF_GRID;
	std::map<const t_grid *, t_grid::t_node> INIT_OF_GRID;
	std::map<const t_scal *, t_scal::t_data> TEMP_OF_SCAL;
	std::map<const t_scal *, t_scal::t_data> INIT_OF_SCAL;
	//...
	std::set<h_grid_scatter> GRID_OF_FLUX;
	std::set<h_scal_regular> SCAL_OF_FLUX;
	//...
	h_data_dynamic METH_DYN;
	h_data_inverse METH_INV;
	//...
	t_flow::t_geom<cond, t_real> GEOM;
	//...
	c_flow_plane2d FLOW;
	//...
	t_real dt1, dt2, dt3, dt6;
};

t_bool t_meth::run_meth_dynamic(const c_flow_plane2d &FLOW) {

	//Определяем параметры текущего потока:
	t_flow::t_cond cond = FLOW->cond();

	//Определяем параметры адвекции:
	auto data = get_data_dynamic(FLOW);
	if (!data) __ERR_METH("You must first call t_meth::set_data_dynamic(...)!");

	auto dati = get_data_inverse(FLOW);
	if (!dati) __ERR_METH("You must first call t_meth::set_data_inverse(...)!");

	if (cond == t_flow::t_cond::PERIOD0)
		return t_meth_dynamic<t_flow::t_cond::PERIOD0>
		(FLOW, data, dati).run();
	if (cond == t_flow::t_cond::PERIODX)
		return t_meth_dynamic<t_flow::t_cond::PERIODX>
		(FLOW, data, dati).run();
	if (cond == t_flow::t_cond::PERIODY)
		return t_meth_dynamic<t_flow::t_cond::PERIODY>
		(FLOW, data, dati).run();
	if (cond == t_flow::t_cond::PERIOD2)
		return t_meth_dynamic<t_flow::t_cond::PERIOD2>
		(FLOW, data, dati).run();
	//...
	return false;

}

//Инициализирует адвектор:
t_bool t_meth::set_data_dynamic(const c_flow &FLOW, t_real lent) {
	auto dat = get_data_of_flow(FLOW);
	if (dat == nullptr) return false;
	//...
	dat->DYNAMIC =
	new t_data_dynamic(lent);
	return true;
}

}
