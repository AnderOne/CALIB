#include <calib/calib.hpp>

namespace CALIB {

//Строит карту изолиний на прямоугольной сетке:
c_scal_isoline t_meth::run_meth_convert(const c_scal_regular &FROM, const t_scal::t_rect &RECT) {

	//Определяем параметры потока:
	t_flow::t_cond cond = FROM->grid()->flow()->cond();

	//Определяем параметры сетки:
	t_real minx = FROM->grid()->rect().minx(), miny = FROM->grid()->rect().miny();
	t_real maxx = FROM->grid()->rect().maxx(), maxy = FROM->grid()->rect().maxy();
	t_real lenx = FROM->grid()->rect().lenx(), leny = FROM->grid()->rect().leny();
	t_size numx = FROM->grid()->rect().numx(), numy = FROM->grid()->rect().numy();
	t_real minz = RECT.minz(), maxz = RECT.maxz();
	t_size numz = RECT.numz();

	t_scal::t_data DATA(FROM->data());

	std::vector<t_real> BUFF(numz);
	t_real *valz = BUFF.data(); for (int i = 0; i < numz; ++ i) valz[i] = RECT.valz(i);

	//Класс итератора вершины:
	struct t_node_iter {
		inline t_node_iter(int _ix, int _iy): ix(_ix), iy(_iy) {}
		int ix, iy;
	};

	//Класс итератора ребра:
	struct t_edge_iter {
		inline t_edge_iter(int _ix, int _iy, int _dr):
		                   ix(_ix), iy(_iy), dr(_dr) {}
		inline t_node_iter get_head() const {
			return t_node_iter(ix + !(dr & 1), iy + (dr & 1));
		}
		inline t_node_iter get_tail() const {
			return t_node_iter(ix, iy);
		}
		int ix, iy;
		int dr;
	};

	//Класс итератора ячейки:
	struct t_cell_iter {
		inline t_cell_iter(int _ix, int _iy): ix(_ix), iy(_iy) {}
		inline t_edge_iter get_edge(int _id) const {
			static const int sx[] = {0, 1, 0, 0};
			static const int sy[] = {0, 0, 1, 0};
			static const int dr[] = {0, 1, 0, 1};
			return t_edge_iter(ix + sx[_id],
			                   iy + sy[_id],
			                   dr[_id]);
		}
		int ix, iy;
	};

	//Вспомогательные массивы (для контуров и узлов):
	std::vector<t_size> CONT_I, CONT_N; std::vector<t_bool> CONT_S;
	std::vector<t_real> NODE_X, NODE_Y;
	t_size numk, numc;

	//Определяем тип краевых условий:
	t_size ncx = numx - 1, ncy = numy - 1;	//Число ячеек по X и по Y;
	t_size mx, my;
	switch (cond) {
	case t_flow::t_cond::PERIOD2:
		mx = numx - 1; my = numy - 1;
		break;
	case t_flow::t_cond::PERIODX:
		mx = numx - 1; my = numy;
		break;
	case t_flow::t_cond::PERIODY:
		mx = numx; my = numy - 1;
		break;
	default:
		mx = numx; my = numy;
	}

	//Заполняем вспомогательные массивы:
	std::vector<t_long> HEAD_E(2 * mx * my), TAIL_E(2 * mx * my);	//Хранит индексы минимального и максимального уровня, связанного с текущим ребром;
	std::vector<t_long> NODE_E(2 * mx * my);	//Хранит индекс точки (в массиве NODE), лежащей на линии минимального уровня, которая проходит через текущее ребро;
	std::vector<t_long> STAT_E(2 * mx * my);	//Хранит информацию о направлении ребра;
	std::vector<t_long> NEXT_N, LAST_N;	//Хранит порядок обхода новых узлов;
	std::vector<t_long> STAT_N;	//Хранит соответствия между узлами и уровнями квантования;
	std::vector<t_long> NEXT_C, LAST_C;	//Хранит порядок обхода контуров;
	std::vector<t_long> HEAD_C(numz);	//Хранит индекс начального контура на текущем уровне;
	std::vector<t_long> STAT_C;	//Хранит соответствия между контурами и уровнями квантования;

	t_real x0, y0, z0, x1, y1, z1, x2, y2, z2, dx, dy, dz, s;
	t_long ix, iy, ix1, iy1, ix2, iy2, l1, l2, k, l, dr;
	//Корректируем сеточное поле, чтобы избежать вырожденных случаев:
	for (ix = 0; ix < mx; ++ ix)
	for (iy = 0; iy < my; ++ iy) if (std::binary_search(valz, valz + numz, DATA.valz(ix, iy))) DATA.valz(ix, iy) += 1.e-14;

	//Находит все пересечения с текущим ребром:
	auto cross = [&] (t_long ix, t_long iy, t_long dr) {
		t_edge_iter ed = t_edge_iter(ix, iy, dr); k = dr * (my) * (mx) + ix * (my) + iy;
		t_node_iter n1 = ed.get_tail();
		t_node_iter n2 = ed.get_head();
		ix1 = n1.ix; iy1 = n1.iy;
		ix2 = n2.ix; iy2 = n2.iy;
		z1 = DATA.valz((ix1 + mx) % mx, (iy1 + my) % my); z2 = DATA.valz((ix2 + mx) % mx, (iy2 + my) % my);
		dz = z2 - z1;
		if (dz < 0) {
			l1 = std::lower_bound(valz, valz + numz, z2) - valz;
			l2 = std::upper_bound(valz, valz + numz, z1) - valz;
			STAT_E[k] = - 1;
		}
		else
		if (dz > 0) {
			l1 = std::lower_bound(valz, valz + numz, z1) - valz;
			l2 = std::upper_bound(valz, valz + numz, z2) - valz;
			STAT_E[k] = + 1;
		}
		else {
			STAT_E[k] = 0;
			l1 = l2 = 0;
		}
		//Запоминаем диапазон уровней и направление ребра:
		NODE_E[k] = numk;
		TAIL_E[k] = l1;
		HEAD_E[k] = l2;
		//Запоминаем точки пересечения:
		x1 = minx + ix1 * lenx; y1 = miny + iy1 * leny;
		x2 = minx + ix2 * lenx; y2 = miny + iy2 * leny;
		dx = x2 - x1; dy = y2 - y1;
		for (l = l1; l < l2; ++ l) {
			s = (valz[l] - z1) / dz;
			NODE_X.push_back(x1 + dx * s);
			NODE_Y.push_back(y1 + dy * s);
			NEXT_N.push_back(-1);
			LAST_N.push_back(-1);
			STAT_N.push_back(l);
			++ numk;
		}
	};
	//Запоминаем точки пересечения:
	for (numk = k = 0, dr = 0; dr < 2; ++ dr)
	for (ix = 0; ix < ncx; ++ ix)
	for (iy = 0; iy < ncy; ++ iy) {
		cross(ix, iy, dr);
	}
	if ((cond != t_flow::t_cond::PERIOD2) &&
	    (cond != t_flow::t_cond::PERIODY)) {
	for (ix = 0; ix < ncx; ++ ix)
		cross(ix, ncy, 0);
	}
	if ((cond != t_flow::t_cond::PERIOD2) &&
	    (cond != t_flow::t_cond::PERIODX)) {
	for (iy = 0; iy < ncy; ++ iy)
		cross(ncx, iy, 1);
	}

	//Формируем списки точек:
	for (ix = 0; ix < ncx; ++ ix) for (iy = 0; iy < ncy; ++ iy) {
		static const int NORM[] = {+1, +1, -1, -1};
		static int L1[4], L2[4], IK[4], DR[4];
		t_cell_iter c1 = t_cell_iter(ix, iy);
		for (int i = 0; i < 4; ++ i) {
			t_edge_iter e1 = c1.get_edge(i);
			ix1 = (e1.ix + mx) % mx;
			iy1 = (e1.iy + my) % my;
			k = e1.dr * (my) * (mx) + ix1 * (my) + iy1;
			DR[i] = (STAT_E[k] == NORM[i]);
			L1[i] = TAIL_E[k];
			L2[i] = HEAD_E[k];
			IK[i] = NODE_E[k];
		}
		static int i1, i2, n1, l1, l2;
		for (int i = 0; i < 4; ++ i) {
			if (!DR[i]) continue;
			for (; L1[i] < L2[i]; ++ L1[i])
			for (int j = 1; j < 4; ++ j) {
				k = (i + j) % 4;
				if (DR[k] || (L2[k] <= L1[i])
				          || (L1[k] > L1[i])) {
					continue;
				}
				i1 = IK[i];
				i2 = IK[k] + L1[i] - L1[k];
				if (LAST_N[i2] != - 1) {
					continue;
				}
				//Формируем сегмент контура:
				LAST_N[i2] = i1;
				NEXT_N[i1] = i2;
				++ IK[i];
				break;
			}
			//assert(L1[i] == L2[i]);
		}
	}

	//Осуществляем сборку контуров (путем последовательного обхода узлов):
	std::vector<std::vector<std::pair<int, t_size> > > HEAD_S(numz);
	t_grid::t_node NODE(numk);
	t_long k0 = 0, k1 = 0;
	numc = 0;
	//Собираем незамкнутые контуры:
	for (t_long k, ik = 0; ik < numk; ++ ik) {
		if ((LAST_N[ik] != - 1) || (STAT_N[ik] == - 1)) continue;
		HEAD_S[STAT_N[ik]].push_back(std::make_pair(ik, numc));
		NODE.valx(k1) = NODE_X[ik]; NODE.valy(k1) = NODE_Y[ik];
		STAT_N[ik] = - 1;
		++ k1;
		for (k = NEXT_N[ik]; k != - 1; k = NEXT_N[k]) {
			NODE.valx(k1) = NODE_X[k];
			NODE.valy(k1) = NODE_Y[k];
			STAT_N[k] = - 1;
			++ k1;
		}
		//assert(k1 - k0 >= 3);
		CONT_N.push_back(k1 - k0);
		CONT_I.push_back(k0);
		CONT_S.push_back(0);
		++ numc;
		k0 = k1;
	}
	//Собираем замкнутые контуры:
	for (t_long k, ik = 0; ik < numk; ++ ik) {
		if (STAT_N[ik] == - 1) continue;
		HEAD_S[STAT_N[ik]].push_back(std::make_pair(ik, numc));
		NODE.valx(k1) = NODE_X[ik]; NODE.valy(k1) = NODE_Y[ik];
		STAT_N[ik] = - 1;
		dx = dy = 0;
		++ k1;
		for (k = NEXT_N[ik]; k != ik; k = NEXT_N[k]) {
			NODE.valx(k1) = NODE_X[k];
			NODE.valy(k1) = NODE_Y[k];
			STAT_N[k] = - 1;
			++ k1;
		}
		//assert(k1 - k0 >= 3);
		CONT_N.push_back(k1 - k0);
		CONT_I.push_back(k0);
		CONT_S.push_back(1);
		++ numc;
		k0 = k1;
	}
	//Сортируем массивы контуров:
	t_grid::t_cont CONT(numc);
	t_grid::t_step STEP(numz);
	for (t_long c0 = 0, c1 = 0, i = 0; i < numz; ++ i) {
		for (t_long j = 0; j < HEAD_S[i].size(); ++ j) {
			t_long c = HEAD_S[i][j].second;
			CONT.head(c1) = CONT_I[c];
			CONT.tail(c1) = CONT_N[c];
			CONT.stat(c1) = CONT_S[c];
			HEAD_S[i][j].second = c1;
			++ c1;
		}
		STEP.tail(i) = c1 - c0;
		STEP.head(i) = c0;
		c0 = c1;
	}
	//Формируем новое поле:
	c_grid_isoline GRID = t_meth::new_grid_isoline(
		FROM->grid()->flow(),
		std::move(NODE),
		std::move(CONT),
		std::move(STEP)
	);
	c_scal_isoline SCAL = t_meth::new_scal_isoline(
		GRID, minz, maxz
	);
	//Копируем параметры:
	cpy_data(FROM, SCAL);
	//...
	return SCAL;
}

}
