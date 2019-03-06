#include <calib/calib.hpp>

#include <valarray>
#include <vector>
#include <list>

namespace CALIB {

#define eps3 (eps * eps * eps)
#define eps (1.e-14)

struct t_cell_sparse;
struct t_cell_dense;

//Шаблонный класс, содержащий реализацию процедуры хирургии контуров для различных типов краевых условий;
//Операции слияния/разбиения контуров объединены в операцию связывания;
//Для поиска близких узлов используется сеточный (клеточный) поиск;
//Диапазон ячеек определяется в пределах прямоугольной области;
//Поиск ведется до обнаружения наиближайшего сегмента контура;
template <typename t_cell, t_flow::t_cond cond> struct t_grid_correct_surgery: public t_meth {

	//Класс итератора для обхода близлежащих узлов:
	struct t_cell_iter {

		inline t_cell_iter(t_grid_correct_surgery &_SUR, t_long i): SUR(_SUR) {
			//Определяем текущую ячейку хирургической сетки:
			t_real cx = (SUR.NODE.valx(i) - SUR.minx) / SUR.hsx;
			t_real cy = (SUR.NODE.valy(i) - SUR.miny) / SUR.hsy;
			//Обзор узлов в близлежащих ячейках:
			ix1 = (t_long) (cx - SUR.rsx - eps);
			iy1 = (t_long) (cy - SUR.rsy - eps);
			ix2 = (t_long) (cx + SUR.rsx + eps);
			iy2 = (t_long) (cy + SUR.rsy + eps);
			if (!IS_PERIODX(cond) && !std::is_same<t_cell, t_cell_sparse>::value) {
				ix1 = std::min(std::max(ix1, 0), SUR.nsx - 1); ix2 = std::min(std::max(ix2, 0), SUR.nsx - 1);
			}
			if (!IS_PERIODY(cond) && !std::is_same<t_cell, t_cell_sparse>::value) {
				iy1 = std::min(std::max(iy1, 0), SUR.nsy - 1); iy2 = std::min(std::max(iy2, 0), SUR.nsy - 1);
			}
			ix = ix1 - 1;
			iy = iy1;
			k = - 1;
			next();
		}
		inline void next() {
			if (k >= 0) {
				k = SUR.TAILI(k); if (k >= 0) return;
			}
			++ ix;
			while (iy <= iy2) {
				while (ix <= ix2) {
					k = SUR.GRIDI(
						IS_PERIODX(cond)? ((ix + SUR.nsx) % SUR.nsx): ix,
						IS_PERIODY(cond)? ((iy + SUR.nsy) % SUR.nsy): iy
					);
					if (k >= 0) return;
					++ ix;
				}
				ix = ix1;
				++ iy;
			}
		}
		inline operator t_long() const { return k; }
	private:
		t_grid_correct_surgery &SUR;
		t_long ix1, iy1, ix2, iy2;
		t_long ix, iy;
		t_long k;
	};

	friend struct t_cell_iter;

	//...
	t_grid_correct_surgery(const c_scal_isoline &_from, t_real minl, t_real maxl):
	                       FROM(_from), NODE(_from->grid()->node()),
	                       DIFF(_from->grid()->node().size()) {

		//Параметры перераспределения:
		//delta = 2.0 * minl; ell = 0.25 * (maxl * maxl) / delta; mu = maxl / ell;
		delta = 2.0 * minl; mu = 8.0 * (minl / maxl); ell = maxl / mu;
		//...
		minx = FROM->grid()->flow()->rect().minx();
		miny = FROM->grid()->flow()->rect().miny();
		maxx = FROM->grid()->flow()->rect().maxx();
		maxy = FROM->grid()->flow()->rect().maxy();
		numk = FROM->grid()->node().size();
		numc = FROM->grid()->cont().size();
		numz = FROM->grid()->step().size();
		//...
		GEOM = t_flow::t_geom<cond, t_real> (minx, maxx, miny, maxy);
		//...
		DISTI.resize(numk);
		//...
		HEADI.resize(numk); TAILI.resize(numk);	//Списки узлов в ячейках "хирургической" сетки;
		NEXTI.resize(numk); LASTI.resize(numk);	//Списки узлов на контурах;
		STATI.resize(numk);
		//...
		//Оптимальный масштаб хирургической сетки:
		rad = sqrt(delta * delta + (mu * ell) * (mu * ell));
		if (FROM->grid()->flow()->rect().finx()) {
			nsx = (t_long) ((maxx - minx) / rad);
			hsx = (maxx - minx) / (nsx - 1);
		}
		else {
			minx = 0; hsx = rad;
		}
		if (FROM->grid()->flow()->rect().finy()) {
			nsy = (t_long) ((maxy - miny) / rad);
			hsy = (maxy - miny) / (nsy - 1);
		}
		else {
			miny = 0; hsy = rad;
		}
		GRIDI.resize(nsx, nsy);
		LISTI.resize(numz);
	}

	//Обрабатывает отдельно взятый уровень:
	void level(t_long iz) {

		const t_grid::t_step &STEP = FROM->grid()->step();
		const t_grid::t_cont &CONT = FROM->grid()->cont();

		//Заполняем вспомогательные массивы:
		GRIDI.clear();
		t_real r = 0; t_long c1, c2; c2 = (c1 = STEP.head(iz)) + STEP.tail(iz);
		for (t_long i0, i1, i2, ic = c1; ic < c2; ++ ic) {
			i2 = (i0 = CONT.head(ic)) + CONT.tail(ic);
			LISTI[iz].push_back(i0);
			//Вычисляем расстояния между узлами:
			t_real x2, y2;
			if (CONT.stat(ic)) { x2 = NODE.valx(i0); y2 = NODE.valy(i0); }
			else {
				x2 = NODE.valx(i2 - 1);
				y2 = NODE.valy(i2 - 1);
				DIFF.valx(i2 - 1) = 0;
				DIFF.valy(i2 - 1) = 0;
				DISTI(i2 - 1) = 0;
			}
			for (t_long ik = i2 - 1 - !CONT.stat(ic); ik >= i0; -- ik) {
				t_real x1 = NODE.valx(ik), tx1 = GEOM.subx(x2, x1);
				t_real y1 = NODE.valy(ik), ty1 = GEOM.suby(y2, y1);
				DIFF.valx(ik) = tx1;
				DIFF.valy(ik) = ty1;
				DISTI(ik) = tx1 * tx1 + ty1 * ty1;
				r = std::max(r, DISTI(ik));
				x2 = x1;
				y2 = y1;
			}
			//Заполняем списки узлов:
			for (t_long i = i0; i < i2; ++ i) {
				LASTI(i) = i - 1; NEXTI(i) = i + 1;
				HEADI(i) = TAILI(i) = - 1;
				insert_node(i);
				STATI(i) = 1;
			}
			if (CONT.stat(ic)) {
				LASTI(NEXTI(i2 - 1) = i0) = i2 - 1;
			}
			else {
				NEXTI(i2 - 1) = - 1;
				LASTI(i0) = - 1;
			}
		}
		//Определяем радиус хирургии:
		rad = sqrt(delta * delta + r) + eps;
		rsx = rad / hsx;
		rsy = rad / hsy;
		//Удаляем тонкие нити (филаменты):
		for (t_long i0, i1, i2, i, c = c1; c < c2; ++ c) {
			i = i0 = CONT.head(c);
			do {
				//Сокращаем тонкую нить:
				i1 = LASTI(i); i2 = NEXTI(i);
				shorten(i1, i, i2, 1);
				if (!STATI(i0)) i0 = i;
				i = NEXTI(i);
			} while (
				(i != i0) &&
				(i >= 0) &&
				STATI(i)
			);
			//Добавляем обработанную часть контура в список начальных узлов:
			if (STATI(i0)) {
				LISTI[iz].push_back(i0);
			}
		}
		//Поиск всех возможных сближений контуров +
		//связывание + удаление хирургических швов:
		//for (t_size n = 1, m = 1; (n != 0) && (m <= 1/*maxm*/); ++ m) {
		t_size n = 0, m = 1;
		//Обход массива исходных контуров:
		for (t_long c = c1; c < c2; ++ c) {
			t_long ii1, ii2;
			ii2 = (ii1 = CONT.head(c)) + CONT.tail(c);
			for (t_long ii = ii1; ii < ii2; ++ ii) {
				t_real rr, tx1, ty1, tx2, ty2;
				t_long j1, j2, i1, i2;
				//Если текущий узел был удален ранее:
				if (STATI(ii) != m) continue;
				t_long i = ii;
				i1 = LASTI(i);
				i2 = NEXTI(i);
				if ((i1 < 0) ||
				    (i2 < 0)) continue;
				//Ведем поиск близлежащего узла:
				t_real dd = 2.0 * delta;
				t_long jj = - 1;
				t_cell_iter ITER(*this, i);
				for (t_long j = ITER; j >= 0; ITER.next(), j = ITER) {
					//Исключаем из поиска смежные узлы:
					if ((j == i) || (j == i1) || (j == i2) || (LASTI(j) == i2)) continue;
					//Исключаем крайние узлы:
					if ((NEXTI(j) < 0) || (LASTI(j) < 0)) continue;
					//Вычисляем разность узлов:
					t_real tx1 = DIFF.valx(i);
					t_real ty1 = DIFF.valy(i);
					t_real tx2 = GEOM.subx(NODE.valx(i), NODE.valx(j));
					t_real ty2 = GEOM.suby(NODE.valy(i), NODE.valy(j));
					//Проверяем расстояние между узлами:
					t_real r = sqrt(tx2 * tx2 + ty2 * ty2);
					if (r > rad) continue;
					//Проверяем расстояние от узла до отрезка:
					t_real d = tx1 * tx2 + ty1 * ty2;
					d = (d * (d + DISTI(i)) <= 0)?
					    (std::abs(tx2 * ty1 - tx1 * ty2) / (sqrt(DISTI(i)) + eps)):
					    (r);
					if (d > delta) continue;
					//Определяем наиближайший узел:
					if ((d < dd) && (d != 0)) {
						dd = d; rr = r; jj = j;
					}
					/*//Когда таких узлов несколько:
					else
					if ((d == dd) && (j > jj)) {
						dd = d; rr = r; jj = j;
					}*/
				}
				//Если близкий узел не был найден:
				if (jj < 0) continue;
				t_long j = jj;
				n = n + 1;
				//Определяем ближайший узел текущего сегмента:
				tx2 = GEOM.subx(NODE.valx(i2), NODE.valx(j));
				ty2 = GEOM.suby(NODE.valy(i2), NODE.valy(j));
				dd = sqrt(tx2 * tx2 + ty2 * ty2);
				if ((dd < rr) && (dd != 0)) { tx1 = tx2; ty1 = ty2; i = i2; }
				//Выполняем операцию связывания контуров:
				i1 = LASTI(i); i2 = NEXTI(i);
				j1 = LASTI(j); j2 = NEXTI(j);
				connect(i1, i, i2, j1, j, j2, m + 1);
				//Обрабатываем хирургические швы на отделенной части контура:
				j1 = LASTI(j); j2 = NEXTI(j);
				shorten(j1, j, j2, m + 1);
				//Обрабатываем хирургические швы на текущей части контура:
				i1 = LASTI(i); i2 = NEXTI(i);
				shorten(i1, i, i2, m + 1);
				//Добавляем новые части контуров в список начальных узлов:
				if (STATI(i)) LISTI[iz].push_back(i);
				if (STATI(j)) LISTI[iz].push_back(j);
			}
		}
		//}
	}

	//Выполняет сборку результирующих контуров для получения новой карты изолиний:
	c_scal_isoline get() {
		//Выполняем процедуру сборки контуров:
		numk = 0;
		numc = 0;
		for (t_long iz = 0; iz < numz; ++ iz) {
			for (auto it = LISTI[iz].begin(); it != LISTI[iz].end(); ) {
				t_long i0 = *it;
				//Удаляем дубликаты ранее просмотренных контуров:
				if (!STATI(i0)) {
					it = LISTI[iz].erase(it); continue;
				}
				//Обходим узлы контура:
				t_long i = i0; t_size n = 0;
				do {
					//Помечаем все просмотренные узлы:
					STATI(i) = 0;
					i = NEXTI(i);
					++ n;
				} while ((i >= 0) && (i != i0));
				//Если встретили незамкнутый контур:
				if (i < 0) {
					//Выполняем обход назад:
					i = LASTI(i0);
					while (i >= 0) {
						STATI(i) = 0;
						i0 = i;
						i = LASTI(i);
						++ n;
					}
					*it = i0;
				}
				//Удаляем вырожденные части контуров:
				if (n < 3) {//mink) {
					it = LISTI[iz].erase(it);
					continue;
				}
				numk += n;
				numc += 1;
				++ it;
			}
		}
		//Выполняем процедуру сортировки контуров:
		if (!numk || !numc) {
			std::cerr << "WARNING: Result of surgery is empty!" << std::endl;
			//...
			t_grid::t_step STEP(numz); for (t_long i = 0; i < numz; ++ i) STEP.tail(i) = 0;
			return t_meth::new_scal_isoline(
			t_meth::new_grid_isoline(
				FROM->grid()->flow(), t_grid::t_node(),
				t_grid::t_cont(), STEP
			),
			FROM->rect().minz(),
			FROM->rect().maxz()
			);
		}
		t_grid::t_node NEW_NODE(numk);
		t_grid::t_cont NEW_CONT(numc);
		t_grid::t_step NEW_STEP(numz);
		for (t_long ik = 0, ic = 0, iz = 0; iz < numz; ++ iz) {
			//Обход новых контуров:
			NEW_STEP.head(iz) = ic;
			for (auto it = LISTI[iz].begin(); it != LISTI[iz].end(); ++ it) {
				t_long i0 = *it;
				NEW_CONT.head(ic) = ik;
				t_long i = i0;
				do {
					NEW_NODE.valx(ik) = NODE.valx(i);
					NEW_NODE.valy(ik) = NODE.valy(i);
					i = NEXTI(i);
					++ ik;
				} while ((i >= 0) && (i != i0));
				NEW_CONT.tail(ic) = ik - NEW_CONT.head(ic);
				NEW_CONT.stat(ic) = (i >= 0);
				++ ic;
			}
			NEW_STEP.tail(iz) = ic - NEW_STEP.head(iz);
		}
		//Формируем новое поле:
		c_grid_isoline GRID = t_meth::new_grid_isoline(
			FROM->grid()->flow(), std::move(NEW_NODE),
			std::move(NEW_CONT), std::move(NEW_STEP)
		);
		c_scal_isoline SCAL = t_meth::new_scal_isoline(
			GRID, FROM->rect().minz(),
			FROM->rect().maxz()
		);
		return SCAL;
	}

	c_scal_isoline run() {
		//Заполняем вспомогательные массивы:
		for (t_long iz = 0; iz < numz; ++ iz) level(iz);
		return get();
	}
private:
	//Производит связывание контуров:
	void connect(t_long &i1, t_long &i, t_long &i2, t_long &j1, t_long &j, t_long &j2, t_long m) {

		static t_real x1, y1, x2, y2, tx1, ty1, tx2, ty2;
		if ((i1 < 0) || (i2 < 0) || (j1 < 0) || (j2 < 0)) return;

		//Смещаем узлы:
		tx1 = GEOM.subx(NODE.valx(j), NODE.valx(i));
		ty1 = GEOM.suby(NODE.valy(j), NODE.valy(i));
		x1 = GEOM.movx(NODE.valx(i) + tx1 / 2.0);
		y1 = GEOM.movy(NODE.valy(i) + ty1 / 2.0);
		x2 = x1; y2 = y1;
		//Обновляем положение узлов:
		remove_node(i);
		remove_node(j);
		NODE.valx(i) = x1; NODE.valy(i) = y1; NODE.valx(j) = x2; NODE.valy(j) = y2;
		insert_node(i);
		insert_node(j);
		//Обновляем вспомогательные массивы:
		tx1 = GEOM.subx(x1, NODE.valx(i1)); ty1 = GEOM.suby(y1, NODE.valy(i1));
		DIFF.valx(i1) = tx1; DIFF.valy(i1) = ty1;
		DISTI(i1) = tx1 * tx1 + ty1 * ty1;
		tx1 = GEOM.subx(NODE.valx(j2), x1); ty1 = GEOM.suby(NODE.valy(j2), y1);
		DIFF.valx(i) = tx1; DIFF.valy(i) = ty1;
		DISTI(i)  = tx1 * tx1 + ty1 * ty1;
		tx2 = GEOM.subx(x2, NODE.valx(j1)); ty2 = GEOM.suby(y2, NODE.valy(j1));
		DIFF.valx(j1) = tx2; DIFF.valy(j1) = ty2;
		DISTI(j1) = tx2 * tx2 + ty2 * ty2;
		tx2 = GEOM.subx(NODE.valx(i2), x2); ty2 = GEOM.suby(NODE.valy(i2), y2);
		DIFF.valx(j) = tx2; DIFF.valy(j) = ty2;
		DISTI(j) = tx2 * tx2 + ty2 * ty2;
		//Меняем порядок следования узлов:
		LASTI(NEXTI(i) = j2) = i;
		LASTI(NEXTI(j) = i2) = j;
		//Помечаем смещенные сегменты:
		STATI(i1) = STATI(i) = m;
		STATI(j1) = STATI(j) = m;
	}

	//Сокращает тонкую нить:
	void shorten(t_long &i1, t_long &i, t_long &i2, t_long m) {

		static t_real tx1, ty1, tx2, ty2, x1, y1, d1, d2, d;
		static t_long i01, i02, k1, k2, j1, j2;
		//Если текущий узел был удален ранее:
		if (!STATI(i)) return;
		while (true) {
			if ((i1 < 0) || (i2 < 0) || (i1 == i2)) break;
			//КРАЙНИЕ УЗЛЫ НЕ ТРОГАЕМ!!!
			if ((LASTI(i1) < 0) || (NEXTI(i2) < 0)) break;
			//...
			tx1 = DIFF.valx(i1); tx2 = DIFF.valx(i);
			ty1 = DIFF.valy(i1); ty2 = DIFF.valy(i);
			x1 = NODE.valx(i); y1 = NODE.valy(i);
			d1 = DISTI(i1); d2 = DISTI(i);
			//Проверяем наличие острого угла:
			if (tx1 * tx2 + ty1 * ty2 >= 0) break;
			//Проверяем расстояние между узлами:
			d = std::abs(tx2 * ty1 - tx1 * ty2) / (std::max(sqrt(d1), sqrt(d2)) + eps);
			if (d * (d - delta) >= 0) break;
			//Определяем смещение углового узла:
			x1 = GEOM.movx(NODE.valx(i1) + (tx1 + tx2) / 2.0);
			y1 = GEOM.movy(NODE.valy(i1) + (ty1 + ty2) / 2.0);
			//Обновляем положение узлов:
			remove_node(i);
			NODE.valx(i) = x1; NODE.valy(i) = y1;
			insert_node(i);
			//Удаляем лишние узлы:
			remove_node(i1); STATI(i1) = 0;
			remove_node(i2); STATI(i2) = 0;
			//Если на контуре осталось 3 узла:
			if (NEXTI(i2) == i1) {
				remove_node(i2);
				STATI(i2) = 0;
				i2 = i1;
				break;
			}
			//Меняем порядок следования узлов:
			i1 = LASTI(i1); i2 = NEXTI(i2);
			if (i1 >= 0) NEXTI(i1) = i;
			if (i2 >= 0) LASTI(i2) = i;
			LASTI(i) = i1;
			NEXTI(i) = i2;
			//Смещаем сегменты контура:
			if (i1 >= 0) {
				tx1 = GEOM.subx(x1, NODE.valx(i1));
				ty1 = GEOM.suby(y1, NODE.valy(i1));
				d1  = tx1 * tx1 + ty1 * ty1;
				DIFF.valx(i1) = tx1;
				DIFF.valy(i1) = ty1;
				DISTI(i1) = d1;
			}
			if (i2 >= 0) {
				tx2 = GEOM.subx(NODE.valx(i2), x1);
				ty2 = GEOM.suby(NODE.valy(i2), y1);
				d2  = tx2 * tx2 + ty2 * ty2;
				DIFF.valx(i) = tx2;
				DIFF.valy(i) = ty2;
				DISTI(i) = d2;
			}
			else {
				DIFF.valx(i) = 0;
				DIFF.valy(i) = 0;
				DISTI(i) = 0;
			}
			//Помечаем смещенные сегменты:
			if (i1 >= 0) STATI(i1) = m;
			STATI(i) = m;
		}
		//Если остаточный контур вырожденный:
		if ((i1 < 0) && (i1 == i2)) {
			remove_node(i); STATI(i) = 0; LASTI(i) = -1; NEXTI(i) = -1;
			return;
		}
		if (i1 == i2) {
			if (i1 != i) {
				if (STATI(i1)) { remove_node(i1); STATI(i1) = 0; }
			}
			if (STATI(i)) { remove_node(i); STATI(i) = 0; }
			i1 = i2 = i;
			return;
		}
	}

	//Добавление узла в список:
	void insert_node(t_long i) {
		t_long ix = (t_long) ((NODE.valx(i) - minx) / hsx);
		t_long iy = (t_long) ((NODE.valy(i) - miny) / hsy);
		//assert((ix >= 0) && (ix < nsx) && (iy >= 0) && (iy < nsy));
		t_long ni = GRIDI(ix, iy); GRIDI(ix, iy) = i;
		if (ni >= 0) HEADI(TAILI(i) = ni) = i;
	}

	//Удаление узла из списка:
	void remove_node(t_long i) {
		t_long ix = (t_long) ((NODE.valx(i) - minx) / hsx);
		t_long iy = (t_long) ((NODE.valy(i) - miny) / hsy);
		//assert((ix >= 0) && (ix < nsx) && (iy >= 0) && (iy < nsy));
		t_long i1 = HEADI(i), i2 = TAILI(i);
		if (i1 < 0) GRIDI(ix, iy) = i2;
		else {
			TAILI(i1) = i2;
		}
		if (i2 >= 0) {
			HEADI(i2) = i1;
		}
		TAILI(i) = - 1;
		HEADI(i) = - 1;
	}

	t_real minx, miny, maxx, maxy, delta, ell, mu;
	t_size numk, numc, numz;
	//Массивы, содержащие характеристики контуров:
	std::vector<std::list<t_long> > LISTI;
	t_data<t_real> DISTI;
	t_data<t_long> STATI, NEXTI, LASTI;
	//Массивы, формирующие хирургическую сетку:
	t_data<t_long> HEADI, TAILI;
	t_cell GRIDI;
	//Параметры хирургической сетки:
	t_real rsx, rsy, rad;
	t_real hsx, hsy;
	t_long nsx, nsy;
	//...
	t_flow::t_geom<cond, t_real> GEOM;
	t_grid::t_node NODE, DIFF;
	//...
	c_scal_isoline FROM;
};

//Шаблонный класс, содержащий реализацию процедуры перераспределения для различных типов краевых условий;
//Контур аппроксимируется локальными кубическими сплайнами [Dritschel, 1989];
//Плотность распределения узлов зависит от кривизны [Dritschel, 1997];
//Угловые точки на контуре фиксируются;
template <t_flow::t_cond cond> struct t_meth_correct_renoder:
public t_meth {

	c_scal_isoline run(const c_scal_isoline &FROM) const {

		//Определяем параметры области:
		c_flow_plane2d FLOW = t_meth::get_flow_plane2d(FROM->grid()->flow());
		t_real minx = FLOW->rect().minx(), maxx = FLOW->rect().maxx();
		t_real miny = FLOW->rect().miny(), maxy = FLOW->rect().maxy();
		t_real minz = FROM->rect().minz(), maxz = FROM->rect().maxz();
		const t_grid::t_node &NODE = FROM->grid()->node();
		const t_grid::t_cont &CONT = FROM->grid()->cont();
		const t_grid::t_step &STEP = FROM->grid()->step();
		t_size numk = FROM->grid()->node().size();
		t_size numc = FROM->grid()->cont().size();
		t_size numz = FROM->grid()->step().size();

		//Параметры перераспределения:
		t_real delta = 2.0 * minl, mu = 8.0 * (minl / maxl), ell = maxl / mu;
		//...
		t_flow::t_geom<cond, t_real> GEOM(minx, maxx, miny, maxy);
		//...
		std::valarray<t_real> CUBSA(numk), CUBSB(numk), CUBSC(numk);
		std::valarray<t_real> DIFFX(numk), DIFFY(numk);
		std::valarray<t_real> DISTI(numk), CURVI(numk);
		std::vector<t_bool>   CORNI(numk);
		//...
		std::vector<t_real> NODEX, NODEY;
		std::vector<t_size> CONTI, CONTN;
		std::vector<t_bool> CONTS;
		std::vector<t_size> STEPI, STEPN;

		//Выполняем перераспределение узлов:
		for (t_long iz = 0; iz < numz; ++ iz) {
			t_real x2, y2; t_long i0, i1, i2, c0, c1, c2; c2 = (c1 = STEP.head(iz)) + STEP.tail(iz);
			c0 = CONTI.size();
			for (t_long ic = c1; ic < c2; ++ ic) {
				i2 = (i0 = CONT.head(ic)) + CONT.tail(ic);
				//Вычисляем расстояния между узлами:
				if (CONT.stat(ic)) { x2 = NODE.valx(i0); y2 = NODE.valy(i0); }
				else {
					x2 = NODE.valx(i2 - 1);
					y2 = NODE.valy(i2 - 1);
				}
				for (t_long ik = i2 - 1 - !CONT.stat(ic); ik >= i0; -- ik) {
					t_real x1 = NODE.valx(ik), tx1 = GEOM.subx(x2, x1);
					t_real y1 = NODE.valy(ik), ty1 = GEOM.suby(y2, y1);
					DIFFX[ik] = tx1;
					DIFFY[ik] = ty1;
					DISTI[ik] = tx1 * tx1 + ty1 * ty1 + eps;
					x2 = x1;
					y2 = y1;
				}
				//Вычисляем кривизну в узлах:
				t_real tx1 = DIFFX[i2 - 1], ty1 = DIFFY[i2 - 1], D1 = DISTI[i2 - 1];
				for (t_long i1 = i2 - 1, ik = i0; ik < i2; ++ ik) {
					t_real tx2 = DIFFX[ik], ty2 = DIFFY[ik], D2 = DISTI[ik];
					t_real W1 = tx1 * D2 + tx2 * D1;
					t_real W2 = ty1 * D2 + ty2 * D1;
					//Определяем наличие острого угла:
					if (tx1 * tx2 + ty1 * ty2 >= 0) {
						CURVI[ik] = 2.0 * (tx1 * ty2 - tx2 * ty1) /
						                   sqrt(W1 * W1 + W2 * W2 + eps3);
						CORNI[ik] = false;
					}
					else {
						CURVI[i1] = CURVI[ik] = 0;
						CORNI[ik] = true;
					}
					tx1 = tx2; ty1 = ty2; D1 = D2;
					i1 = ik;
				}
				//Концы незамкнутых контуров помечаем как угловые (не изменяемые) точки:
				if (!CONT.stat(ic)) CORNI[i0] = CORNI[i2 - 1] = 1;
				if (CORNI[i0]) CURVI[i2 - 1] = 0;
				//Вычисляем коэф-ты для кубической интерполяции:
				t_real K2 = CURVI[i0];
				for (t_long i1 = i0, ik = i2 - 1; ik >= i0; -- ik) {
					if (!CORNI[ik] && !CORNI[i1]) {
						t_real K1 = CURVI[ik], D1 = sqrt(DISTI[ik]);
						CUBSA[ik] = - D1 * K1 / 3.0 - D1  * K2 / 6.0;
						CUBSB[ik] =   D1 * K1 / 2.0;
						CUBSC[ik] =   D1 * (K2 - K1) / 6.0;
						K2 = K1;
						i1 = ik;
						continue;
					}
					//Вблизи острых углов используется
					//линейная интерполяция:
					CUBSA[ik] = CUBSB[ik] = CUBSC[ik] = 0;
					K2 = CURVI[ik];
					i1 = ik;
				}
				//Вычисляем кривизну дуги:
				t_real e = 1.0 / (ell * ell);
				K2 = CURVI[i0];
				for (t_long ik = i2 - 1; ik >= i0; -- ik) {
					t_real K1 = CURVI[ik]; t_real r = (K1 + K2) / 2.0;
					CURVI[ik] = sqrt(e + r * r);
					K2 = K1;
				}
				//Вычисляем взвешенную кривизну:
				t_real d = 4.0 * delta * delta;
				D1 = DISTI[i2 - 1];
				t_real W1 = sqrt(D1) / (D1 + d);
				t_real K1 = CURVI[i2 - 1];
				for (t_long ik = i0; ik < i2; ++ ik) {
					t_real D2 = DISTI[ik], W2 = sqrt(D2) / (D2 + d);
					K2 = CURVI[ik];
					CURVI[ik] = (W1 * K1 + W2 * K2) / (W1 + W2);
					W1 = W2;
					K1 = K2;
				}
				//Вычисляем необходимое число узлов на каждом интервале:
				d = 2.0 / delta;
				e = mu * sqrt(ell);
				K2 = CURVI[i0];
				for (t_long ik = i2 - 1; ik >= i0; -- ik) {
					K1 = CURVI[ik]; K2 = (K1 + K2) / 2.0;
					//DISTI[ik] = sqrt(DISTI[ik]) / (mu * ell); //Сохраняет постоянный шаг вдоль контура;
					DISTI[ik] = sqrt(DISTI[ik]) *
					            fmin(d, sqrt(K2) / e + K2);
					K2 = K1;
				}
				//Выполняем перераспределение узлов на контурах:
				t_long n2 = CONT.tail(ic);
				t_long ik = i0;
				//Определяем индекс начального угла:
				t_long i1 = i0;
				//Выполняем обход старых узлов:
				t_long k0 = NODEX.size();
				t_long n1 = 0;
				t_long i = i1;
				while (n1 < n2) {
					//Добавляем начальную точку:
					NODEX.push_back(NODE.valx(i));
					NODEY.push_back(NODE.valy(i));
					t_long m1 = n1;
					//Определяем следующий угол:
					t_real s = 0;
					do { s += DISTI[i]; i = i0 + ((i1 + (++ n1)) - i0) % n2; } while ((n1 < n2) && !CORNI[i]);
					if (s < eps) continue;
					//Добавляем новые узлы там, где это требуется:
					t_long j = i0 + ((i1 + m1) - i0) % n2; t_long m2 = (n1 < n2)? n1: n2;
					t_long n = (t_long) (s) + 1, m = NODEX.size() + n - 1;
					t_real d, p, h = n / s;
					j = i0 + ((i1 + m1 - 1) - i0) % n2;
					s = 0;
					while (NODEX.size() < m) {
						while (s < 1) {
							if (m1 < m2) { j = i0 + (i1 + (m1 ++) - i0) % n2; }
							s += (d = h * DISTI[j]);
						}
						p = 1 - (s -= 1) / d;
						//Интерполяция локальным кубическим сплайном [Dritschel, 1989]:
						t_real tx1 = DIFFX[j];
						t_real ty1 = DIFFY[j];
						t_real eta = p * (CUBSA[j] + p * (CUBSB[j] + p * CUBSC[j]));
						t_real x = NODE.valx(j) + p * tx1 - eta * ty1;
						t_real y = NODE.valy(j) + p * ty1 + eta * tx1;
						NODEX.push_back(GEOM.movx(x));
						NODEY.push_back(GEOM.movy(y));
					}
				}
				//Добавляем часть контура на текущем уровне:
				if (NODEX.size() - k0 < 3) {	//min(numk)!!!
					NODEX.resize(k0);
					NODEY.resize(k0);
					continue;
				}
				CONTN.push_back(NODEX.size() - k0);
				CONTI.push_back(k0);
				CONTS.push_back(CONT.stat(ic));
			}
			STEPN.push_back(CONTI.size() - c0);
			STEPI.push_back(c0);
		}

		//Формируем новый набор изолиний:
		t_grid::t_node NODE_(NODEX.size());
		for (t_long i = 0; i < NODE_.size(); ++ i) {
			NODE_.valx(i) = NODEX[i];
			NODE_.valy(i) = NODEY[i];
		}
		t_grid::t_cont CONT_(CONTI.size());
		for (t_long i = 0; i < CONT_.size(); ++ i) {
			CONT_.head(i) = CONTI[i];
			CONT_.tail(i) = CONTN[i];
			CONT_.stat(i) = CONTS[i];
		}
		t_grid::t_step STEP_(STEPI.size());
		for (t_long i = 0; i < STEP_.size(); ++ i) {
			STEP_.head(i) = STEPI[i];
			STEP_.tail(i) = STEPN[i];
		}
		//...
		c_grid_isoline GRID =
		t_meth::new_grid_isoline(
			FLOW,
			std::move(NODE_),
			std::move(CONT_),
			std::move(STEP_)
		);
		//...
		c_scal_isoline SCAL =
		t_meth::new_scal_isoline(
			GRID, minz, maxz
		);
		//...
		return SCAL;
	}

	inline t_meth_correct_renoder(t_real _minl, t_real _maxl): minl(_minl), maxl(_maxl) {}

private:
	t_real minl, maxl;
};

//Структура разряженной таблицы поиска:
struct t_cell_sparse {
	//...
	inline t_long &operator() (t_long i, t_long j) { return DATA[i][j]; }
	inline void resize(t_size n, t_size m) {}
	inline void clear() { DATA.clear(); }
	struct t_item {
		inline operator t_long &() { return val; }
	private:
		t_long val = - 1;
	};
private:
	std::map<t_long, std::map<t_long, t_item> > DATA;
};

//Структура плотной таблицы поиска:
struct t_cell_dense {
	//...
	inline t_long &operator() (t_long i, t_long j) { return DATA(i, j); }
	inline void resize(t_size n, t_size m) {
		DATA.resize(n, m);
	}
	inline void clear() {
		DATA() = - 1;
	}
private:
	t_data<t_long> DATA;
};

//Выполняет коррекцию контуров:
c_scal_isoline t_meth::run_meth_correct(const c_scal_isoline &FROM, t_type_correct type) {

	//Определяем параметры текущего потока:
	t_flow::t_cond cond = FROM->grid()->flow()->cond();

	//Определяем параметры коррекции:
	auto datc = get_data_correct(FROM);
	if (!datc) __ERR_METH("You must first call t_meth::set_data_correct(...)!");
	t_real minl = datc->minl;
	t_real maxl = datc->maxl;

	c_scal_isoline SCAL;
	//Выполняет коррекцию топологии контуров (хирургию):
	if (type == t_type_correct::SURGERY) {
		//Отдельный случай для бесконечной области:
		if (!FROM->grid()->flow()->rect().fin2()) {
			if (cond == t_flow::t_cond::PERIOD0)
				SCAL = t_grid_correct_surgery<t_cell_sparse, t_flow::t_cond::PERIOD0>
				(FROM, minl, maxl).run();
			if (cond == t_flow::t_cond::PERIODX)
				SCAL = t_grid_correct_surgery<t_cell_sparse, t_flow::t_cond::PERIODX>
				(FROM, minl, maxl).run();
			if (cond == t_flow::t_cond::PERIODY)
				SCAL = t_grid_correct_surgery<t_cell_sparse, t_flow::t_cond::PERIODY>
				(FROM, minl, maxl).run();
			if (cond == t_flow::t_cond::PERIOD2)
				SCAL = t_grid_correct_surgery<t_cell_sparse, t_flow::t_cond::PERIOD2>
				(FROM, minl, maxl).run();
		}
		else {
			if (cond == t_flow::t_cond::PERIOD0)
				SCAL = t_grid_correct_surgery<t_cell_dense, t_flow::t_cond::PERIOD0>
				(FROM, minl, maxl).run();
			if (cond == t_flow::t_cond::PERIODX)
				SCAL = t_grid_correct_surgery<t_cell_dense, t_flow::t_cond::PERIODX>
				(FROM, minl, maxl).run();
			if (cond == t_flow::t_cond::PERIODY)
				SCAL = t_grid_correct_surgery<t_cell_dense, t_flow::t_cond::PERIODY>
				(FROM, minl, maxl).run();
			if (cond == t_flow::t_cond::PERIOD2)
				SCAL = t_grid_correct_surgery<t_cell_dense, t_flow::t_cond::PERIOD2>
				(FROM, minl, maxl).run();
		}
	}
	else
	//Выполняет перераспределение узлов на контурах:
	if (type == t_type_correct::RENODER) {
		if (cond == t_flow::t_cond::PERIOD0)
			SCAL = t_meth_correct_renoder<t_flow::t_cond::PERIOD0>
			(minl, maxl).run(FROM);
		if (cond == t_flow::t_cond::PERIODX)
			SCAL = t_meth_correct_renoder<t_flow::t_cond::PERIODX>
			(minl, maxl).run(FROM);
		if (cond == t_flow::t_cond::PERIODY)
			SCAL = t_meth_correct_renoder<t_flow::t_cond::PERIODY>
			(minl, maxl).run(FROM);
		if (cond == t_flow::t_cond::PERIOD2)
			SCAL = t_meth_correct_renoder<t_flow::t_cond::PERIOD2>
			(minl, maxl).run(FROM);
	}
	else {
		__ERR_METH("Unknown type of corrector!");
	}
	//Копируем параметры:
	cpy_data(FROM, SCAL);

	return SCAL;
}

//Применяет процедуру ко всем трассерам:
t_bool t_meth::run_meth_correct(const c_flow_plane2d &FLOW, t_type_correct type) {

	std::set<c_flux> FLUX; for (auto it: get_flux(FLOW)) FLUX.insert(it.second);
	while (FLUX.size()) {
		auto it = FLUX.begin();
		c_flux_isoline _flx = get_flux_isoline(*it); FLUX.erase(it);
		if (_flx == nullptr)
			continue;
		c_scal_isoline _scl = run_meth_correct(_flx->scal(), type);
		_flx->hand().mov(new_flux_isoline(_scl)->hand());
	}
	return true;
}

//Инициализирует корректор:
t_bool t_meth::set_data_correct(const c_scal &SCAL, t_real minl, t_real maxl) {
	t_hand<t_data_correct> DATA(
		new t_data_correct(minl, maxl)
	);
	set_data_correct(SCAL, DATA);
	return true;
}

}
