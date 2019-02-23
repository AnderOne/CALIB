#include <calib/calib.hpp>
#include <fftw/fftw3.h>

#include <complex>

#define PI (3.14159265358979323846)

namespace CALIB {

//Абстрактный класс солвера для уравнения Пуассона:
struct t_meth_poisson  {
public:
	virtual t_scal::t_data get(const t_scal::t_data &RVORT) = 0;
	virtual ~t_meth_poisson() {}
};

//Класс солвера для случая двойной периодичности:
struct t_meth_poisson_period2:
public t_meth_poisson {

	t_meth_poisson_period2(const t_grid::t_rect &RECT) {
		minx = RECT.minx(); maxx = RECT.maxx();
		numx = RECT.numx() - 1;
		halx = numx / 2 + 1;
		miny = RECT.miny(); maxy = RECT.maxy();
		numy = RECT.numy() - 1;
		haly = numy / 2 + 1;
		//Выделяем память под вспомогательные массивы:
		DATA1.resize(halx, numy); GREEN.resize(halx, numy); GRIDZ.resize(numx, numy); GRIDS.resize(numx, numy);
		//Заполняем план преобразования:
		plan_r2c = fftw_plan_dft_r2c_2d(numy, numx, GRIDZ.data(), (fftw_complex *) DATA1.data(), FFTW_ESTIMATE);
		plan_c2r = fftw_plan_dft_c2r_2d(numy, numx, (fftw_complex *) DATA1.data(), GRIDZ.data(), FFTW_ESTIMATE);
		//Заполняем матрицу коэффициентов:
		t_real sx = 2.0 * PI / (maxx - minx), sy = 2.0 * PI / (maxy - miny), sx2 = sx * sx, sy2 = sy * sy;
		for (t_long iy1, iy = 0; iy < numy; ++ iy) for (t_long ix1, ix = 0; ix < halx; ++ ix) {
			ix1 = ix; iy1 = (2 * iy < numy)? (iy): (numy - iy);
			if ((ix1 != 0) || (iy1 != 0)) {
				GREEN(ix, iy) = - 1.0 / (sx2 * (ix1 * ix1) + sy2 * (iy1 * iy1));
			}
			else {
				GREEN(ix, iy) = 0;
			}
		}
	}

	t_scal::t_data get(const t_scal::t_data &RVORT) {
		//Исключаем значения граничных узлов:
		GRIDZ(t_sect::all(), t_sect::all()) = RVORT(t_sect(0, numx - 1), t_sect(0, numy - 1));
		//Выполняем спектральную инверсию:
		fftw_execute_dft_r2c(plan_r2c, GRIDZ.data(), (fftw_complex *) DATA1.data());
		DATA1() *= GREEN();
		fftw_execute_dft_c2r(plan_c2r, (fftw_complex *) DATA1.data(), GRIDS.data());
		GRIDS() /= (numx * numy);
		//Восстанавливаем граничные значения:
		t_scal::t_data SFUNC(numx + 1, numy + 1);
		SFUNC(t_sect(0, numx - 1), t_sect(0, numy - 1)) =
		GRIDS(t_sect::all(), t_sect::all());
		SFUNC(t_sect::all(), numy) =
		SFUNC(t_sect::all(), 0);
		SFUNC(numx, t_sect::all()) =
		SFUNC(0, t_sect::all());
		//...
		return SFUNC;
	}

	~t_meth_poisson_period2() {
		fftw_destroy_plan(plan_r2c); fftw_destroy_plan(plan_c2r);
	}
private:
	t_data<std::complex<t_real> > DATA1;
	t_data<t_real> GREEN, GRIDS, GRIDZ;
	fftw_plan plan_r2c, plan_c2r;
	t_real minx, miny, maxx, maxy; t_size numx, numy, halx, haly;
};

//Класс солвера для периодического канала:
struct t_meth_poisson_periodx:
public t_meth_poisson {

	t_meth_poisson_periodx(const t_grid::t_rect &RECT) {
		//...
		minx = RECT.minx(); maxy = RECT.maxy(); numx = RECT.numx() - 1; halx = numx / 2 + 1;
		miny = RECT.miny(); maxx = RECT.maxx(); numy = RECT.numy() - 2;
		//Выделяем память под вспомогательные массивы:
		DATA1.resize(halx, numy); DATA2.resize(numx, numy);
		GREEN.resize(halx, numy);
		GRIDZ.resize(numx, numy); GRIDS.resize(numx, numy);
		TMPX1.resize(numy); TMPX2.resize(numy);
		TMPY1.resize(halx); TMPY2.resize(numx);
		//Заполняем план преобразования:
		plan_r2c_x = fftw_plan_dft_r2c_1d(numx, TMPY2.data(), (fftw_complex *) TMPY1.data(), FFTW_ESTIMATE);
		plan_c2r_x = fftw_plan_dft_c2r_1d(numx, (fftw_complex *) TMPY1.data(), TMPY2.data(), FFTW_ESTIMATE);
		plan_r2r_y = fftw_plan_r2r_1d(numy, TMPX1.data(), TMPX2.data(), FFTW_RODFT00, FFTW_ESTIMATE);
		//Заполняем матрицу коэффициентов:
		t_real sx = 2.0 * PI / (maxx - minx), sy = PI / (maxy - miny), sx2 = sx * sx, sy2 = sy * sy;
		for (t_long iy1, iy = 0; iy < numy; ++ iy) for (t_long ix1, ix = 0; ix < halx; ++ ix) {
			ix1 = ix; iy1 = iy + 1;
			GREEN(ix, iy) = - 1.0 / (sx2 * (ix1 * ix1) + sy2 * (iy1 * iy1));
		}
	}

	t_scal::t_data get(const t_scal::t_data &RVORT) {
		//Выполняем прямое преобразование Фурье:
		for (int ix = 0; ix < numx; ++ ix) {
			for (int i = 0; i < numy; ++ i)
				TMPX1(i) = RVORT(ix, i + 1);
			fftw_execute_r2r(plan_r2r_y, TMPX1.data(), TMPX2.data());
			for (int i = 0; i < numy; ++ i)
				DATA2(ix, i) = TMPX2(i);
		}
		for (int iy = 0; iy < numy; ++ iy) {
			for (int i = 0; i < numx; ++ i)
				TMPY2(i) = DATA2(i, iy);
			fftw_execute_dft_r2c(plan_r2c_x, TMPY2.data(), (fftw_complex *) TMPY1.data());
			for (int i = 0; i < halx; ++ i)
				DATA1(i, iy) = TMPY1(i);
		}
		//Вычисляем функцию тока:
		DATA1() *= GREEN();
		for (int iy = 0; iy < numy; ++ iy) {
			for (int i = 0; i < halx; ++ i)
				TMPY1(i) = DATA1(i, iy);
			fftw_execute_dft_c2r(plan_c2r_x, (fftw_complex *) TMPY1.data(), TMPY2.data());
			for (int i = 0; i < numx; ++ i)
				DATA2(i, iy) = TMPY2(i);
		}
		DATA2() /= numx;
		for (int ix = 0; ix < numx; ++ ix) {
			for (int i = 0; i < numy; ++ i)
				TMPX2(i) = DATA2(ix, i);
			fftw_execute_r2r(plan_r2r_y, TMPX2.data(), TMPX1.data());
			for (int i = 0; i < numy; ++ i)
				GRIDS(ix, i) = TMPX1(i);
		}
		GRIDS() /= (2 * (numy + 1));
		//Восстанавливаем граничные значения:
		t_scal::t_data SFUNC(numx + 1, numy + 2);
		SFUNC(t_sect(0, numx - 1), t_sect(1, numy)) =
		GRIDS(t_sect::all(), t_sect::all());
		SFUNC(t_sect::all(), numy + 1) = 0;
		SFUNC(t_sect::all(), 0) = 0;
		SFUNC(numx, t_sect::all()) =
		SFUNC(0, t_sect::all());
		//...
		return SFUNC;
	}

	~t_meth_poisson_periodx() {
		fftw_destroy_plan(plan_r2r_y); fftw_destroy_plan(plan_r2c_x);
		fftw_destroy_plan(plan_c2r_x);
	}
private:
	t_data<std::complex<t_real> > DATA1, TMPY1;
	t_data<t_real> DATA2, GREEN, GRIDS, GRIDZ;
	t_data<t_real> TMPY2, TMPX1, TMPX2;
	fftw_plan plan_c2r_x, plan_r2r_y;
	fftw_plan plan_r2c_x;
	t_real minx, maxx, miny, maxy;
	t_size numx, numy, halx;
};

struct t_meth_poisson_periody:
public t_meth_poisson_periodx {

	t_meth_poisson_periody(const t_grid::t_rect &RECT):
	t_meth_poisson_periodx(
		t_grid::t_rect{RECT.miny(), RECT.maxy(), RECT.numy(), RECT.minx(), RECT.maxx(), RECT.numx()}
	) {}

	static t_scal::t_data transp(const t_scal::t_data &DATA) {
		t_scal::t_data TEMP(DATA.ncol(), DATA.nrow());
		for (int i = 0; i < DATA.nrow(); ++ i)
		for (int j = 0; j < DATA.ncol(); ++ j) {
			TEMP(j, i) = DATA(i, j);
		}
		return TEMP;
	}

	t_scal::t_data get(const t_scal::t_data &RVORT) {
		return transp(
			t_meth_poisson_periodx::get(
			transp(RVORT)
			)
		);
	}
};

//Класс солвера для замкнутого бассейна:
struct t_meth_poisson_period0:
public t_meth_poisson {

	t_meth_poisson_period0(const t_grid::t_rect &RECT) {
		//...
		minx = RECT.minx(); maxx = RECT.maxx(); numx = RECT.numx() - 2;
		miny = RECT.miny(); maxy = RECT.maxy(); numy = RECT.numy() - 2;
		//Выделяем память под вспомогательные массивы:
		DATA1.resize(numx, numy); DATA2.resize(numx, numy);
		GREEN.resize(numx, numy);
		GRIDZ.resize(numx, numy); GRIDS.resize(numx, numy);
		TMPX1.resize(numy); TMPX2.resize(numy);
		TMPY1.resize(numx); TMPY2.resize(numx);
		//Заполняем план преобразования:
		plan_r2r_y = fftw_plan_r2r_1d(numy, TMPX1.data(), TMPX2.data(), FFTW_RODFT00, FFTW_ESTIMATE);
		plan_r2r_x = fftw_plan_r2r_1d(numx, TMPY1.data(), TMPY2.data(), FFTW_RODFT00, FFTW_ESTIMATE);
		//Заполняем матрицу коэффициентов:
		t_real sx = PI / (maxx - minx), sy = PI / (maxy - miny), sx2 = sx * sx, sy2 = sy * sy;
		for (t_long iy1, iy = 0; iy < numy; ++ iy) for (t_long ix1, ix = 0; ix < numx; ++ ix) {
			ix1 = ix + 1; iy1 = iy + 1;
			GREEN(ix, iy) = - 1.0 / (sx2 * (ix1 * ix1) + sy2 * (iy1 * iy1));
		}
	}

	t_scal::t_data get(const t_scal::t_data &RVORT) {
		//Выполняем прямое преобразование Фурье:
		for (int ix = 0; ix < numx; ++ ix) {
			for (int i = 0; i < numy; ++ i)
				TMPX1(i) = RVORT(ix + 1, i + 1);
			fftw_execute_r2r(plan_r2r_y, TMPX1.data(), TMPX2.data());
			for (int i = 0; i < numy; ++ i)
				DATA2(ix, i) = TMPX2(i);
		}
		for (int iy = 0; iy < numy; ++ iy) {
			for (int i = 0; i < numx; ++ i)
				TMPY1(i) = DATA2(i, iy);
			fftw_execute_r2r(plan_r2r_x, TMPY1.data(), TMPY2.data());
			for (int i = 0; i < numx; ++ i)
				DATA1(i, iy) = TMPY2(i);
		}
		//Вычисляем функцию тока:
		DATA1() *= GREEN();
		for (int iy = 0; iy < numy; ++ iy) {
			for (int i = 0; i < numx; ++ i)
				TMPY2(i) = DATA1(i, iy);
			fftw_execute_r2r(plan_r2r_x, TMPY2.data(), TMPY1.data());
			for (int i = 0; i < numx; ++ i)
				DATA2(i, iy) = TMPY1(i);
		}
		DATA2() /= (2 * (numx + 1));
		for (int ix = 0; ix < numx; ++ ix) {
			for (int i = 0; i < numy; ++ i)
				TMPX2(i) = DATA2(ix, i);
			fftw_execute_r2r(plan_r2r_y, TMPX2.data(), TMPX1.data());
			for (int i = 0; i < numy; ++ i)
				GRIDS(ix, i) = TMPX1(i);
		}
		GRIDS() /= (2 * (numy + 1));
		//Восстанавливаем граничные значения:
		t_scal::t_data SFUNC(numx + 2, numy + 2);
		SFUNC(t_sect(1, numx), t_sect(1, numy)) =
		GRIDS(t_sect::all(), t_sect::all());
		SFUNC(t_sect::all(), numy + 1) = 0;
		SFUNC(t_sect::all(), 0) = 0;
		SFUNC(numx + 1, t_sect::all()) = 0;
		SFUNC(0, t_sect::all()) = 0;
		//...
		return SFUNC;
	}

	~t_meth_poisson_period0() {
		fftw_destroy_plan(plan_r2r_y); fftw_destroy_plan(plan_r2r_x);
	}
private:
	t_data<t_real> DATA1, DATA2, GREEN, GRIDS, GRIDZ;
	t_data<t_real> TMPY1, TMPY2, TMPX1, TMPX2;
	fftw_plan plan_r2r_x, plan_r2r_y;
	t_real minx, maxx, miny, maxy;
	t_size numx, numy;
};

//...

#define MAX_ITER (100)
#define EPS (1.e-7)

struct t_meth_inverse_of_casl:
public t_data_inverse,
public t_meth {

	inline explicit t_meth_inverse_of_casl(const c_grid_regular &_grid, t_flow::t_cond _cond): GRID(_grid), cond(_cond) {

		if (cond == t_flow::t_cond::PERIOD0) METH = new t_meth_poisson_period0(GRID->rect());
		if (cond == t_flow::t_cond::PERIODX) METH = new t_meth_poisson_periodx(GRID->rect());
		if (cond == t_flow::t_cond::PERIODY) METH = new t_meth_poisson_periody(GRID->rect());
		if (cond == t_flow::t_cond::PERIOD2) METH = new t_meth_poisson_period2(GRID->rect());
		//...
	}
	//...
	c_scal_regular run(const c_scal_regular &_scal) { return t_meth::run_meth_regular_jacobian(SFUNC, _scal); }

	c_vect_scatter run(const c_grid_scatter &_grid) { return t_meth::run_meth_convert(VELOC, _grid); }

	t_bool run(t_real dt) {

		c_scal_isoline PVORT = t_meth::get_flux_isoline(GRID->flow(), "PVORT")->scal();
		//Вычисляем поле относительной завихренности:
		RVORT = t_meth::run_meth_convert(PVORT, GRID);
		if (DEPTH != nullptr) { get_hand(RVORT)->valz() *= DEPTH->valz(); }
		if (FBETA != nullptr) { get_hand(RVORT)->valz() -= FBETA->valz(); }
		//Решаем уравнение Пуассона (без учета топографии):
		if (DEPTH == nullptr) {
			VELOC = run_meth_regular_velocity(SFUNC = new_scal_regular(GRID, METH->get(RVORT->data())));
			return true;
		}
		//Решаем модифицированное уравнение Пуассона:
		if (DIFFH == nullptr) { DIFFH = run_meth_regular_gradient(DEPTH); }
		t_scal::t_data RINIT(RVORT->data()); RINIT() *= DEPTH->valz();
		//...
		VELOC = t_meth::run_meth_regular_velocity(SFUNC = new_scal_regular(GRID, METH->get(RINIT)));
		get_hand(VELOC)->valz() /= DEPTH->valz();
		//Итерационный процесс:
		t_scal::t_data STEMP(SFUNC->data());
		t_scal::t_data RTEMP(RINIT);
		static int i; t_real r = 0;
		for (i = 0; i < MAX_ITER; ++ i) {
			RTEMP() = RINIT() + (DIFFH->valx() * VELOC->valy() - DIFFH->valy() * VELOC->valx());
			//Вычисляем значения поля скорости в узлах сетки:
			VELOC = run_meth_regular_velocity(
			SFUNC = new_scal_regular(GRID, METH->get(RTEMP))
			);
			get_hand(VELOC)->valz() /= DEPTH->valz();
			//Проверяем условие остановки:
			r = max(abs(SFUNC->valz() - STEMP()));
			if (r <= EPS * max(abs(STEMP()))) {
				break;
			}
			STEMP() = SFUNC->valz();
		}
		//TEST: For debug only!
		//std::cout << "[" << i << "] :: " << r << std::endl;
		return true;
	}

	t_bool run() {
		return true;
	}
protected:
	friend struct t_meth;
	t_hand<t_meth_poisson> METH;
	c_scal_regular RVORT;
	c_scal_regular SFUNC;
	c_vect_regular VELOC;
	c_scal_regular FBETA;
	c_scal_regular DEPTH;
	c_vect_regular DIFFH;
	t_real f0, bt;
	c_grid_regular GRID;
	t_flow::t_cond cond;
};

t_bool t_meth::set_data_inverse(const c_flow_plane2d &FLOW, const c_scal_isoline &PVORT,
                                                            const c_scal_regular &DEPTH, t_real f0, t_real bt) {

	set_data_inverse(FLOW, PVORT, DEPTH->grid(), f0, bt);

	t_hand<t_meth_inverse_of_casl> METH = get_data_inverse(FLOW).get<t_meth_inverse_of_casl> ();

	if ((f0 != 0) || (bt != 0)) {
		const t_grid::t_rect &RECT = DEPTH->grid()->rect();
		t_scal::t_data DATA(RECT.numx(), RECT.numy());
		for (t_size ix = 0; ix < RECT.numx(); ++ ix)
		for (t_size iy = 0; iy < RECT.numy(); ++ iy) {
			DATA(ix, iy) = f0 + bt * RECT.valy(iy);
		}
		METH->FBETA = t_meth::new_scal_regular(
			DEPTH->grid(), std::move(DATA)
		);
	}
	METH->DEPTH = DEPTH;

	return true;
}

t_bool t_meth::set_data_inverse(const c_flow_plane2d &FLOW, const c_scal_isoline &PVORT,
                                                            const c_grid_regular &GRID, t_real f0, t_real bt) {

	if ((PVORT->grid()->flow() != GRID->flow()) || (GRID->flow() != FLOW)) {
		__ERR_METH("Incompatible parameters of inversion! (objects belonging different flows)");
	}
	//Определяем параметры потока:
	t_flow::t_cond cond = FLOW->cond();
	//...
	if (IS_PERIODY(cond) && (bt != 0)) {
		__ERR_METH("This version not support y-periodic flows on a beta-plane!");
	}
	//Проверяем идентификаторы:
	if (get_flux(FLOW, "PVORT") != nullptr) {
		__ERR_METH("Special name of the flux \"PVORT\" is busy!");
	}
	t_meth::new_flux_isoline(
		PVORT, "PVORT"
	);
	//...
	t_hand<t_data_inverse> METH =
	new t_meth_inverse_of_casl(
		GRID, cond
	);
	set_data_inverse(
		FLOW, METH
	);
	//...
	return true;
}

//...

}
