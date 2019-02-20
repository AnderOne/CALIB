#include <calib/calib.hpp>

namespace CALIB {

inline static c_vect_regular _new(const c_grid_regular &GRID, t_vect::t_data &&DATA) { return t_meth::new_vect_regular(GRID, std::move(DATA)); }

inline static c_scal_regular _new(const c_grid_regular &GRID, t_scal::t_data &&DATA) { return t_meth::new_scal_regular(GRID, std::move(DATA)); }

template <typename T_INP, typename T_OUT>
struct t_meth_convert_analith_to_regular:
public t_meth {

	inline static t_hand<const T_OUT> run(const t_hand<const T_INP> &FROM, const c_grid_regular &GRID) {

		t_real time = GRID->flow()->time();
		t_size numx = GRID->rect().numx();
		t_size numy = GRID->rect().numy(); typename T_OUT::t_data DATA(numx, numy);
		//...
		for (int ix = 0; ix < numx; ++ ix)
		for (int iy = 0; iy < numy; ++ iy) {
			DATA(ix, iy) = FROM->valz(
				GRID->rect().valx(ix), GRID->rect().valy(iy), time
			);
		}
		return _new(
			GRID, std::move(DATA)
		);
	}
};

c_vect_regular t_meth::run_meth_convert(const c_vect_analith &FROM, const c_grid_regular &GRID) {
	c_vect_regular VECT =
	t_meth_convert_analith_to_regular<
	t_vect_analith, t_vect_regular
	>::run(FROM, GRID);
	t_meth::cpy_data(FROM, VECT);
	return VECT;
}

c_scal_regular t_meth::run_meth_convert(const c_scal_analith &FROM, const c_grid_regular &GRID) {
	c_scal_regular SCAL =
	t_meth_convert_analith_to_regular<
	t_scal_analith, t_scal_regular
	>::run(FROM, GRID);
	t_meth::cpy_data(FROM, SCAL);
	return SCAL;
}

}
