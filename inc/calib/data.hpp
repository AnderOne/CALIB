#ifndef __INCLUDE_DATA_H
#define __INCLUDE_DATA_H

#include <blitz/array.h>
#include <algorithm>
#include <vector>
#include <map>

namespace CALIB {

typedef unsigned char t_bool;
typedef unsigned char t_byte;
typedef uint32_t t_size;
typedef int32_t t_long;
typedef double t_real; static const t_real INF = std::numeric_limits<t_real>::infinity();

template <typename TYPE, t_size SIZE> using t_item = blitz::TinyVector<TYPE, SIZE>;

template <typename K, typename V> using t_dict = std::map<K, V>;

template <typename V> using t_list = std::vector<V>;

struct t_name: public std::string {
	template <typename ... T> t_name(T ... args): std::string(args ...) {}
	inline t_name tolower() const {
		t_name s(*this);
		std::transform(s.begin(), s.end(), s.begin(), ::tolower);
		return s;
	}
	inline t_name toupper() const {
		t_name s(*this);
		std::transform(s.begin(), s.end(), s.begin(), ::toupper);
		return s;
	}
};

typedef blitz::Range t_sect;

template <typename TYPE,
t_size SIZE = 1>
struct t_data: public t_data<CALIB::t_item<TYPE, SIZE> > {

	typedef CALIB::t_item<TYPE, SIZE> t_item;

	template <typename ... T>
	t_data(T ... args):
	t_data<t_item>
	(args ...) {}
};

template <typename TYPE>
struct t_data<TYPE, 1> {

	inline t_data(const t_data &_rhs): DATA(new TYPE[_rhs.num]), num(_rhs.num), row(_rhs.row), col(_rhs.col) {
		std::copy(_rhs.DATA, _rhs.DATA + num, DATA);
	}

	inline t_data(t_data &&_rhs): DATA(_rhs.DATA), num(_rhs.num), row(_rhs.row), col(_rhs.col) {
		_rhs.DATA = nullptr;
		_rhs.num = _rhs.row = _rhs.col = 0;
	}

	inline explicit t_data(t_size _row, t_size _col):
	                DATA(new TYPE[_row * _col]), num(_row * _col), row(_row), col(_col) {}

	inline explicit t_data(t_size _num):
	                DATA(new TYPE[_num]), num(_num), row(_num), col(1) {}

	inline explicit t_data():
	                DATA(nullptr), num(0), row(0), col(0) {}

	inline virtual ~t_data() {
		delete[] DATA;
	}

	//Assignment operators:
	inline t_data &operator=(const t_data &_rhs) {
		delete[] DATA; DATA = new TYPE[_rhs.num];
		num = _rhs.num;
		row = _rhs.row;
		col = _rhs.col;
		std::copy(_rhs.DATA, _rhs.DATA + num, DATA);
		return *this;
	}

	inline t_data &operator=(t_data &&_rhs) {
		delete[] DATA;
		num = _rhs.num;
		row = _rhs.row;
		col = _rhs.col;
		DATA = _rhs.DATA; _rhs.DATA = nullptr;
		_rhs.num = 0;
		_rhs.row = 0;
		_rhs.col = 0;
		return *this;
	}

	//...
	inline void resize(t_size _row, t_size _col) {
		delete[] DATA;
		DATA = new TYPE[_row * _col];
		num = _row * _col;
		row = _row;
		col = _col;
	}

	inline void resize(t_size _num) {
		delete[] DATA;
		DATA = new TYPE[_num];
		num = _num;
		row = _num;
		col = 1;
	}

	//Nested types:
	typedef blitz::Array<TYPE, 2>
	t_part_2d;
	typedef blitz::Array<TYPE, 1>
	t_part_1d;
	typedef TYPE
	t_item;

	//Slicers:
	inline const t_part_2d operator()(t_size _row, const t_sect &_col) const {
		return t_part_2d(DATA, blitz::shape(row, col), blitz::neverDeleteData, blitz::ColumnMajorArray<2>())(t_sect(_row, _row), _col);
	}
	inline t_part_2d operator()(t_size _row, const t_sect &_col) {
		return t_part_2d(DATA, blitz::shape(row, col), blitz::neverDeleteData, blitz::ColumnMajorArray<2>())(t_sect(_row, _row), _col);
	}
	inline const t_part_2d operator()(const t_sect &_row, t_size _col) const {
		return t_part_2d(DATA, blitz::shape(row, col), blitz::neverDeleteData, blitz::ColumnMajorArray<2>())(_row, t_sect(_col, _col));
	}
	inline t_part_2d operator()(const t_sect &_row, t_size _col) {
		return t_part_2d(DATA, blitz::shape(row, col), blitz::neverDeleteData, blitz::ColumnMajorArray<2>())(_row, t_sect(_col, _col));
	}
	inline const t_part_2d operator()(const t_sect &_row, const t_sect &_col) const {
		return t_part_2d(DATA, blitz::shape(row, col), blitz::neverDeleteData, blitz::ColumnMajorArray<2>())(_row, _col);
	}
	inline t_part_2d operator()(const t_sect &_row, const t_sect &_col) {
		return t_part_2d(DATA, blitz::shape(row, col), blitz::neverDeleteData, blitz::ColumnMajorArray<2>())(_row, _col);
	}
	inline const t_part_1d operator()(const t_sect &_ind) const {
		return t_part_1d(DATA, num, blitz::neverDeleteData)(_ind);
	}
	inline t_part_1d operator()(const t_sect &_ind) {
		return t_part_1d(DATA, num, blitz::neverDeleteData)(_ind);
	}
	inline const t_part_1d operator()() const {
		return t_part_1d(DATA, num, blitz::neverDeleteData);
	}
	inline t_part_1d operator()() {
		return t_part_1d(DATA, num, blitz::neverDeleteData);
	}
	//Items:
	inline const TYPE &operator() (t_size _row, t_size _col) const {
		return DATA[_row + _col * row];
	}
	inline TYPE &operator() (t_size _row, t_size _col) {
		return DATA[_row + _col * row];
	}
	inline const TYPE &operator() (t_size _ind) const {
		return DATA[_ind];
	}
	inline TYPE &operator() (t_size _ind) {
		return DATA[_ind];
	}
	inline const TYPE *data() const {
		return DATA;
	}
	inline TYPE *data() {
		return DATA;
	}
	//...
	inline t_size size() const {
		return num;
	}
	inline t_size nrow() const {
		return row;
	}
	inline t_size ncol() const {
		return col;
	}
private:
	t_size num, row, col;
	TYPE *DATA;
};

//...

}

#endif //__INCLUDE_DATA_H
