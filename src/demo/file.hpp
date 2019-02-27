#ifndef __INCLUDE_FILE_H
#define __INCLUDE_FILE_H

#include <calib/calib.hpp>

#include <cstdio>

namespace CALIB {

template <typename TYPE> struct t_file {

	inline void write(const t_hand<const TYPE> &_ptr) {

		if (numt == 0) {
			data_write(_ptr->grid()->flow()->rect().minx()); data_write(_ptr->grid()->flow()->rect().miny());
			data_write(_ptr->grid()->flow()->rect().maxx()); data_write(_ptr->grid()->flow()->rect().maxy());
			head = ftell(fid); data_write(numt); data_write(_ptr);
			curr = ftell(fid);
			fseek(fid, head, SEEK_SET);
		}
		else {
			fseek(fid, curr, SEEK_SET); data_write(_ptr);
			curr = ftell(fid);
			fseek(fid, head, SEEK_SET); data_write(numt);
		}
		fflush(fid);
		++ numt;
	}

	inline t_file(FILE *&&_fid): fid(_fid) {}

	inline ~t_file() { fclose(fid); }

private:
	template <typename T>
	inline void data_write(const T *_dat, size_t _num) { fwrite(_dat, sizeof(T) * _num, 1, fid); }
	template <typename T>
	inline void data_write(T _val) { fwrite(&_val, sizeof(T), 1, fid); }
	//...
	void data_write(const t_hand<const TYPE> &_ptr);
	void data_read(const t_hand<const TYPE> &_ptr);
	//...
	long long head, curr;
	t_size numt = 0;
	FILE *fid;
};

template <> void t_file<t_scal_isoline>::data_write(const c_scal_isoline &SCAL) {

	const t_grid::t_node &NODE = SCAL->grid()->node();
	const t_grid::t_cont &CONT = SCAL->grid()->cont();
	const t_grid::t_step &STEP = SCAL->grid()->step();

	t_data<t_real, 1> DATA; t_data<t_size, 1> SIZE; t_data<t_byte, 1> STAT;

	t_size num_k = NODE.size(), num_c = CONT.size(), num_s = STEP.size();

	DATA.resize(num_k);
	data_write(num_k);
	DATA() = NODE.valx(); data_write(DATA.data(), num_k);
	DATA() = NODE.valy(); data_write(DATA.data(), num_k);
	//...
	SIZE.resize(num_c); STAT.resize(num_c);
	data_write(num_c);
	SIZE() = CONT.head(); data_write(SIZE.data(), num_c);
	SIZE() = CONT.tail(); data_write(SIZE.data(), num_c);
	STAT() = CONT.stat(); data_write(STAT.data(), num_c);
	//...
	SIZE.resize(num_s);
	data_write(num_s);
	SIZE() = STEP.head(); data_write(SIZE.data(), num_s);
	SIZE() = STEP.tail(); data_write(SIZE.data(), num_s);
}

template <> void t_file<t_scal_regular>::data_write(const c_scal_regular &SCAL) {

	if (numt == 0) {
		data_write(SCAL->grid()->rect().numx());
		data_write(SCAL->grid()->rect().numy());
	}
	data_write(
	SCAL->valz().data(), SCAL->valz().size()
	);
}

typedef t_file<t_scal_isoline> f_scal_isoline;
typedef t_file<t_scal_regular> f_scal_regular;

//...

}

#endif //__INCLUDE_FILE_H
