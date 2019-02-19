#ifndef __INCLUDE_HAND_H
#define __INCLUDE_HAND_H

#include <type_traits>
#include <cstddef>
#include <memory>
#include <set>

using std::nullptr_t;

namespace CALIB {

/** Объявление шаблонного класса для linked_ptr **/

struct t_buff {
	inline explicit t_buff(void *_dat): data(_dat) {}
	std::set<struct t_link *> hard;	//Список сильных ссылок;
	std::set<struct t_link *> weak;	//Список слабых ссылок;
	void *data;
};

struct t_link {
	inline void mov(t_buff *_new) {
		if (_new == nullptr) { del(); return; }
		t_buff *_buf = buff;
		if (_buf == nullptr) {
			set(_new);
			return;
		}
		for (auto &_tmp = _buf->weak; _tmp.size() != 0; (*_tmp.begin())->set(_new));
		for (auto &_tmp = _buf->hard; _tmp.size() != 1; (*_tmp.begin())->set(_new));
		(
		*(_buf->hard).begin()
		)->set(_new);
	}
	inline void del() {
		if (buff == nullptr) return;
		t_buff *_buf = buff;
		for (auto &_tmp = _buf->weak; _tmp.size() != 0; (*_tmp.begin())->off());
		for (auto &_tmp = _buf->hard; _tmp.size() != 1; (*_tmp.begin())->off());
		(
		*(_buf->hard).begin()
		)->off();
	}
	virtual void set(t_buff *_buf) = 0;
	virtual void off() = 0;
protected:
	t_buff *buff = nullptr;
};

template <typename R, typename T>
using CONVERT = typename std::enable_if<std::is_convertible<R *, T *>::value>::type;

template <typename T>
struct t_hand;

template <typename T>
struct t_weak;

template <typename T>
struct t_hand:
private
t_link {

	inline operator  const  T *() const { return (buff != nullptr)? static_cast<T *> (buff->data): (nullptr); }
	inline const T *operator->() const { return static_cast<T *> (buff->data); }
	inline const T &operator*() const { return *static_cast<T *> (buff->data); }
	inline T *operator->() { return static_cast<T *> (buff->data); }
	inline T &operator*() { return *static_cast<T *> (buff->data); }

	template <typename> friend class t_hand;
	template <typename> friend class t_weak;

	template <typename R, typename = CONVERT<R, T> >
	inline t_hand &operator=(const t_hand<R> &_rhs) {
		set(_rhs.buff);
		return *this;
	}
	template <typename R, typename = CONVERT<R, T> >
	inline t_hand &operator=(const t_weak<R> &_rhs) {
		set(_rhs.buff);
		return *this;
	}
	inline t_hand &operator=(const t_hand &_rhs) {
		set(_rhs.buff);
		return *this;
	}
	template <typename R, typename = CONVERT<R, T> >
	inline t_hand &operator=(R *&&_dat) {
		set(new t_buff(_dat));
		return *this;
	}
	inline t_hand &operator=(nullptr_t) {
		off();
		return *this;
	}

	template <typename H>
	inline t_hand<H> get() const {
		static_assert(!std::is_const<T>::value ||
		               std::is_const<H>::value, "Violation of a const qualifier!");
		if (!buff || !dynamic_cast<H *> (static_cast<T *> (buff->data))) return nullptr;
		t_hand<H> ptr;
		ptr.set(buff);
		return ptr;
	}

	inline t_hand<T> out() {
		static_assert(!std::is_const<T>::value, "Violation of a const qualifier!");
		t_hand<T> ptr;
		if (buff) { ptr.set(new t_buff(buff->data)); buff->data = nullptr; }
		del();
		return ptr;
	}

	template <typename R, typename = CONVERT<R, T> >
	inline void mov(const t_hand<R> &_rhs) {
		static_assert(!std::is_const<T>::value &&
		              !std::is_const<R>::value, "Violation of a const qualifier!");
		t_link::mov(_rhs.buff);
	}
	template <typename R, typename = CONVERT<R, T> >
	inline void mov(const t_weak<R> &_rhs) {
		static_assert(!std::is_const<T>::value &&
		              !std::is_const<R>::value, "Violation of a const qualifier!");
		t_link::mov(_rhs.buff);
	}
	template <typename R, typename = CONVERT<R, T> >
	inline void mov(R *&&_dat) {
		static_assert(!std::is_const<T>::value &&
		              !std::is_const<R>::value, "Violation of a const qualifier!");
		mov(t_hand<R>(std::move(_dat)));
	}

	template <typename R>
	inline void set(const R &_rhs) { (*this) = _rhs; }
	template <typename R>
	inline void set(R &_rhs) { (*this) = _rhs; }
	template <typename R>
	inline void set(R &&_rhs) {
		(*this) = std::move(_rhs);
	}

	inline void del() {
		static_assert(!std::is_const<T>::value, "Violation of a const qualifier!");
		t_link::del();
	}
	inline void off() {
		if (buff == nullptr) return;
		//assert(buff->hard.size());
		buff->hard.erase(this);
		if (buff->hard.size() == 0) {
			for (auto &_tmp = buff->weak; _tmp.size(); (*_tmp.begin())->off());
			delete
			static_cast<T *> (
				buff->data
			);
			delete buff;
		}
		buff = nullptr;
	}

	template <typename R, typename = CONVERT<R, T> >
	inline t_hand(const t_hand<R> &_rhs): t_hand() { (*this) = _rhs; }
	template <typename R, typename = CONVERT<R, T> >
	inline t_hand(const t_weak<R> &_rhs): t_hand() { (*this) = _rhs; }

	inline t_hand(const t_hand &_rhs): t_hand() { (*this) = _rhs; }

	template <typename R, typename = CONVERT<R, T> >
	inline t_hand(R *&&_dat): t_hand() {
		(*this) = std::move(_dat);
	}

	inline t_hand(nullptr_t): t_link() {}

	inline t_hand(): t_link() {}

	virtual ~t_hand() {
		off();
	}
private:
	void set(t_buff *_buf) {
		if (buff != _buf) { off(); buff = _buf; }
		if (buff) {
			buff->hard.insert(this);
		}
	}
};

template <typename T>
struct t_weak:
private
t_link {

	inline operator  const  T *() const { return (buff != nullptr)? static_cast<T *> (buff->data): (nullptr); }
	inline const T *operator->() const { return static_cast<T *> (buff->data); }
	inline const T &operator*() const { return *static_cast<T *> (buff->data); }
	inline T *operator->() { return static_cast<T *> (buff->data); }
	inline T &operator*() { return *static_cast<T *> (buff->data); }

	template <typename> friend class t_hand;
	template <typename> friend class t_weak;

	template <typename R, typename = CONVERT<R, T> >
	inline t_weak &operator=(const t_hand<R> &_rhs) {
		set(_rhs.buff);
		return *this;
	}
	template <typename R, typename = CONVERT<R, T> >
	inline t_weak &operator=(const t_weak<R> &_rhs) {
		set(_rhs.buff);
		return *this;
	}
	inline t_weak &operator=(const t_weak &_rhs) {
		set(_rhs.buff);
		return *this;
	}
	inline t_weak &operator=(nullptr_t) {
		off();
		return *this;
	}

	template <typename H>
	inline t_weak<H> get() const {
		static_assert(!std::is_const<T>::value ||
		               std::is_const<H>::value, "Violation of a const qualifier!");
		t_weak<H> ptr;
		if (!buff || !dynamic_cast<H *> (static_cast<T *> (buff->data))) return ptr;
		ptr.set(buff);
		return ptr;
	}

	inline t_hand<T> out() {
		static_assert(!std::is_const<T>::value, "Violation of a const qualifier!");
		t_hand<T> ptr;
		if (buff) { ptr.set(new t_buff(buff->data)); buff->data = nullptr; }
		del();
		return ptr;
	}

	template <typename R, typename = CONVERT<R, T> >
	inline void mov(const t_hand<R> &_rhs) {
		static_assert(!std::is_const<T>::value &&
		              !std::is_const<R>::value, "Violation of a const qualifier!");
		t_link::mov(_rhs.buff);
	}
	template <typename R, typename = CONVERT<R, T> >
	inline void mov(const t_weak<R> &_rhs) {
		static_assert(!std::is_const<T>::value &&
		              !std::is_const<R>::value, "Violation of a const qualifier!");
		t_link::mov(_rhs.buff);
	}
	template <typename R, typename = CONVERT<R, T> >
	inline void mov(R *&&_dat) {
		static_assert(!std::is_const<T>::value &&
		              !std::is_const<R>::value, "Violation of a const qualifier!");
		mov(t_hand<R>(std::move(_dat)));
	}

	template <typename R>
	inline void set(const R &_rhs) { (*this) = _rhs; }
	template <typename R>
	inline void set(R &_rhs) { (*this) = _rhs; }
	template <typename R>
	inline void set(R &&_rhs) {
		(*this) = std::move(_rhs);
	}

	inline void del() {
		static_assert(!std::is_const<T>::value, "Violation of a const qualifier!");
		t_link::del();
	}
	inline void off() {
		if (buff != nullptr) { buff->weak.erase(this); buff = nullptr; }
	}

	template <typename R, typename = CONVERT<R, T> >
	inline t_weak(const t_hand<R> &_rhs): t_weak() { (*this) = _rhs; }
	template <typename R, typename = CONVERT<R, T> >
	inline t_weak(const t_weak<R> &_rhs): t_weak() { (*this) = _rhs; }

	inline t_weak(const t_weak &_rhs): t_weak() { (*this) = _rhs; }

	inline t_weak(nullptr_t): t_link() {}

	inline t_weak(): t_link() {}

	virtual ~t_weak() {
		off();
	}
private:
	void set(t_buff *_buf) {
		if (buff != _buf) { off(); buff = _buf; }
		if (buff) {
			buff->weak.insert(this);
		}
	}
};

//...

}

#endif //__INCLUDE_HAND_H
