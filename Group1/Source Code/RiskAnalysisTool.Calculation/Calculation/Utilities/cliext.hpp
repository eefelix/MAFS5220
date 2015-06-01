#pragma once

#include <utility>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

namespace cliext {
#ifdef _M_CEE

	template <typename T, template <class> class NativeImpl = boost::shared_ptr >
	public ref class shared_ptr {
	public:
		shared_ptr() {
			impl_ = new NativeImpl < T > ;
		}

		shared_ptr(const NativeImpl<T> &ptr) {
			impl_ = new NativeImpl<T>(ptr);
		}

		shared_ptr(const shared_ptr<T> %ptr) {
			impl_ = new NativeImpl<T>(*(ptr.impl_));
		}

		shared_ptr(T *ptr) {
			impl_ = new NativeImpl<T>(ptr);
		}

		~shared_ptr() {
			(*this).!shared_ptr();
		}

		!shared_ptr() {
			delete impl_;
		}

		const NativeImpl<T> &native() {
			return *impl_;
		}

		void reset() {
			delete impl_;
			impl_ = new NativeImpl<T>();
		}

		void reset(const NativeImpl<T> &other) {
			delete impl_;
			impl_ = new NativeImpl<T>(other);
		}


	public:
		T* operator ->() {
			return (*impl_).get();;
		}

		shared_ptr<T> % operator = (const NativeImpl<T> &other) {
			delete impl_;
			impl_ = new NativeImpl<T>(other);
			return *this;
		}

		shared_ptr<T> % operator = (const shared_ptr<T> %other) {
			delete impl_;
			impl_ = new NativeImpl<T>(*(other.impl_));
			return *this;
		}

	private:
		NativeImpl<T> *impl_;
	};


	template <typename T, typename ... Args>
	shared_ptr<T> make_shared(Args && ... args) {
		return shared_ptr<T>(boost::make_shared<T>(::std::forward(args) ...));
	}

#endif
}

