#include <iostream>
#include <utility>
#include <type_traits>
#include <cstring>
#include <cassert>

/** Generic copy for any type */
template <typename T>
T* copy(T* first, T* last, T* out, std::false_type) {
  for ( ; first != last; ++first, ++out)
    *out = *first;
  return out;
}

/** Copy optimized for types with trivial copy */
template <typename T>
T* copy(T* first, T* last, T* out, std::true_type) {
  memmove(out, first, (last-first)*sizeof(T));
  return out + (first - last);
}


template <typename T, typename use_trivial_copy = std::is_trivial<T> >
class Vector {
 public:
  typedef T value_type;
  typedef unsigned size_type;

  Vector()
    : v_(), n_(0), cap_(0) {
  }
  explicit Vector(size_type n)
    : v_(new T[n]), n_(n), cap_(n) {
  }
  ~Vector() {
    delete[] v_;
  }
  size_type size() const {
    return n_;
  }
  void push_back(const T& x) {
    if (n_ >= cap_) {
      cap_ = (cap_ == 0 ? 16 : 2*cap_);
      T* new_v = new T[cap_];
      copy(v_, v_ + n_, new_v, use_trivial_copy());
      delete[] v_;
      v_ = new_v;
    }
    v_[n_] = x;
    ++n_;
  }
  const T& operator[](size_type i) const {
    return v_[i];
  }

 private:
  T* v_;
  size_t n_;
  size_t cap_;
};



int main()
{
  struct my_type {
    unsigned long long x[2];
  };

  Vector<my_type> v(3);
  for (unsigned long long i = 100; i < 5000; ++i) {
    v.push_back(my_type{{i,i}});
  }

  std::cout << v[v.size()-1].x[0] << std::endl;
}
