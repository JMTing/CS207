#include <iostream>
#include <vector>
#include <cassert>

// USE_EXPR_TEMP == 0  --  Default to the trivial op* and op+
// USE_EXPR_TEMP == 1  --  Use lazy evaluation to prevent temporaries
#define USE_EXPR_TEMP 1

/** Curiously Recurring Template Pattern (CRTP) for Matrix Expressions */
template <typename E>
struct MatExpr {
  const E& derived() const {
    return static_cast<const E&>(*this);
  }
  double operator()(unsigned i, unsigned j) const {
    return derived()(i,j);
  }
  unsigned rows() const {
    return derived().rows();
  }
  unsigned cols() const {
    return derived().cols();
  }
  // No member data! Only interface.
};

/** Lazy evaluation for a scaled matrix expression */
template <typename E1>
struct MatScale : public MatExpr<MatScale<E1>> {
  MatScale(double b, const MatExpr<E1>& a)
      : alpha_(b), a_(a.derived()) {
  }
  MatScale(const MatExpr<E1>& a, double b)
      : alpha_(b), a_(a.derived()) {
  }
  unsigned rows() const {
    return a_.rows();
  }
  unsigned cols() const {
    return a_.cols();
  }
  double operator()(unsigned i, unsigned j) const {
    return alpha_ * a_(i,j);
  }
private:
  double alpha_;
  const E1& a_;
};

/** Lazy evaluation for the sum of two matrix expressions */
template <typename E1, typename E2>
struct MatAdd : public MatExpr<MatAdd<E1, E2>> {
  MatAdd(const MatExpr<E1>& a, const MatExpr<E2>& b)
      : a_(a.derived()), b_(b.derived()) {
    assert(a_.rows() == b.rows() && a.cols() == b.cols());
  }
  unsigned rows() const {
    return a_.rows();
  }
  unsigned cols() const {
    return b_.cols();
  }
  double operator()(unsigned i, unsigned j) const {
    return a_(i, j) + b_(i, j);
  }
private:
  const E1& a_;
  const E2& b_;
};

/** Lazy evaluation for the product of two matrix expressions */
template <typename E1, typename E2>
struct MatMult : public MatExpr<MatMult<E1, E2>> {
  MatMult(const MatExpr<E1>& a, const MatExpr<E2>& b)
      : a_(a.derived()), b_(b.derived()) {
    assert(a_.cols() == b_.rows());
  }
  unsigned rows() const {
    return a_.rows();
  }
  unsigned cols() const {
    return b_.cols();
  }
  double operator()(unsigned i, unsigned j) const {
    double result = 0;
    for (unsigned k = 0; k < a_.cols(); ++k)
      result += a_(i, k) * b_(k, j);
    return result;
  }
private:
  const E1& a_;
  const E2& b_;
};

/** Concrete Matrix of values */
struct Matrix : public MatExpr<Matrix> {
  Matrix(unsigned rows, unsigned cols, double val = 0)
      : rows_(rows), cols_(cols), v_(rows_*cols_, val) {
  }
  template <typename E>
  Matrix(const MatExpr<E>& m)
      : rows_(m.rows()), cols_(m.cols()) {
    v_.reserve(rows() * cols());
    for (unsigned i = 0; i < rows(); ++i)
      for (unsigned j = 0; j < cols(); ++j)
        v_.push_back(m(i, j));
  }
  unsigned rows() const {
    return rows_;
  }
  unsigned cols() const {
    return cols_;
  }
  const double& operator()(unsigned i, unsigned j) const {
    return v_[i*rows() + j];
  }
  double& operator()(unsigned i, unsigned j) {
    return v_[i*rows() + j];
  }
private:
  unsigned rows_, cols_;
  std::vector<double> v_;
};


/** Output the values of any matrix expression */
template <typename E>
std::ostream& operator<<(std::ostream& os, const MatExpr<E>& m) {
  for (unsigned i = 0; i < m.rows(); ++i) {
    for (unsigned j = 0; j < m.cols(); ++j)
      os << m(i, j) << "\t";
    os << "\n";
  }
  return os;
}


#if USE_EXPR_TEMP

/** Add any two matrix expressions */
template <typename E1, typename E2>
MatAdd<E1, E2> operator+(const MatExpr<E1>& a, const MatExpr<E2>& b) {
   return MatAdd<E1, E2>(a, b);
}

/** Multiply any two matrix expressions */
template <typename E1, typename E2>
MatMult<E1, E2> operator*(const MatExpr<E1>& a, const MatExpr<E2>& b) {
  return MatMult<E1, E2>(a, b);
}

/** Scale any matrix expression */
template <typename E>
MatScale<E> operator*(const MatExpr<E>& a, double b) {
  return MatScale<E>(a, b);
}
template <typename E>
MatScale<E> operator*(double b, const MatExpr<E>& a) {
  return MatScale<E>(b, a);
}

#else

// Classic op+ for two Matrices
Matrix operator+(const Matrix& a, const Matrix& b) {
  Matrix r(a.rows(), b.cols());
  for (unsigned i = 0; i < r.rows(); ++i)
    for (unsigned j = 0; j < r.cols(); ++j)
      r(i,j) = a(i,j) + b(i,j);
  return r;
}

// Classic op* for two Matrices
Matrix operator*(const Matrix& a, const Matrix& b) {
  Matrix r(a.rows(), b.cols());
  for (unsigned i = 0; i < r.rows(); ++i)
    for (unsigned j = 0; j < r.cols(); ++j)
      for (unsigned k = 0; k < a.cols(); ++k)
	r(i,j) += a(i,k) * b(k,j);
  return r;
}

// Classic op* for a scalar and a Matrix
Matrix operator*(const Matrix& a, double b) {
  Matrix r(a.rows(), a.cols());
  for (unsigned i = 0; i < r.rows(); ++i)
    for (unsigned j = 0; j < r.cols(); ++j)
      r(i,j) = b * a(i,j);
  return r;
}
Matrix operator*(double b, const Matrix& a) {
  Matrix r(a.rows(), a.cols());
  for (unsigned i = 0; i < r.rows(); ++i)
    for (unsigned j = 0; j < r.cols(); ++j)
      r(i,j) = b * a(i,j);
  return r;
}

#endif


#include <chrono>

class Clock {
 public:
  typedef std::chrono::high_resolution_clock clock;
  typedef typename clock::time_point         time_point;
  typedef typename clock::duration           duration;
  typedef typename duration::rep             tick_type;
  // Construct the Clock
  Clock() {
    start();
  }
  // Start the Clock
  void start() {
    start_ = clock::now();
  }
  // Get the duration (in seconds) on this clock
  double elapsed() const {
    typedef std::chrono::duration<double> units;
    return std::chrono::duration_cast<units>(clock::now() - start_).count();
  }
 private:
  time_point start_;
};

int main() {
  unsigned N = 512;
  Matrix A(N, N, 1);
  Matrix B(N, N, 2);
  Matrix C(N, N, 5.76);
  double alpha = 3.14;

  Clock c;

  c.start();
  Matrix D = A*B + alpha*B + alpha*A + C;
  double time = c.elapsed();
  std::cout << "D(0,0) = " << D(0,0)
            << " computed in " << time << " secs"<< std::endl;
}
