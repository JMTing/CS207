/**
 * @file mtl_test.cpp
 * Test script for interfacing with MTL4 and it's linear solvers.
 */

<<<<<<< HEAD
#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

// Define a IdentityMatrix that interfaces with MTL
struct IdentityMatrix{

	IdentityMatrix(int s): s(s*s) {}

	/** Helper function to perform multiplication. Allows for delayed
	* evaluation of results and various assignment operations such * as +=, -=, and =.
	* @pre @a size(v) == size(w) */
	template <typename VectorIn , typename VectorOut , typename Assign > 
	void mult(const VectorIn& v, VectorOut& w, Assign) const{
		assert(size(v) == size(w));
		for (int i = 0; i < s; i++)
			Assign::apply(w[i], v[i]); 
	}	

	/** Matrix-vector multiplication forwards to MTL's lazy mat_cvec_multiplier*/
	template <typename VectorIn>
	mtl::vector::mat_cvec_multiplier<IdentityMatrix, VectorIn>
	operator*(const VectorIn& v) const{
		return mtl::vector::mat_cvec_multiplier<IdentityMatrix, VectorIn>(*this, v);
	}

	int s;
	private:
	// Empty!
};

/** The number of elements in the matrix */
inline std::size_t size(const IdentityMatrix& A){
	return A.s * A.s;
}

/** The number of rows in the matrix */
inline std::size_t num_rows(const IdentityMatrix& A){
	return A.s;
}	

/** The number of columns in the matrix */
inline std::size_t num_cols(const IdentityMatrix& A){
	return A.s;
}	

/** traits that MTL uses to determine properties of our IdentityMatrix. */
namespace mtl {namespace ashape{

/** Define IdentityMAtrix to be a non-scalar type */
template<>
struct ashape_aux<IdentityMatrix>{
	typedef nonscal type;
};

}	// end namespace stl

template<>
struct Collection<IdentityMatrix>{
	typedef double 		value_type;
	typedef unsigned 	size_type;
};

}	// end namespace stl

using namespace mtl;
using namespace itl;

int main(int, char**)
{
  	// Construct an IdentityMatrix and "solve" Ix = b
  	// using MTL's conjugate gradient solver

    const int size = 10, N = size * size;

    typedef IdentityMatrix matrix_type;
    matrix_type                               I(size);
    itl::pc::identity<matrix_type>            P(I);

    mtl::dense_vector<double>                 x(N, 1.0), b(N);

    b = I * x;
    x= 0;

    itl::cyclic_iteration<double>             iter(b, 100, 1.e-11, 0.0, 5);
    cg(I, x, b, P, iter);
=======
// HW3: Need to install/include Boost and MTL in Makefile
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

// HW3: YOUR CODE HERE
// Define a IdentityMatrix that interfaces with MTL

int main()
{
  // HW3: YOUR CODE HERE
  // Construct an IdentityMatrix and "solve" Ix = b
  // using MTL's conjugate gradient solver
>>>>>>> cs207/master

  return 0;
}
