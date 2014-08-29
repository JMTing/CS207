/**
 * @file poisson.cpp
 * Test script for treating the Graph as a MTL Matrix
 * and solving a Poisson equation.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles.
<<<<<<< HEAD
 * Second file: Edges (one per line) defined by 2 indices into the point list
=======
 * Second file: Eges (one per line) defined by 2 indices into the point list
>>>>>>> cs207/master
 *              of the first file.
 *
 * Launches an SDLViewer to visualize the solution.
 */

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"

#include "Graph.hpp"

#include <fstream>

<<<<<<< HEAD
#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include <math.h>

typedef Graph<int, int> GraphType;

=======
// HW3: YOUR CODE HERE
>>>>>>> cs207/master
// Define a GraphSymmetricMatrix that maps
// your Graph concept to MTL's Matrix concept. This shouldn't need to copy or
// modify Graph at all!

<<<<<<< HEAD
struct GraphSymmetricMatrix{

  // Constructor to wrap Graph
  GraphSymmetricMatrix(const GraphType& g): g_(g), s_(g.size()){}

  /** Helper function to perform multiplication. Allows for delayed
  * evaluation of results and various assignment operations such * as +=, -=, and =.
  * @pre @a size(v) == size(w) */
  template <typename VectorIn , typename VectorOut , typename Assign > 
  void mult(const VectorIn& u, VectorOut& b, Assign) const{
    assert(size(u) == size(b));
    for (uint i = 0; i < s_; i++){
      double val = 0;
      //std::cout << u[i] << '\n';

      // handles i = j 
      if (g_.node(i).value() == 1)
        val += u[i];
      else
        val += -double(g_.node(i).degree()) * u[i];

      // handles i != j
      for (auto it = g_.node(i).edge_begin(); it != g_.node(i).edge_end(); ++it){
        if ((*it).node1().value() != 1 && (*it).node2().value() != 1){
          if ((*it).node1().index() != i)
            val += u[i];
          else
            val += u[i];
        }
      }
      Assign::apply(val, b[i]);
    }
  } 

  /** Matrix-vector multiplication forwards to MTL's lazy mat_cvec_multiplier*/
  template <typename VectorIn>
  mtl::vector::mat_cvec_multiplier<GraphSymmetricMatrix, VectorIn>
  operator*(const VectorIn& v) const{
    return mtl::vector::mat_cvec_multiplier<GraphSymmetricMatrix, VectorIn>(*this, v);
  }

  GraphType g_;
  uint s_;
  private:
  // Empty!
};

/** The number of elements in the matrix */
inline std::size_t size(const GraphSymmetricMatrix& A){
  return A.s_ * A.s_;
}

/** The number of rows in the matrix */
inline std::size_t num_rows(const GraphSymmetricMatrix& A){
  return A.s_;
} 

/** The number of columns in the matrix */
inline std::size_t num_cols(const GraphSymmetricMatrix& A){
  return A.s_;
} 

/** traits that MTL uses to determine properties of our GraphSymmetricMatrix. */
namespace mtl {namespace ashape{

/** Define IdentityMAtrix to be a non-scalar type */
template<>
struct ashape_aux<GraphSymmetricMatrix>{
  typedef nonscal type;
};

} // end namespace stl

template<>
struct Collection<GraphSymmetricMatrix>{
  typedef double    value_type;
  typedef unsigned  size_type;
};

} // end namespace stl

namespace itl{

  template <class Real, class OStream = std::ostream>
  class visual_iteration : public cyclic_iteration<Real>
  {
    typedef cyclic_iteration<Real> super;
    typedef visual_iteration self;

    public:
      visual_iteration(const mtl::dense_vector<double>& r0, int max_iter_, Real tol_, Real atol_, int cycle_, 
        const mtl::dense_vector<double>& u, const GraphType& g, CS207::SDLViewer& viewer, OStream& out = std::cout)
        : super(r0, max_iter_, tol_, atol_, cycle_), u_(u), g_(g), viewer_(viewer), out(out) {
        };

      bool finished(){
        auto node_map = viewer_.empty_node_map(g_);
        viewer_.add_nodes(g_.node_begin(), g_.node_end(), NodeColor(u_), PoissonPosition(u_), node_map);
        viewer_.add_edges(g_.edge_begin(), g_.edge_end(), node_map);
        return super::finished();
      }

      template <typename T>
      bool finished(const T& r){
        auto node_map = viewer_.empty_node_map(g_);
        viewer_.add_nodes(g_.node_begin(), g_.node_end(), NodeColor(u_), PoissonPosition(u_), node_map);
        viewer_.add_edges(g_.edge_begin(), g_.edge_end(), node_map);
        bool ret= super::finished(r);
        super::print_resid();
        return ret;
      }

      mtl::dense_vector<double> u_;
      GraphType g_;
      CS207::SDLViewer viewer_;
    protected:
      OStream&   out;
  };

} // namespace itl

// checks if node is on boundary and apply boundary conditions
double boundary(Point p){
  if (norm_inf(p) == 1.0)
    return 0;
  else if (norm_inf(p-Point(-0.6,-0.6,0)) < 0.2 || norm_inf(p-Point(0.6,-0.6,0)) < 0.2 || 
    norm_inf(p-Point(-0.6,0.6,0)) < 0.2 || norm_inf(p-Point(0.6,0.6,0)) < 0.2)
    return -0.2;
  else if (p.x > -0.6 && p.x < 0.6 && p.y > -0.2 && p.y < 0.2 && p.z > -1.0 && p.z < 1.0)
    return 1;
  // indicates node is not on the boundary
  else
    return -1;
}

// forcing function
double force(Point p){
  return 5 * cos(norm_1(p));
}

// color functor
template <typename dv>
struct NodeColor {
  dv u_;
  NodeColor(mtl::dense_vector<double> u): u_(u){
  };
  template <typename NODE>
  CS207::Color operator()(const NODE& n) {
    return CS207::Color::make_heat(static_cast<float> (u_[n.index()]));
  }
};

// node position functor
struct PoissonPosition {
  mtl::dense_vector<double> u_;
  PoissonPosition(mtl::dense_vector<double> u): u_(u){
  };
  template <typename NODE>
  Point operator()(const NODE& n) {
    return Point(n.position().x, n.position().y, u_[n.index()]);
  }
};

using namespace mtl;
using namespace itl;

=======
>>>>>>> cs207/master
int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE EDGES_FILE\n";
    exit(1);
  }

<<<<<<< HEAD
  GraphType graph;
  std::vector<GraphType::node_type> nodes;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  while (CS207::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,2> t;
  while (CS207::getline_parsed(tets_file, t))
    graph.add_edge(nodes[t[0]], nodes[t[1]]);

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  viewer.launch();

  // Construct a GraphSymmetricMatrix and "solve" Ix = b
  // using MTL's conjugate gradient solver
  const uint size = graph.size();

  typedef GraphSymmetricMatrix matrix_type;
  matrix_type                               A(graph);

  dense_vector<double>                 u(size, 1.0), b(size);

  // defines b and flags boundaries with values as 1 and non-boundaries as 0
  // all edge length are equal and treated as h = 1
  for (uint i = 0; i < size; i++){
    if (boundary(graph.node(i).position()) != -1){
      b[i] = boundary(graph.node(i).position());
      graph.node(i).value() = 1;
    }

    else{
      bool hf = true;
      for (auto it = graph.node(i).edge_begin(); it != graph.node(i).edge_end(); ++it){
        if (hf) {
            b[i] = (*it).length() * force(graph.node(i).position());
            hf = false;
          }
        if ((*it).node1().index() != i && boundary((*it).node1().position()) != -1)
          b[i] -= boundary((*it).node1().position());
        else if ((*it).node2().index() != i && boundary((*it).node2().position()) != -1)
          b[i] -= boundary((*it).node2().position());
      }
      graph.node(i).value() = 0;
    }
  }
  u = 0;

  cyclic_iteration<double>             iter(b, 175, 1.e-10, 0.0, 50);
  visual_iteration<double>             iter2(b, 175, 1.e-10, 0.0, 50, u, graph, viewer);
  cg(A, u, b, iter);
=======
  // HW3: YOUR CODE HERE
  // Construct the Graph
  // Construct the GraphSymmetricMatrix
  // Define b
  // Solve Au = b using MTL
>>>>>>> cs207/master

  return 0;
}
