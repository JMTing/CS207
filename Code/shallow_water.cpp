/**
 * @file shallow_water.cpp
 * Implementation of a shallow water system using Mesh
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D point list (one per line) defined by three doubles
 * Second file: Triangles (one per line) defined by 3 indices into the point list
 */

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"
#include "Point.hpp"
#include <fstream>
#include <math.h> 

#include "Mesh.hpp"

// Standard gravity (average gravity at Earth's surface) in meters/sec^2
static constexpr double grav = 9.80665;

/** Water column characteristics */
struct QVar {
  double h;	  // Height of fluid
  double hu;	// Height times average x velocity of column
  double hv;	// Height times average y velocity of column

  /** Default constructor.
   *
   * A default water column is 1 unit high with no velocity. */
  QVar()
    : h(1), hu(0), hv(0) {
  }

  /** Construct the given water column. */
  QVar(double h_, double hu_, double hv_)
    : h(h_), hu(hu_), hv(hv_) {
  }

  /** Add QVar @a b to this QVar */
  QVar& operator+=(const QVar b) {
    h += b.h;
    hu += b.hu;
    hv += b.hv;
    return *this;
  }

  /** Subtract QVar @a b from this QVar */
  QVar& operator-=(const QVar b) {
    h -= b.h;
    hu -= b.hu;
    hv -= b.hv;
    return *this;
  }

  /** Scale this QVar up by scalar @a b */
  QVar& operator*=(double b) {
    h *= b;
    hu *= b;
    hv *= b;
    return *this;
  }

  /** Divide this QVar by scalar @a b */
  QVar& operator/=(double b) {
    h /= b;
    hu /= b;
    hv /= b;
    return *this;
  }

};

/** Add QVar @a a to @a b QVar */
QVar operator+(QVar a, const QVar& b) {
  return a += b;
}

/** Subtract QVar @a a from @a b QVar */
QVar operator-(QVar a, const QVar& b) {
  return a -= b;
}

/** Multiply QVar @a a by @a b scalar */
QVar operator*(QVar a, double b) {
  return a *= b;
}

/** Multiply QVar @a a by @a b scalar */
QVar operator*(double b, QVar a) {
  return a *= b;
}

/** Divide QVar @a a by @a b QVar */
QVar operator/(QVar a, double b) {
  return a /= b;
}

/** Write a QVar to an output stream */
std::ostream& operator<<(std::ostream& s, const QVar& a) {
  return (s << a.h << ' ' << a.hu << ' ' << a.hv);
}

/** Water column characteristics */
struct floatObj {
  double F;   // Force in the z-direction 
  double A;   // Area of 
  double ro;  // fluid density

  /** Default constructor.
   *
   * A default floating Obj has no force on water */
  floatObj()
    : F(0), A(1), ro(1000) {
  }

  /** Construct the given water column. */
  floatObj(double F_, double A_, double ro_)
    : F(F_), A(A_), ro(ro_) {
  }
};

/** Function object for calculating shallow-water flux.
 *          |e
 *   T_k    |---> n = (nx,ny)   T_m
 *   QBar_k |                   QBar_m
 *          |
 * @param[in] nx, ny Defines the 2D outward normal vector n = (@a nx, @a ny)
 *            from triangle T_k to triangle T_m. The length of n is equal to the
 *            the length of the edge, |n| = |e|.
 * @param[in] dt The time step taken by the simulation. Used to compute the
 *               Lax-Wendroff dissipation term.
 * @param[in] qk The values of the conserved variables on the left of the edge.
 * @param[in] qm The values of the conserved variables on the right of the edge.
 * @return The flux of the conserved values across the edge e
 */
struct EdgeFluxCalculator {
  QVar operator()(double nx, double ny, double dt,
                  const QVar& qk, const QVar& qm, floatObj obj = floatObj()) {
    double e_length = sqrt(nx*nx + ny*ny);
    nx /= e_length;
    ny /= e_length;

    // The velocities normal to the edge
    double wm = (qm.hu*nx + qm.hv*ny) / qm.h;
    double wk = (qk.hu*nx + qk.hv*ny) / qk.h;

    // Lax-Wendroff local dissipation coefficient
    double vm = sqrt(grav*qm.h) + sqrt(qm.hu*qm.hu + qm.hv*qm.hv) / qm.h;
    double vk = sqrt(grav*qk.h) + sqrt(qk.hu*qk.hu + qk.hv*qk.hv) / qk.h;
    double a  = dt * std::max(vm*vm, vk*vk);

    // Helper values which accounts for floating objects
    double scale = 0.5 * e_length;
    double gh2   = (0.5 * grav + obj.F/(obj.A * obj.ro)) * (qm.h*qm.h + qk.h*qk.h);

    // Simple flux with dissipation for stability
    return QVar(scale * (wm*qm.h  + wk*qk.h)           - a * (qm.h  - qk.h),
                scale * (wm*qm.hu + wk*qk.hu + gh2*nx) - a * (qm.hu - qk.hu),
                scale * (wm*qm.hv + wk*qk.hv + gh2*ny) - a * (qm.hv - qk.hv));
  }
};

// defines an arbitrary force for floating objects-default values for now
struct ForceCalculator{
  floatObj operator()(){
    return floatObj();
  }
};

/** Node position function object for use in the SDLViewer. */
struct NodePosition {
  template <typename NODE>
  Point operator()(const NODE& n) {
    return Point(n.position().x, n.position().y, n.value().q.h);
  }
};

// Value types for the mesh
struct NodeData{
  QVar q;
  double b;
};

struct EdgeData{
  float flux;
};

struct TriData{
  QVar q_bar;
  QVar q_bar2;
  QVar S;
};

typedef Mesh<NodeData,EdgeData,TriData> MeshType;

/** Integrate a hyperbolic conservation law defined over the mesh m
 *  with flux functor f by dt in time.
 */
template <typename MESH, typename FLUX, typename FORCE>
double hyperbolic_step(MESH& m, FLUX& f, FORCE& force, double t, double dt) {
  for (auto tri_it = m.triangle_begin(); tri_it != m.triangle_end(); ++tri_it) {
    QVar sum = QVar(0,0,0);

    // Defines the force of the floating objects
    floatObj obj = force();

    Point norm_vec = (*tri_it).normal((*tri_it).node1(), (*tri_it).node2());
    // Calcuate flux using adjacent triangle's Q
    if ((*tri_it).adj_triangle((*tri_it).node1(), (*tri_it).node2()).is_valid())
      sum += f(norm_vec.x, norm_vec.y, dt, (*tri_it).value().q_bar, 
      		(*tri_it).adj_triangle((*tri_it).node1(), (*tri_it).node2()).value().q_bar, obj);
    // Boundary conditions
    else
      sum += f(norm_vec.x, norm_vec.y, dt, (*tri_it).value().q_bar, 
      		QVar((*tri_it).value().q_bar.h, 0, 0), obj);

    norm_vec = (*tri_it).normal((*tri_it).node2(), (*tri_it).node3());
    // Calcuate flux using adjacent triangle's Q
    if ((*tri_it).adj_triangle((*tri_it).node2(), (*tri_it).node3()).is_valid()) 
      sum += f(norm_vec.x, norm_vec.y, dt, (*tri_it).value().q_bar, 
      		(*tri_it).adj_triangle((*tri_it).node2(), (*tri_it).node3()).value().q_bar, obj); 
    // Boundary conditions
    else 
      sum += f(norm_vec.x, norm_vec.y, dt, (*tri_it).value().q_bar, 
      		QVar((*tri_it).value().q_bar.h, 0, 0), obj);

    norm_vec = (*tri_it).normal((*tri_it).node3(), (*tri_it).node1());
    // Calcuate flux using adjacent triangle's Q
    if ((*tri_it).adj_triangle((*tri_it).node3(), (*tri_it).node1()).is_valid())
      sum += f(norm_vec.x, norm_vec.y, dt, (*tri_it).value().q_bar, 
      		(*tri_it).adj_triangle((*tri_it).node3(), (*tri_it).node1()).value().q_bar, obj); 
    // Boundary conditions
    else
      sum += f(norm_vec.x, norm_vec.y, dt, (*tri_it).value().q_bar, 
      		QVar((*tri_it).value().q_bar.h, 0, 0), obj);

    (*tri_it).value().q_bar2 = ((*tri_it).value().S * dt/(*tri_it).area()) + (*tri_it).value().q_bar - 
                               (sum * dt/(*tri_it).area());
  }

  // In order to use original (as opposed to update) q_bar values in the equation above, we have to update
  // q_bar values for all triangles after the q_bar values have been calculated.
  for (auto tri_it = m.triangle_begin(); tri_it != m.triangle_end(); ++tri_it)
    (*tri_it).value().q_bar = (*tri_it).value().q_bar2;

  return t + dt;
}

/** Convert the triangle-averaged values to node-averaged values for viewing. */
template <typename MESH>
void post_process(MESH& m) {
  for (auto n_it = m.node_begin(); n_it != m.node_end(); ++n_it) {
    double sum_area = 0;
    QVar value = QVar(0,0,0);

    for (auto tri_it = m.adj_triangle_begin((*n_it).uid()); tri_it != m.adj_triangle_end((*n_it).uid()); ++tri_it) {
      value += (*tri_it).value().q_bar * (*tri_it).area();
      sum_area += (*tri_it).area();
    }

    (*n_it).value().q = value * 1.0/(sum_area);
  }
}

/* Initial conditions functor for still water*/
struct Still{
  QVar operator()(Point p) {
    return QVar(1.0,0,0);
  }
};

/* Initial conditions functor*/
struct Pebble{
  QVar operator()(Point p) {
    return QVar(1.0 - 0.75*exp(-80*(pow((p.x-0.75),2) + pow(p.y,2))), 0, 0);
  }
};

/* Initial conditions functor*/
struct Column {
  QVar operator()(Point p) {
    if (pow(p.x-0.75,2) + pow(p.y,2) - pow(0.15,2) < 0)
      return QVar(1.75,0,0);
    else
      return QVar(1.0,0,0);
  }
};

/* Initial conditions functor*/
struct DamBreak {
  QVar operator()(Point p) {
    if (p.x < 0)
      return QVar(1.75,0,0);
    else
      return QVar(1.0,0,0);
  }
};

/* Initial conditions functor for bathymetry-default 0 for now */
struct bathyMetry {
  double operator()(Point p, double factor) {
    return p.x/(factor*10);
    // return 0;
  }
};

int main(int argc, char* argv[]) {
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: shallow_water NODES_FILE TRIS_FILE\n";
    exit(1);
  }

  MeshType mesh;
  std::vector<typename MeshType::node_type> mesh_node;

  // Read all Points and add them to the Mesh
  std::ifstream nodes_file(argv[1]);
  Point p;
  while (CS207::getline_parsed(nodes_file, p)) {
    mesh_node.push_back(mesh.add_node(p));
  }

  // Read all mesh triangles and add them to the Mesh
  std::ifstream tris_file(argv[2]);
  std::array<int,3> t;
  while (CS207::getline_parsed(tris_file, t)) {
    mesh.add_triangle(mesh_node[t[0]], mesh_node[t[1]], mesh_node[t[2]]);
  }

  // Print out the stats
  std::cout << mesh.num_nodes() << " "
            << mesh.num_edges() << " "
            << mesh.num_triangles() << std::endl;


  /* Set the initial conditions */ 
  // Set the initial values of the nodes and get the maximum height double
  double max_h = 0;
  double dx = 0;
  double dy = 0;
  auto init_cond = Still();
  auto b_init_cond = bathyMetry(); 

  // Find the maximum height and apply initial conditions to nodes
  for (auto it = mesh.node_begin(); it != mesh.node_end(); ++it) { 
    auto n = *it;
    n.value().q = init_cond(n.position());
    n.value().b = b_init_cond(n.position(), mesh.num_nodes());
    max_h = std::max(max_h, n.value().q.h);
  } 

  // Set the initial values of the triangles to the average of their nodes and finds S 
  for (auto it = mesh.triangle_begin(); it != mesh.triangle_end(); ++it) {
    auto t = *it; 
    t.value().q_bar = (t.node1().value().q + 
                       t.node2().value().q + 
                       t.node3().value().q) / 3.0;
    t.value().q_bar2 = t.value().q_bar;

    double b_avg = (t.node1().value().b + 
                    t.node2().value().b + 
                    t.node3().value().b) / 3.0;
    // finds the max dx and dy to calculate Source
    dx = std::max(dx, fabs(t.node1().position().x - t.node2().position().x));
    dx = std::max(dx, fabs(t.node2().position().x - t.node3().position().x));
    dx = std::max(dx, fabs(t.node3().position().x - t.node1().position().x));
    dy = std::max(dy, fabs(t.node1().position().y - t.node2().position().y));
    dy = std::max(dy, fabs(t.node2().position().y - t.node3().position().y));
    dy = std::max(dy, fabs(t.node3().position().y - t.node1().position().y));
    t.value().S = QVar(0, -grav * t.value().q_bar.h * b_avg / dx, -grav * t.value().q_bar.h * b_avg / dy);
  }

  // Calculate the minimum edge length
  double min_edge_length = std::numeric_limits<double>::max();
  for (auto it = mesh.edge_begin(); it != mesh.edge_end(); ++it) { 
    min_edge_length = std::min(min_edge_length, (*it).length());
  } 
	
  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  viewer.launch();

  auto node_map = viewer.empty_node_map(mesh);
  viewer.add_nodes(mesh.node_begin(), mesh.node_end(),
                   CS207::DefaultColor(), NodePosition(), node_map);
  viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);
  viewer.center_view();

  // CFL stability condition requires dt <= dx / max|velocity|
  // For the shallow water equations with u = v = 0 initial conditions
  //   we can compute the minimum edge length and maximum original water height
  //   to set the time-step
  // Compute the minimum edge length and maximum water height for computing dt
  double dt = 0.25 * min_edge_length / (sqrt(grav * max_h));
  double t_start = 0;
  double t_end = 10;

  // Preconstruct a Flux functor
  EdgeFluxCalculator f;
  // Preconstruct a Force functor
  ForceCalculator Force;
  // Begin the time stepping
  for (double t = t_start; t < t_end; t += dt) {
    // Step forward in time with forward Euler
    hyperbolic_step(mesh, f, Force, t, dt);

    // Update node values with triangle-averaged values
    post_process(mesh);

    // Update the viewer with new node positions
    viewer.add_nodes(mesh.node_begin(), mesh.node_end(),
                     CS207::DefaultColor(), NodePosition(), node_map);
    viewer.set_label(t);

    // These lines slow down the animation for small meshes.
    if (mesh.num_nodes() < 100)
      CS207::sleep(0.05);
  }
  return 0;
}
