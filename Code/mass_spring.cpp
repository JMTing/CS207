/**
 * @file mass_spring.cpp
 * Implementation of mass-spring system using Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <fstream>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"

#include "Graph.hpp"
#include "Point.hpp"


// Gravity in meters/sec^2
static constexpr double grav = 9.81;

// holds the node's mass and velocity in its value as the graph type
struct node_data{
  double m_;
  Point v_;
  node_data(double m, const Point& v): m_(m), v_(v){
  };
};

// holds the edge's spring constant and resting length in its value as the graph type
struct edge_data{
  double k_;
  double l_;
  edge_data(double k, double l): k_(k), l_(l){
  };
};

typedef Graph<node_data, edge_data> GraphType;
typedef GraphType::Node Node;
typedef GraphType::Edge Edge;

/* Defines the gravitational force on the node, which is defined as mi*(0,0,-g)
** Complexity O(1)
*/
struct GravityForce{
  template <typename NODE>
  Point operator()(NODE n){
    Point f_g = {0, 0, -grav};
    return n.value().m_*f_g;
  }
};

/* Defines the spring force, which is defined as SUM(-k*(xi-xj)/|xi-xj|*(|xi-xj|-L))
 * K is the spring constant, L is the spring rest length, and |.| is the euclidean distance
 * Complexity O(d), where d is the number of adjacent nodes
*/ 
struct MassSpringForce{
  template <typename NODE>
  Point operator()(NODE n){
    Point f_s = {0, 0, 0};
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it)
      f_s += -(*it).value().k_ * (((*it).node1().position() - (*it).node2().position()) /
          (norm((*it).node1().position() - (*it).node2().position()))) * 
          (norm((*it).node1().position() - (*it).node2().position()) - (*it).length());
    return f_s;
  }
};

/* Defines the damping force, which is defined as -c*vi
 * vi is the velocity of the node i and c is the damping constant set by user
 * Complexity O(1), where d is the number of adjacent nodes
*/ 
struct DampingForce{
  double c_;
  DampingForce(double c):c_(c){
  };
  template <typename NODE>
  Point operator()(NODE n){
    return -c_ * n.value().v_;
  }
};

/* Defines a constraint for the nodes by a user defined z-plane
 * Fixes nodes that violates this constrain by setting the position to the 
 *  nearest point on the plane and setting the z-component of the Node v to 0
 * Complexity O(1)
*/ 
struct z_plane_constraint{
  double p_;
  z_plane_constraint(double p):p_(p){
  };
  template <typename NODE>
  void operator()(NODE n){
    Point constraint = {0, 0, 1};
    if (dot(n.position(), constraint) < p_){
      Point close = {n.position().x, n.position().y, p_};
      n.set_position(close);
      n.value() = node_data(n.value().m_, 0);
    }
  }
};

/* Defines a constraint for the nodes by a user defined sphere with center and radius
 * Fixes nodes that violates this constraint by setting the position to the 
 *  nearest point on the surface of the sphere and setting the component of
 *  the velocity that is normal to the sphere's surface to 0
 * Complexity O(1)
*/ 
struct sphere_constraint{
  Point c_;
  double r_;
  sphere_constraint(Point c, double r):c_(c), r_(r){
  };
  template <typename NODE>
  void operator()(NODE n){
    if (norm(n.position()-c_) < r_){
      // finds the unit vector for fixing nodes that violates constraints
      Point R = (n.position()-c_)/(norm(n.position()-c_));
      n.set_position(R * r_);
      R = n.value().v_ - dot(n.value().v_, R)*R;
      n.value() = node_data(n.value().m_, R);
    }
  }
};

/* Defines a constraint for the nodes that falls within it's shortest edge length r
 * Fixes nodes that violates this constraint by reseting the positon and velocity
 * stores 
 * Complexity O(d)
*/ 
struct self_collision{
  Point original_pos_;
  Point original_vel_;
  double r_;
  // stores original values and find min edge length
  template <typename NODE>
  self_collision(NODE n){
    original_pos_ = n.position();
    original_vel_ = n.value().v_;
    r_ = std::numeric_limits<double>::infinity();
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
        if ((*it).value().l_ < r_)
          r_ = (*it).value().l_;
    }
  }

  template <typename NODE>
  void operator()(NODE n){
      for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
        if ((*it).value().l_ < r_){
          n.set_position(original_pos_);
          n.value() = node_data(n.value().m_, original_vel_);
        }
      }
    }
 };


/* Combines 2 constraints
 * @pre constraints are already defined
*/ 
template <typename NODE, typename C1, typename C2>
void make_combined_constraint(NODE n, C1 c1, C2 c2){
    c1(n);
    c2(n); 
}


/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g Graph
 * @param[in] t The current time (useful for time-dependent forces)
 * @param[in] dt The time step
 * @param[in] force Function object defining the force per node
 * @pre G::node_value_type supports node_data struct
 * @return the next time step (usually @a t + @a dt)
 *
 * @a force is called as @a force(n, @a t), where n is a node of the graph
 * and @a t is the current time parameter. @a force must return a Point
 * representing the node's force at time @a t.
 * Time complexity with the no collision constrain: O(N*d)
 */
template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
  // point for sphere
  //Point a = {0.5, 0.5, -0.5};
  for (auto it = g.node_begin(); it != g.node_end(); ++it){
    // stores the original values for the self-collision contraint
    self_collision c(*it);

    // updates the position
    Point x_n1 = (*it).position() + (*it).value().v_ * dt;
    (*it).set_position(x_n1);

    // finds the force for updating the velocity
    GravityForce f_g;
    MassSpringForce f_s;
    DampingForce f_d(1.0/g.num_nodes());
    Point v_n1 = (*it).value().v_ + (force(f_g(*it), f_s(*it), f_d(*it))*dt)/(*it).value().m_;
    (*it).value() = node_data((*it).value().m_, v_n1);

    // calls the constraint
    c(*it);
    //z_plane_constraint c1(-0.75);
    //sphere_constraint c2(a, 0.15);
    //make_combined_constraint((*it), c1, c2);
  }
  return t + dt;
}

/* Combines 2-3 forces
 * @pre forces are already calculated
*/ 
struct make_combined_force{
  Point operator()(Point f1, Point f2, const Point& f3 = Point({0,0,0})){
    return f1 + f2 + f3;
  }
};

/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n) {
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      return Point (0 ,0 ,0);

    // gravitational force
    Point f_g = {0, 0, -grav};

    // spring forces
    Point f_s = {0, 0, 0};
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
      f_s += -(*it).value().k_ * (((*it).node1().position() - (*it).node2().position()) /
          (norm((*it).node1().position() - (*it).node2().position()))) * 
          (norm((*it).node1().position() - (*it).node2().position()) - (*it).length());
      }

    return n.value().m_*f_g + f_s;
  }
};


int main(int argc, char* argv[]) {
  // Check arguments
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  GraphType graph;
  std::vector<typename GraphType::node_type> nodes;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  while (CS207::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CS207::getline_parsed(tets_file, t)) {
    for (unsigned i = 1; i < t.size(); ++i) {
      graph.add_edge(nodes[t[0]], nodes[t[1]]);
      graph.add_edge(nodes[t[0]], nodes[t[2]]);

      // Diagonal edges: include as of HW2 #2
      graph.add_edge(nodes[t[0]], nodes[t[3]]);
      graph.add_edge(nodes[t[1]], nodes[t[2]]);

      graph.add_edge(nodes[t[1]], nodes[t[3]]);
      graph.add_edge(nodes[t[2]], nodes[t[3]]);
    }
  }

  /* Set initial conditions for your nodes
   * Sets mass to 1/N where N is the number of nodes in the graph (constant density)
   * Sets L to initial length of the edges
   * Sets K to 100
  */
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it)
    (*it).value().m_ = 1.0/graph.num_nodes();

  for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it){
    (*it).value().l_ = (*it).length();
    (*it).value().k_ = 100.0;
  }

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  auto node_map = viewer.empty_node_map(graph);
  viewer.launch();

  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  viewer.center_view();

  // Begin the mass-spring simulation
  double dt = 0.0001;
  double t_start = 0;
  double t_end = 5.0;

  for (double t = t_start; t < t_end; t += dt) {
    //std::cout << "t = " << t << std::endl;
    symp_euler_step(graph, t, dt, make_combined_force());
    // Update viewer with nodes' new positions
    viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
    viewer.set_label(t);

    // These lines slow down the animation for small graphs, like grid0_*.
    // Feel free to remove them or tweak the constants.
    if (graph.size() < 100)
      CS207::sleep(0.001);
  }

  return 0;
}
