/**
 * Jason Ting
 * Computer Science 207 Final Project
 *
 * @file final_project.cpp
 * Implementation of a ball falling into shallow water system using mesh
 * User can pause the simulation by right clicking
 *
 * @brief Reads in four files specified on the command line.
 * First file: 3D point list (one per line) defined by three doubles for water
 * Second file: Triangles (one per line) defined by 3 indices into the point list for water
 * Third file: 3D point list (one per line) defined by three doubles for ball
 * Fourth file: Triangles (one per line) defined by 3 indices into the point list for ball
 */


#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"
#include "Point.hpp"
#include <fstream>
#include <math.h>
#define _USE_MATH_DEFINES 

#include "Mesh.hpp"

// Standard gravity (average gravity at Earth's surface) in meters/sec^2
static constexpr double grav = 9.80665;

// Total mass of the ball mesh
static constexpr double total_mass = 375.0;

// Scale the volume of ball so it fits shallow water
static constexpr double scale = 0.20;

// interaction factor for wind force
static constexpr double wind_const = 0.00009;

// wind velocity constant
static const Point wind_velocity_const = Point(2.0, 3.0, 0.0);

// spring constant
static constexpr double spring_const = 35.0;

// spring constant
static constexpr double gas_const = 30;

// Height the center of the ball starts at
static constexpr double start_height = 2.0;

// plane position constant
static constexpr double plane_const = 0.5;

// Density of the fluids
static constexpr double density = 1000;

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

/** Floating object characteristics */
struct floatObj {
  double F;   // Force in the z-direction 
  double A;   // Area of 
  double ro;  // Fluid density

  /** Default constructor.
   *
   * A default floating Obj has no force on water */
  floatObj()
    : F(0), A(1), ro(density) {
  }

  /** Construct the given floating objext. */
  floatObj(double F_, double A_, double ro_)
    : F(F_), A(A_), ro(ro_) {
  }
};

// Value types for the mesh
struct NodeData{
  QVar q;
  double b;
  double mass;
  Point velocity;
};

struct EdgeData{
  double spring_constant;
  double initial_length;
};

struct TriData{
  QVar q_bar;
  QVar q_bar2;
  QVar S;
  // multiplier holds a -1 or 1, this will be set before the simulation
  // so that every time normal vectors to triangle surfaces are 
  // calculated, we multiply by this, so that we always get the 
  // normal vector that is pointing to the outside of the shape.
  double multiplier;
};

typedef Mesh<NodeData,EdgeData,TriData> MeshType;
typedef MeshType::Node Node;
typedef MeshType::Triangle Triangle;
typedef MeshType::size_type size_type;


/** Calculates radius of cross sectional radius of ball submerged
*  @pre 0 <= dh <= 2 * radius
*/
double cross_radius(double radius, double dh){
  assert(dh <= 2 * radius);
  return sqrt(radius*radius - (radius-dh)*(radius-dh));
}

/** Calculate the point at the center of a mesh. 
 *  Calculated center = sum(node positions)/num_nodes() 
 *  @param[in,out] m Mesh
 *  @return a point that is the calculated center of the mesh
 */
template <typename MESH>
Point get_center(MESH& m)
{
  Point center = Point(0,0,0);
  for(auto it=m.node_begin(); it != m.node_end(); ++ it)
    center += (*it).position();
  center /= m.num_nodes();

  return center;
}

/** Given a triangle, set the triangle's multiplier value
 *  to -1 or 1, depending on whether the calculated normal
 *  vector is pointing towards the center of the mesh object
 *  or away from it.  
 *  @pre The mesh models a convex shape
 *  @param[in,out] i Triangle
 *  @param[in] center Point at the center of the mesh
 *  @post i.value().multiplier = 1 if the calculated
 *  normal vector points away from the inside of the modeled
 *  shape and -1 if it points towards the inside
 */
void set_normal_direction(Triangle i, Point center){
  i.value().multiplier = 1.0;
  Point e01 = i.node3().position() - i.node1().position();
  Point e02 = i.node2().position() - i.node1().position();
  
  Point vec_normal = cross(e01, e02);
  Point vec_center = i.node1().position() - center;
  if (dot(vec_normal, vec_center) < 0)
    i.value().multiplier = -1.0;
}

/** Find the normal vector to a triangle that points towards
 *  the outside of the original shape
 *  @param[in] i Triangle
 *  @returns A normal vector pointing away from the "inside" of
 *  the original shape
 */
Point get_normal_surface(Triangle& i){
  Point e01 = i.node2().position() - i.node1().position();
  Point e02 = i.node3().position() - i.node1().position();
  
  Point vec_normal = cross(e01, e02);
  
  vec_normal *= i.value().multiplier;
  vec_normal /= norm(vec_normal);
  return vec_normal;
}

// Force function object for assessing the force due to gravity on the nodes of the mesh. 
struct GravityForce {
  /** Return the force due to gravity applying to @a n. 
  *  @pre n is templated with the node class from mesh.hpp
  *  @param[in] m The mesh of containing the nodes
  *  @param[in] n The node whose forces we are calculating
  *  @param[in] t The time in the simulation. This is not actually used
  *  here, but we take it as a parameter for consistency since other
  *  forces may use time as a parameter.
  *  @post Returns a force due to gravity, where 
  *  F=Point(0,0,n.value().mass*-grav)
  */
  template <typename MESH, typename NODE>
  Point operator()(MESH& m, NODE n, double t) {
    (void) t;
    (void) m;
    return Point(0,0,n.value().mass*-grav);
  }
};

/** Force function object for assessing the wind force on 
 *  the surfaces of the mesh. */
struct WindForce {
 private:
   Point wind_velocity_;  
 public:
   // public constructor: given wind velocity
   WindForce(const Point& wind_velocity) : wind_velocity_(wind_velocity) {}
   // public constructor: use the default wind velocity
   WindForce() : wind_velocity_(wind_velocity_const) {}

   // Wind Force function on the node based on the wind force
   // on each triangle surface surrounding the node
   template <typename MESH, typename NODE>
   Point operator() (MESH& m, NODE n, double t){
     Point node_velocity = n.value().velocity;
     Point node_normal(0.0, 0.0, 0.0);
     for (auto it = m.adj_triangle_begin(n.uid()); it != m.adj_triangle_end(n.uid()); ++it){
       auto tri = (*it);
       // approximate node normal vector by the sum of normal vectors of its neighboring faces
       node_normal +=  get_normal_surface(tri);
     }
     // calculte wind force
     Point force = wind_const * dot(wind_velocity_ - node_velocity, node_normal) * node_normal;
     (void) t;
     return force;
   }
};

/** Force function object for assessing the pressure force on 
 *  the interior surfaces of the mesh. */
struct PressureForce {
  // constructors
  PressureForce(){
  }

  PressureForce(double p) : pressure(p) {}
  
  // Pressure force function returns the pressure on each node
  // based on the calculated pressures on all of that node's 
  // surrounding triangle faces
  template <typename MESH, typename NODE>
  Point operator()(MESH& m, NODE n, double t) {
    (void) t;
    Point pressure_force = Point(0.0, 0.0, 0.0);

    for(auto it=m.adj_triangle_begin(n.uid()); it != m.adj_triangle_end(n.uid()); ++ it) {
      auto tri = (*it);
      Point normal = get_normal_surface(tri);
      double area = tri.area();
      pressure_force += normal * area * pressure;

    }
    return pressure_force; 
  }  

  void set_pressure(double p) {
    pressure = p;
  }

 private:
  double pressure;
};

/** Force function object for assessing the force due to buoyancy */
struct BuoyantForce {
  /** Return the force due to springs applying to @a n. 
  *  @pre n is templated with the node class from mesh.hpp
  *  @param[in] m The mesh of containing the nodes
  *  @param[in] n The node whose forces we are calculating
  *  @param[in] t The time in the simulation. This is not actually used
  *  here, but we take it as a parameter for consistency since other
  *  forces may use time as a parameter.
  *  @post Returns a force due to amount of sphere submerged
  *  Force = (0,0,PI*dh_/6*(3*a_*a_+dh_*dh_)
  */
  // constructors
  BuoyantForce(): dh_(0), a_(0){}

  BuoyantForce(double dh, double a) : dh_(dh), a_(a) {}

  template <typename MESH, typename NODE>
  Point operator()(MESH& m, NODE n, double t) {
    (void) t;
    (void) m;
    (void) n;
    return Point(0,0,M_PI*dh_/6*(3*a_*a_+dh_*dh_));
  }
  private:
    double dh_;
    double a_;
};

/** Force function object for assessing the force due to springs on the nodes of the mesh. */
struct MassSpringForce {
  /** Return the force due to springs applying to @a n. 
  *  @pre n is templated with the node class from mesh.hpp
  *  @pre the Mesh class in mesh.hpp contains an edge incident iterator for
  *  edges incident to each node
  *  @param[in] n The node whose forces we are calculating
  *  @param[in] t The time in the simulation. This is not actually used
  *  here, but we take it as a parameter for consistency since other
  *  forces may use time as a parameter.
  *  @post Returns a force due to springs
  */

  template <typename MESH, typename NODE>
  Point operator()(MESH& m, NODE n, double t) {
    (void) t;
    (void) m;
    Point force = Point(0,0,0);
    Node node2;
    double distance; 
    double initial_spring_length;
    double spring_const_this_edge;
    
    // iterate through each neighboring node and add spring forces
    for(auto it = n.edge_begin(); it != n.edge_end(); ++it){
        node2 = (*it).node2();
        
        spring_const_this_edge = (*it).value().spring_constant; 
        initial_spring_length  = (*it).value().initial_length;
        distance = norm(n.position()-node2.position()); // Euclidean distance
        force += (-1)*spring_const_this_edge*(n.position()-node2.position())*(distance-initial_spring_length)/distance;
    }
    
    return force;
  }
};

/** Force function object for assessing the damping force on the nodes of the mesh. */
struct DampingForce {
  /** Return the force due to dampening applying to @a n. Dampening 
   *  force is defined as force=-c*velocity where c is a constant equal to 
   *  1/(num_nodes) 
   *  @pre n is templated with the node class from mesh.hpp
   *  @param[in] n The node whose forces we are calculating
   *  @param[in] t The time in the simulation. This is not actually used
   *  here, but we take it as a parameter for consistency since other
   *  forces may use time as a parameter.
   *  @post Returns a force due to damping where 
   *  F=n.value().velocity/num_nodes
   */
  // constructors
  DampingForce(){
  }

  DampingForce(int num)
  {
      num_nodes = num;
  }
  template <typename MESH, typename NODE>
  Point operator()(MESH& m, NODE n, double t) {
    (void) t;
    (void) m;
    Point velocity = n.value().velocity;
    Point return_force = -velocity/num_nodes;
    return Point(return_force.x,return_force.y,0);
  }
 private:
  int num_nodes; //holds the number of nodes in the mesh
};

/** This structure is templated so that it can take two force structs
 *  such as DampingForce, GravityForce or MassSpringForce above and add
 *  them together. When the () operator is called, this structure will
 *  return the two templated forces.*/
template <typename F1, typename F2>
struct combined_force {
  // constuctors  
  combined_force() {}

  combined_force(F1& force1, F2& force2) : f1(force1), f2(force2) {}
  /** returns the combination of the two templated forces. 
   *  @pre F1 and F2 are both structs that contain an overloaded 
   *  () operator with a templated first parameter, a double for the
   *  second parameter and an int for the third parameter. This 
   *  operator on F1 and F2 must return a force_type value.
   *  @pre n is templated with the node class from mesh.hpp
   *  @param[in] n The node whose forces we are calculating
   *  @param[in] t The time in the simulation
   *  @post Returns a force which is the result of combining the forces
   *  from the passed in structures @a F1 and @a F2 
   */
  template <typename MESH, typename NODE>
  Point operator()(MESH& m, NODE n, double t) {
      return f1(m,n,t)+f2(m,n,t);
  }

 private:
  F1 f1;
  F2 f2;
};

/** returns the combination of the two templated forces. 
 *  @pre F1 and F2 are both struct types that contain an overloaded 
 *  () operator with a templated first parameter, a double for the
 *  second parameter and an int for the third parameter. This 
 *  operator on F1 and F2 must return a force_type value. 
 *  @param[in] f1 The first force structure
 *  @param[in] f2 The second force structure
 *  @post Returns a force structure which is the result of combining the 
 *  forces passed in
 */
template<typename F1, typename F2>
combined_force<F1,F2> make_combined_force(F1& f1, F2& f2) {
    auto return_value = combined_force<F1,F2>(f1,f2);
    return return_value;
}

/** returns the combination of the three templated forces. 
 *  @pre F1, F2 and F3 are all struct types that contain an overloaded 
 *  () operator with a templated first parameter, a double for the
 *  second parameter and an int for the third parameter. This 
 *  operator on F1, F2 and F3 must return a force_type value. 
 *  @param[in] f1 The first force structure
 *  @param[in] f2 The second force structure
 *  @param[in] f3 The second force structure
 *  @post Returns a force structure which is the result of combining the 
 *  forces passed in
 */
template<typename F1, typename F2, typename F3>
combined_force<combined_force<F1,F2>,F3> make_combined_force(F1& f1, F2& f2, F3& f3) {
    combined_force<F1,F2> two_combined = make_combined_force(f1,f2);
    auto three_combined = make_combined_force(two_combined,f3);
    return three_combined;
}

/** returns the combination of the four templated forces. 
 *  @pre F1, F2, F3 and F4 are all struct types that contain an overloaded 
 *  () operator with a templated first parameter, a double for the
 *  second parameter and an int for the third parameter. This 
 *  operator on F1, F2, F3 and F4 must return a force_type value. 
 *  @param[in] f1 The first force structure
 *  @param[in] f2 The second force structure
 *  @param[in] f3 The third force structure
 *  @param[in] f4 The fourth force structure
 *  @post Returns a force structure which is the result of combining the 
 *  forces passed in
 */
template<typename F1, typename F2, typename F3, typename F4>
combined_force<combined_force<F1,F2>,combined_force<F3,F4>> make_combined_force(F1& f1, F2& f2, F3& f3, F4& f4) {
    combined_force<F1,F2> two_combined1 = make_combined_force(f1,f2);
    combined_force<F3,F4> two_combined2 = make_combined_force(f3,f4);
    auto four_combined = make_combined_force(two_combined1,two_combined2);
    return four_combined;
}

/** Constraint function object for constraining the nodes of the mesh
 *  so they cannot go below -0 on the z plane. 
 */
struct PlaneConstraint {
 private:
   double plane_;
 public:
   // public constructor: given plane position
   PlaneConstraint(const double& plane) : plane_(plane) {}
   // public constructor: use the default plane position
   PlaneConstraint() : plane_(0) {}
  /** Set these nodes' velocities on to 0.
   *  @pre m is templated with the Mesh class from mesh.hpp
   *  @param[in] m The mesh whose nodes we are constraining
   *  @param[in] t The time in the simulation. This is not actually used
   *  here, but we take it as a parameter for consistency since other
   *  constraint structures may use time as a parameter.
   *  @post Nodes that were below -0.75 on the z plane now have velocities of 0.
   */  
  template <typename MESH>
  void operator()(MESH& m, double z_bottom) {
    if (z_bottom < plane_)
      for(auto it=m.node_begin(); it != m.node_end(); ++it)
        (*it).value().velocity = Point(0,0,0);
  }
};

// Defines an arbitrary force for floating objects-default value
struct ForceObj{
  floatObj operator()(){
    return floatObj();
  }
};

// Defines an arbitrary force for floating objects on to the surface of a triangle
// Suggested forces includes courage, power, and wisdom =D
struct TriForce{
  template <typename TRIANGLE>
  floatObj operator()(const TRIANGLE& tri){
    return floatObj(1,tri.area(),density);
  }
};

/** Node position function object for use in the SDLViewer. */
struct NodePosition {
  template <typename NODE>
  Point operator()(const NODE& n) {
    return Point(n.position().x, n.position().y, n.value().q.h);
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
 * @param[in] obj The object floating on the water
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

    // Helper values which accounts for floating objects (MODIFIED EQN HERE)
    double scale = 0.5 * e_length;
    double gh2   = (0.5 * grav) * (qm.h*qm.h + qk.h*qk.h) + (obj.F/(obj.A * obj.ro)*(qm.h+qk.h));

    // Simple flux with dissipation for stability
    return QVar(scale * (wm*qm.h  + wk*qk.h)           - a * (qm.h  - qk.h),
                scale * (wm*qm.hu + wk*qk.hu + gh2*nx) - a * (qm.hu - qk.hu),
                scale * (wm*qm.hv + wk*qk.hv + gh2*ny) - a * (qm.hv - qk.hv));
  }
};

/** Integrate a hyperbolic conservation law defined over the mesh m
 *  with flux functor f by dt in time.
 */
template <typename MESH, typename FLUX>
double hyperbolic_step(MESH& m, FLUX& f, double t, double dt, Point ball_loc, int water_tris) {
  for (auto tri_it = m.triangle_begin(); (*tri_it).index() != water_tris; ++tri_it) {
    QVar sum = QVar(0,0,0);

    /* Defines the force of the floating objects by checking its submerged 
    *  cross section of the ball intersects with the centers of adjacent triangles 
    */
    floatObj obj = floatObj();
    if (norm((*tri_it).center()-Point(ball_loc.x, ball_loc.y, 0)) < ball_loc.z)
      obj = floatObj(total_mass*grav, M_PI*ball_loc.z*ball_loc.z, density);

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

  /* In order to use original (as opposed to update) q_bar values in the equation above, we have to update
   * q_bar values for all triangles after the q_bar values have been calculated.
   * Updates the Source term
   */
  for (auto tri_it = m.triangle_begin(); (*tri_it).index() != water_tris; ++tri_it){
    (*tri_it).value().S = (*tri_it).value().S/(*tri_it).value().q_bar.h*(*tri_it).value().q_bar2.h;
    (*tri_it).value().q_bar = (*tri_it).value().q_bar2;
  }

  return t + dt;
}

/** Convert the triangle-averaged values to node-averaged values for viewing and 
  * apply simple euler step to move the ball with forces
  */
template <typename MESH, typename FORCE, typename CONSTRAINT>
Point post_process(MESH& m, FORCE force, CONSTRAINT& c, double t, double dt, uint water_nodes) {
  static double ball_bottom = std::numeric_limits<double>::max();
  double water_dis = std::numeric_limits<double>::max();
  double dh = 0;
  static double submerged_height = 0;
  Point bottom_loc;
  Point water_loc;

  for (auto n_it = m.node_begin(); n_it != m.node_end(); ++n_it) {
    // handles the shallow water
  	if ((*n_it).index() < water_nodes) {
	    double sum_area = 0;
	    QVar value = QVar(0,0,0);

	    for (auto tri_it = m.adj_triangle_begin((*n_it).uid()); tri_it != m.adj_triangle_end((*n_it).uid()); ++tri_it) {
	      value += (*tri_it).value().q_bar * (*tri_it).area();
	      sum_area += (*tri_it).area();
	    }

	    (*n_it).value().q = value * 1.0/(sum_area);
	  }
    // handles the ball
    else {
      (*n_it).set_position((*n_it).position() + (*n_it).value().velocity*dt);
      dh = (*n_it).value().q.h - (*n_it).position().z;
      (*n_it).value().q.h = (*n_it).position().z;
      (*n_it).value().velocity += force(m,(*n_it),t)*dt/(*n_it).value().mass;
      if ((*n_it).value().q.h < ball_bottom) {
        ball_bottom = (*n_it).position().z;
        bottom_loc = (*n_it).position();
      }
    }
  }
  // find the water node closest to the bottom of the ball
  for (auto n_it = m.node_begin(); (*n_it).index() != water_nodes; ++n_it) {
    if (norm(Point((*n_it).position().x,(*n_it).position().y,0)-Point(bottom_loc.x,bottom_loc.y,0)) < water_dis){
      water_dis = norm(Point((*n_it).position().x,(*n_it).position().y,0)-Point(bottom_loc.x,bottom_loc.y,0));
      water_loc = Point((*n_it).position().x,(*n_it).position().y,(*n_it).value().q.h);
    }
  }

  // apply contraints of neccessary
  c(m,ball_bottom);

  // determines if the ball fell below shallow water and updates height submerged
  if (bottom_loc.z < water_loc.z)
    submerged_height += dh;
  return Point(bottom_loc.x, bottom_loc.y, submerged_height);
}

/* Initial conditions functor for still water*/
struct Still{
  QVar operator()(Point p) {
    (void) p;
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

/* Initial conditions functor for bathymetry Source 0*/
struct StillSource {
  double operator()(Point p) {
    (void) p;
    return 0;
  }
};

/* Initial conditions functor for bathymetry Source*/
struct SlopeX {
  double operator()(Point p) {
    return fabs(p.x/100000);
  }
};

/* Initial conditions functor for bathymetry Source*/
struct SlopeY {
  double operator()(Point p) {
    return fabs(p.y/100000);
  }
};

/* Initial conditions functor for bathymetry Source*/
struct Cone {
  double operator()(Point p) {
    return sqrt(p.x*p.x + p.y * p.y)/100000;
  }
};

/* Initial conditions functor for bathymetry Source*/
struct ReverseCone {
  double operator()(Point p) {
    return -sqrt(p.x*p.x + p.y * p.y)/100000;
  }
};

// SDL listener- pauses and unpauses the simulation on right clicks
struct my_listener : public CS207::SDL_Listener {
  void handle(SDL_Event e) {   // we are forced to implement this function
    switch (e.type) {
      case SDL_MOUSEBUTTONDOWN: {
        if (e.button.button == SDL_BUTTON_RIGHT ) {
          if (dt_ == original_dt_)
            dt_ = 0;
          else
            dt_ = original_dt_;
        }
      } break;
    }
  }
  
  // constructor
  my_listener(CS207::SDLViewer& viewer, MeshType& mesh, double& dt) 
    : viewer_(viewer), mesh_(mesh), dt_(dt), original_dt_(dt) {};
  private:
   CS207::SDLViewer& viewer_;
   MeshType& mesh_;
   double& dt_;
   double original_dt_;
};

// calls the color map
struct Color {
  uint water_nodes_;
  Color(uint water_nodes) : water_nodes_(water_nodes) {}
  CS207::Color operator()(const MeshType::node_type n) {
    // handles shallow water by making it blue
    if (n.index() < water_nodes_){ 
      float value = (n.value().q.h-0.5)*1.5;
      if (value > 1)
        value = 1;
      else if (value < 0)
        value = 0;
      return CS207::Color::make_hsv(0.666,value,value);
    }
    // handles the ball
    else {
      // change to teal if submerged, otherwise green
      if (n.value().q.h < 1)
        return CS207::Color::make_hsv(0.5,0.75,0.33);
      return CS207::Color::make_hsv(0.333,0.75,0.50);
    }
  }
};





int main(int argc, char* argv[]) {
  // Check arguments
  if (argc < 5) {
    std::cerr << "Usage: final_project NODES_FILE TRIS_FILE ball.nodes ball.tris \n";
    exit(1);
  }

  MeshType mesh;
  std::vector<typename MeshType::node_type> mesh_node;

  // Read all water Points and add them to the Mesh
  std::ifstream nodes_file(argv[1]);
  Point p;
  uint water_nodes = 0;
  while (CS207::getline_parsed(nodes_file, p)) {
    mesh_node.push_back(mesh.add_node(p));
    water_nodes++;
  }

  // Read all water mesh triangles and add them to the Mesh
  std::ifstream tris_file(argv[2]);
  std::array<int,3> t;
  int water_tris = 0;
  while (CS207::getline_parsed(tris_file, t)) {
    mesh.add_triangle(mesh_node[t[0]], mesh_node[t[1]], mesh_node[t[2]]);
    water_tris++;
  }
  uint water_edges = mesh.num_edges();

  std::ifstream nodes_file2(argv[3]);
  double radius = 1 * scale;
  while (CS207::getline_parsed(nodes_file2, p)) {
    p *= scale;
    p.z += + start_height;
    mesh_node.push_back(mesh.add_node(p));
  }

  // Read all ball mesh triangles and add them to the mesh
  std::ifstream tris_file2(argv[4]);
  while (CS207::getline_parsed(tris_file2, t)) {
    mesh.add_triangle(mesh_node[t[0]+water_nodes], mesh_node[t[1]+water_nodes], mesh_node[t[2]+water_nodes]);
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
  auto b_init_cond = Cone(); 

  // Find the maximum height and apply initial conditions to nodes
  for (auto it = mesh.node_begin(); it != mesh.node_end(); ++it) { 
    auto n = *it;
    if (n.index() < water_nodes){
	    n.value().q = init_cond(n.position());
	    n.value().b = b_init_cond(n.position());
	    max_h = std::max(max_h, n.value().q.h);
  	}
  	else {
  		n.value().q = QVar(n.position().z, 0, 0);
      n.value().mass = total_mass/(mesh.num_nodes() - water_nodes);
      n.value().velocity = Point(0.0,0.0,0.0);
  	}
  } 

  // Set the initial values of the triangles to the average of their nodes and finds S
  // Set the triangle direction values so that we can determine which 
  // way to point normal vectors. This part assumes a convex shape
  Point center = get_center(mesh); 
  for (auto it = mesh.triangle_begin(); it != mesh.triangle_end(); ++it) {
    auto t = *it; 
    if (t.index() < water_tris){
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
    else
      set_normal_direction((*it),center);
  }

  // Calculate the minimum edge length and set edge inital condititons
  double min_edge_length = std::numeric_limits<double>::max();
  uint count = 0;
  for (auto it = mesh.edge_begin(); it != mesh.edge_end(); ++it, count++){
    if (count < water_edges)
      min_edge_length = std::min(min_edge_length, (*it).length());
    else {
      (*it).value().spring_constant = spring_const;
      (*it).value().initial_length = (*it).length();
    }
  }
	
  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  viewer.launch();

  auto node_map = viewer.empty_node_map(mesh);
  viewer.add_nodes(mesh.node_begin(), mesh.node_end(),
                   Color(water_nodes), NodePosition(), node_map);
  viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);
  // adds solid color-slows down program significantly
  //viewer.add_triangles(mesh.triangle_begin(), mesh.triangle_end(), node_map);
  viewer.center_view();


  // CFL stability condition requires dt <= dx / max|velocity|
  // For the shallow water equations with u = v = 0 initial conditions
  //   we can compute the minimum edge length and maximum original water height
  //   to set the time-step
  // Compute the minimum edge length and maximum water height for computing dt
  double dt = 0.25 * min_edge_length / (sqrt(grav * max_h));
  double t_start = 0;
  double t_end = 10;
  Point ball_loc = Point(0,0,0);
  double dh = 0;
  // double pressure = gas_const/(4/3*M_PI*radius*radius*radius);

  // add listener
  my_listener* l = new my_listener(viewer,mesh,dt); 
  viewer.add_listener(l);

  // Preconstruct a Flux functor
  EdgeFluxCalculator f;
  // defines the constraint
  PlaneConstraint c = PlaneConstraint(plane_const);

  // Begin the time stepping
  for (double t = t_start; t < t_end; t += dt) {
    // define forces on ball
    GravityForce g_force;
    BuoyantForce b_force = BuoyantForce(dh, ball_loc.z);
    WindForce w_force;
    MassSpringForce ms_force;
    // PressureForce p_force = PressureForce(pressure);
    // DampingForce d_force = DampingForce(mesh.num_nodes());
    auto combined_forces = make_combined_force(g_force, b_force, w_force, ms_force);
    
    // Step forward in time with forward Euler
    hyperbolic_step(mesh, f, t, dt, ball_loc, water_tris);

    // Update node values with triangle-averaged values
    ball_loc = post_process(mesh, combined_forces, c, t, dt, water_nodes);

    // Update the viewer with new node positions
    viewer.add_nodes(mesh.node_begin(), mesh.node_end(), 
                     Color(water_nodes), NodePosition(), node_map);
    // viewer.add_triangles(mesh.triangle_begin(), mesh.triangle_end(), node_map);
    viewer.set_label(t);

    // find radius of cross sectional radius of ball submerged
    dh = ball_loc.z;
    if (dh > 2*radius)
      dh = 2 * radius;
    ball_loc.z = cross_radius(radius, dh);

    // These lines slow down the animation for small meshes.
    if (mesh.num_nodes() < 100)
      CS207::sleep(0.05);
  }
  return 0;
}
