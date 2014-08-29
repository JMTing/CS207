/**
 * @file mass_spring.cpp
 * Implementation of mass-spring system using Mesh
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

#include "Mesh.hpp"
#include "Point.hpp"


// Gravity in meters/sec^2
static constexpr double grav = 9.81;

// spring constant
static constexpr double spring_const = 35.0;

// spring constant
static constexpr double gas_const = 30;

// total mass of the mesh
static constexpr double total_mass = 225;

// interaction factor for wind force
static constexpr double wind_const = 0.00009;

// wind velocity constant
static const Point wind_velocity_const = Point(0.0, 1.0, 2.0);

// plane position constant
static constexpr double plane_const = -2.0;

// the time at which to add gas
static constexpr double t_addgas = 1.0;


//Define structures to pass in for node, edge and triangle values
struct node_values{
    double mass;
    Point velocity;
};

struct edge_values{
    double spring_constant;
    double initial_length;
};

struct triangle_values{
	// multiplier holds a -1 or 1, this will be set before the simulation
	// so that every time normal vectors to triangle surfaces are 
	// calculated, we multiply by this, so that we always get the 
	// normal vector that is pointing to the outside of the shape.
	double multiplier; 	
};

// setup typedefs
typedef Mesh<node_values, edge_values, triangle_values> MeshType;
typedef MeshType::Node Node;
typedef MeshType::Edge Edge;
typedef MeshType::Triangle Triangle;
typedef MeshType::size_type size_type;

// define type for forces. In this implementation, forces will be 
// points so that we can keep track of force in three dimensions
typedef Point force_type;

/** Calculate the point at the center of a mesh. 
 *  Calculated center = sum(node positions)/num_nodes() 
 *  @param[in,out] m Mesh
 *  @return a point that is the calculated center of the mesh
 */
template <typename MESH>
Point get_center(MESH& m)
{
  Point center = Point(0,0,0);
  for(auto it=m.node_begin(); it != m.node_end(); ++ it) {
    center += (*it).position();
  } 
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
 *	normal vector points away from the inside of the modeled
 *	shape and -1 if it points towards the inside
 */
void set_normal_direction(Triangle i, Point center){
  i.value().multiplier = 1.0;
  Point e01 = i.node(1).position() - i.node(0).position();
  Point e02 = i.node(2).position() - i.node(0).position();
  
  Point vec_normal = cross(e01, e02);
  Point vec_center = i.node(0).position() - center;
  if (dot(vec_normal, vec_center) < 0)
    i.value().multiplier = -1.0;
}

/** Find the normal vector to a triangle that points towards
 *	the outside of the original shape
 *  @param[in] i Triangle
 *	@returns A normal vector pointing away from the "inside" of
 *	the original shape
 */
Point get_normal_surface(Triangle& i){
  Point e01 = i.node(1).position() - i.node(0).position();
  Point e02 = i.node(2).position() - i.node(0).position();
  
  Point vec_normal = cross(e01, e02);
  
  vec_normal *= i.value().multiplier;
  vec_normal /= norm(vec_normal);
  return vec_normal;
}

/** Calculate the volume of a 3D mesh shape
 *	@param[in] m Mesh
 *  @return The approximate volume of the shape
 */
template <typename MESH>
double get_volume(MESH& m) {
  double volume = 0.0;
  for(auto it = m.triangle_begin(); it != m.triangle_end(); ++it){
    auto tri = (*it);
    Point normal = get_normal_surface(tri);
    double area = tri.area();
    volume += normal.z * area * (tri.node(0).position().z + tri.node(1).position().z + tri.node(2).position().z) / 3.0;
  }
  return volume;
}

/** Change a mesh's nodes according to a step of the symplectic Euler
 *  method with the given node force.
 *  @param[in,out] m Mesh
 *  @param[in] t The current time (useful for time-dependent forces)
 *  @param[in] dt The time step
 *  @param[in] force Function object defining the force per node
 *  @pre G::node_value_type supports position() function for getting a
 *  node's position as a Point type
 *  @pre G::node_value_type supports value() for getting a structure 
 *  containing mass and velocity
 *  @post all nodes' velocities are updated with the appropriate forces 
 *  on those nodes
 *  @return the next time step (usually @a t + @a dt)
 *
 *  @a force is called as @a force(n, @a t), where n is a node of the mesh
 *  and @a t is the current time parameter. @a force must return a Point
 *  representing the node's force at time @a t.
 */
template <typename M, typename F, typename C>
double symp_euler_step(M& m, double t, double dt, F force, C constraints) {
  double mass;
  
  // iterate through each node in the mesh and update position
  for(auto it=m.node_begin(); it != m.node_end(); ++ it)
  {
    (*it).set_position((*it).position() + (*it).value().velocity*dt);
  }

  // run nodes through constraint functors to correct values for violated
  // constraints
  constraints(m,t);
  

  // iterate through each node in the mesh and update velocity
  for(auto it=m.node_begin(); it != m.node_end(); ++ it)
  {
        mass = (*it).value().mass;
        (*it).value().velocity += force((*it),t)*dt/mass;
  }


  return t + dt;
}


/** Force function object for assessing the force due to gravity on 
 *  the nodes of the mesh. I use a struct rather than a class here
 *  since the functions default to public. */
struct GravityForce {
  /** Return the force due to gravity applying to @a n. 
  *  @pre n is templated with the node class from mesh.hpp
  *  @param[in] n The node whose forces we are calculating
  *  @param[in] t The time in the simulation. This is not actually used
  *  here, but we take it as a parameter for consistency since other
  *  forces may use time as a parameter.
  *  @post Returns a force due to gravity, where 
  *  F=n.value().mass*force_type(0,0,-grav)
  */
  template <typename NODE>
  force_type operator()(NODE n, double t) {
    (void) t;
    force_type force = force_type(0,0,0);
    double mass = n.value().mass;
    return (force + mass * force_type(0,0,-grav));
  }
};

/** Force function object for assessing the force due to springs on 
 *  the nodes of the mesh. I use a struct rather than a class here
 *  since the functions default to public. */
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
  template <typename NODE>
  force_type operator()(NODE n, double t) {
    (void) t;
    force_type force = force_type(0,0,0);
    Node node2;
    double distance; 
    double initial_spring_length;
    double sprint_const_this_edge;
    
    // iterate through each neighboring node and add spring forces
    for(auto it = n.edge_begin(); it != n.edge_end(); ++it)
    {
        node2 = (*it).node2();
        
        sprint_const_this_edge = (*it).value().spring_constant; 
        initial_spring_length = (*it).value().initial_length;
        
        distance = norm(n.position()-node2.position()); // Euclidean distance
        force += (-1)*sprint_const_this_edge*(n.position()-node2.position())*(distance-initial_spring_length)/distance;
    }
    return force;
  }
};

/** Force function object for assessing the damping force on 
 *  the nodes of the mesh. I use a struct rather than a class here
 *  since the functions default to public. */
struct DampingForce {
  // constructors
  DampingForce(){
  }
  
  DampingForce(int num)
  {
      num_nodes = num;
  }
  
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
  template <typename NODE>
  force_type operator()(NODE n, double t) {
    (void) t;
    Point velocity = n.value().velocity;
    Point return_force = -velocity/num_nodes;
    return Point(return_force.x,return_force.y,0);
  }
 
 private:
  int num_nodes; //holds the number of nodes in the mesh
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
   template <typename NODE>
   force_type operator() (NODE n, double t){
     Point node_velocity = n.value().velocity;
     Point node_normal(0.0, 0.0, 0.0);
     for (auto it = n.triangle_begin(); it != n.triangle_end(); ++it){
       auto tri = (*it);
       // approximate node normal vector by the sum of normal vectors of its neighboring faces
       node_normal +=  get_normal_surface(tri);
     }
     // calculte wind force
     force_type force = wind_const * dot(wind_velocity_ - node_velocity, node_normal) * node_normal;
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
  template <typename NODE>
  force_type operator()(NODE n, double t) {
    (void) t;
    force_type pressure_force = Point(0.0, 0.0, 0.0);

    for(auto it=n.triangle_begin(); it != n.triangle_end(); ++ it) {
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
  template <typename NODE>
  force_type operator()(NODE n, double t) {
      return f1(n,t)+f2(n,t);
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


/** returns the combination of the five templated forces. 
 *  @pre F1, F2, F3, F4, and F5 are all struct types that contain an overloaded 
 *  () operator with a templated first parameter, a double for the
 *  second parameter and an int for the third parameter. This 
 *  operator on F1, F2, F3, F4, and F5 must return a force_type value. 
 *  @param[in] f1 The first force structure
 *  @param[in] f2 The second force structure
 *  @param[in] f3 The third force structure
 *  @param[in] f4 The fourth force structure
 *  @param[in] f5 The fifth force structure
 *  @post Returns a force structure which is the result of combining the 
 *  forces passed in
 */
template<typename F1, typename F2, typename F3, typename F4, typename F5>
combined_force<combined_force<F1,F2>, combined_force<combined_force<F3,F4>, F5>> make_combined_force(F1& f1, F2& f2, F3& f3, F4& f4, F5& f5) {
    combined_force<F1,F2> two_combined = make_combined_force(f1,f2);
    auto three_combined = make_combined_force(f3,f4,f5);
    auto five_combined = make_combined_force(two_combined,three_combined);
    return five_combined;
}


/** Constraint function object for constraining the nodes of the mesh
 *  so they cannot go below -0.75 on the z plane. I use a struct rather 
 *  than a class here since the functions default to public. */
struct PlaneConstraint {
 private:
   double plane_;
 public:
   // public constructor: given plane position
   PlaneConstraint(const double& plane) : plane_(plane) {}
   // public constructor: use the default plane position
   PlaneConstraint() : plane_(plane_const) {}

  /** Move nodes that go below -0.75 on the z plane back to the closest
   *  point on the z plane at -0.75 and set these nodes' velocities on 
   *  the z axis to 0.
   *  @pre m is templated with the Mesh class from mesh.hpp
   *  @param[in] m The mesh whose nodes we are constraining
   *  @param[in] t The time in the simulation. This is not actually used
   *  here, but we take it as a parameter for consistency since other
   *  constraint structures may use time as a parameter.
   *  @post Nodes that were below -0.75 on the z plane are now at the 
   *  nearest location on the z plane at -0.75
   *  @post Nodes that were below -0.75 on the z plane now have 
   *  z-component velocities of 0.
   */  
  template <typename MESH>
  void operator()(MESH& m, double t) {
    (void) t;
    for(auto it=m.node_begin(); it != m.node_end(); ++ it)
    {
        if ((*it).position().z < plane_)
        {
            Point current_position = (*it).position();
            Point current_veloc = (*it).value().velocity;
            (*it).set_position(Point(current_position.x,current_position.y,plane_));
            (*it).value().velocity = Point(current_veloc.x,current_veloc.y,-current_veloc.z);
        }
    }
  }
};

/** Constraint function object for constraining the nodes of the graph
 *  so they cannot collied with each other as described below. I use a 
 *  struct rather than a class here since the functions default to 
 *  public. */
struct SelfCollisionConstraint {
  /** Constraint that nodes cannot collide. We will consider two nodes
   *  to be "collided" if the first node is within a distance of r
   *  from the other node, where r is the length of the shortest edge
   *  connected to the first node.
   *  @pre m is templated with the Mesh class from mesh.hpp
   *  @param[in] m The mesh whose nodes we are constraining
   *  @param[in] t The time in the simulation. This is not actually used
   *  here, but we take it as a parameter for consistency since other
   *  constraint structures may use time as a parameter.
   *  @post If a node n1 is found within the boundaries of the sphere
   *  defined for n2 (center = n2.position() and r = length of shortest 
   *  connected edge), the position of the node is set to the nearest 
   *  point on the surface of the sphere
   *  @post If a node n1 is found within the boundaries of the sphere
   *  defined for n2 (center = n2.position() and r = length of shortest 
   *  connected edge), the component of the velocity that is normal to 
   *  the sphere's surface is set to zero
   */      
  template <typename MESH>
  void operator()(MESH& m, double t) {
    (void) t;
    double distance;  
    double r;
    Point directionVector;
    Point pos;
    Point c;
    Point R;
    Point velocity;
    Node node1;
    Node node2;
    
    for(auto it=m.node_begin(); it != m.node_end(); ++ it)
    {
        
        node1 = (*it);
        c = node1.position();
        r = 0.0;
        
        // iterate through connected edges and set r equal to  
        // the length of the shortest connected edge
        for(auto edge_it = node1.edge_begin(); edge_it != node1.edge_end(); ++edge_it)
        {
            if (r == 0.0 || (*edge_it).length() < r)
                r = (*edge_it).length();
        }
        
        // now iterate through all of the other nodes and see if any of
        // them are within distance r
        for(auto it2=m.node_begin(); it2 != m.node_end(); ++it2)
        {
            node2 = (*it2);
            pos = node2.position();
            distance = norm(pos-c);
            if (distance < r && node1 != node2)
            {
                // SET POSITION TO NEAREST POINT ON SURFACE OF SPHERE
                
                // get the direction vector from the center to our node 
                directionVector = pos-c;
                //normalize that vector
                directionVector = directionVector/norm(directionVector);
                    
                // set pos to the direction vector time the sphere's radius
                // which should get us to the point on the sphere closest
                // to the node's current position
                pos = c+directionVector*(r);
                node2.set_position(pos);
                    
                // SET THE COMPONENT OF VELOCITY THAT IS NORMAL TO THE 
                // SPHERE'S SURFACE TO ZERO
                R = (pos-c)/norm(pos-c);
                velocity = node2.value().velocity;
                node2.value().velocity = velocity - inner_prod(velocity,R)*R;
            }
        }
    }
  }
};

/** This structure is templated so that it can take two constraint 
 *  structs such as PlaneConstraint or CornerConstraint above and put 
 *  them together into a combined structure. When the () operator is 
 *  called, this structure will run the () operator on each of the 
 *  combined structs
 */
template <typename C1, typename C2>
struct combined_constraints{
  // constructors
  combined_constraints() {}    
  combined_constraints(C1 a, C2 b) : c1(a), c2(b)  {}
  
  /** Runs the () operators on the two combined constraint structures.
   *  @pre m is templated with the Mesh class from mesh.hpp
   *  @param[in] m The mesh whose nodes we are constraining
   *  @param[in] t The time in the simulation
   *  @post The () operator has been run on the two constraint structs 
   *  that are combined together in this struct. 
   */
  template <typename MESH>
  void operator()(MESH& m, double t) {
      c1(m,t);
      c2(m,t);
  }

 private:
  C1 c1;
  C2 c2;
};

/** returns the combination of two templated constraint structures. 
 *  @pre C1 and C2 are both struct types that contain an overloaded 
 *  () operator with a templated first parameter, and a double for the
 *  second parameter. This operator must not return a value.
 *  @param[in] c1 The first constraint structure
 *  @param[in] c2 The second constraint structure
 *  @post Returns a constraint structure which is the result of 
 *  combining the two constraint structures passed in
 */
template<typename C1, typename C2>
combined_constraints<C1, C2> make_combined_constraints(C1 a, C2 b) {
    auto return_value = combined_constraints<C1,C2>(a,b);
    return return_value;
}

/** returns the combination of three templated constraint structures. 
 *  @pre C1, C2 and C3 are all struct types that contain an overloaded 
 *  () operator with a templated first parameter, and a double for the
 *  second parameter. This operator must not return a value.
 *  @param[in] c1 The first constraint structure
 *  @param[in] c2 The second constraint structure
 *  @param[in] c3 The second constraint structure
 *  @post Returns a constraint structure which is the result of 
 *  combining the three constraint structures passed in
 */
template<typename C1, typename C2, typename C3>
combined_constraints<combined_constraints<C1,C2>,C3> make_combined_constraints(C1 a, C2 b, C3 c) {
    combined_constraints<C1,C2> two_combined = make_combined_constraints(a,b);
    auto three_combined = make_combined_constraints(two_combined,c);
    return three_combined;
}

/** see specs for make_combined_constraints(C1,C2,C3) above. The specs
 *  for this function are the same, but with a 4th constraint added*/
template<typename C1, typename C2, typename C3, typename C4>
combined_constraints<combined_constraints<C1,C2>,combined_constraints<C3,C4>> make_combined_constraints(C1 a, C2 b, C3 c, C4 d) {
    combined_constraints<C1,C2> two_combined1 = make_combined_constraints(a,b);
    combined_constraints<C3,C4> two_combined2 = make_combined_constraints(c,d);
    auto four_combined = make_combined_constraints(two_combined1,two_combined2);
    return four_combined;
}






int main(int argc, char* argv[]) {
  // Check arguments
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  MeshType mesh;
  std::vector<typename MeshType::node_type> nodes;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  
  // Interpret each line of the nodes_file as a 3D Point and add to the Mesh
  Point p;
  while (CS207::getline_parsed(nodes_file, p))
    nodes.push_back(mesh.add_node(p));

  // Create a trianges_file from the second input argument
  std::ifstream triangles_file(argv[2]);

  // Interpret each line of the tets_file as three ints which refer to nodes
  std::array<int,3> t;


  // add in the triangles
  while (CS207::getline_parsed(triangles_file, t))
    for (unsigned i = 1; i < t.size(); ++i)
      for (unsigned j = 0; j < i; ++j)
        for (unsigned k = 0; k < j; ++k)
        {
          mesh.add_triangle(nodes[t[i]], nodes[t[j]], nodes[t[k]]);
      }   

  // Set masses of nodes equal to 1/N where N is the number of nodes
  // and the initial velocities to 0. Also, get the indexes of
  // the nodes at positions (1,0,0) and (0,0,0)
  for(auto it=mesh.node_begin(); it != mesh.node_end(); ++ it)
  {
      (*it).value().mass = total_mass/mesh.num_nodes();
      (*it).value().velocity = Point(0.0,0.0,0.0);
  }

  // Set spring constants for each node equal to spring_const variable
  // and set initial length values equal to lengths of edges prior to
  // running the symplectic Euler steps
  for(auto it=mesh.edge_begin(); it != mesh.edge_end(); ++ it)
  {
      (*it).value().spring_constant = spring_const;
      (*it).value().initial_length = (*it).length();
  }


  // Set the triangle direction values so that we can determine which 
  // way to point normal vectors. This part assumes a convex shape
  Point center = get_center(mesh);
  for(auto it=mesh.triangle_begin(); it != mesh.triangle_end(); ++ it)
  {  
  	//std::cout << (*it).index() << std::endl;
  	set_normal_direction((*it),center);
  }

  // Print out the stats
  std::cout << mesh.num_nodes() << " " << mesh.num_edges() << std::endl;

  std::cout << "Center: " << get_center(mesh) << std::endl;
  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  auto node_map = viewer.empty_node_map(mesh);
  viewer.launch();

  viewer.add_nodes(mesh.node_begin(), mesh.node_end(), node_map);
  viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);

  viewer.center_view();

  // Begin the mass-spring simulation
  double dt = 0.0001;
  double t_start = 0;
  double t_end = 20.0;
  

  // Initialize constraints
  PlaneConstraint c1;
  //SelfCollisionConstraint c2;
  //auto combined_constraints = make_combined_constraints(c1,c2);
  


  for (double t = t_start; t < t_end; t += dt) {
    MassSpringForce ms_force;
    PressureForce p_force = PressureForce(0.0);
    DampingForce d_force = DampingForce(mesh.num_nodes());
    GravityForce g_force;
    WindForce w_force;

    (void) d_force; // prevents compiler from throwing error for unused variable
      
    
    if (t >= t_addgas - dt) {
      p_force.set_pressure(gas_const/get_volume(mesh));
      if (t < t_addgas)
        std::cout << "Adding gas to the ball now..." << std::endl;
    }

    auto combined_forces = make_combined_force(ms_force, p_force, w_force, g_force);

    symp_euler_step(mesh, t, dt, combined_forces, c1);

    // Update viewer with nodes' new positions
    viewer.add_nodes(mesh.node_begin(), mesh.node_end(), node_map);

    // update the viewer's label with the ball center's position on the z axis 
    viewer.set_label(get_center(mesh).z);

  }
  return 0;
}
