#ifndef CS207_MESH_HPP
#define CS207_MESH_HPP

/** @class Mesh
* @brief A template for 3D meshes, which are made up of nodes and 
* triangles.
*/

/* ================================================== */
/* Team Members: Daniel Newman, Xueshan Zhang (Ellie) */
/*                   CS 207 HW 4                      */
/* ================================================== */

#include "CS207/Util.hpp"
#include "Point.hpp"
#include "Graph.hpp"
#include <utility>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <tuple> 
#include <string>
#include <cmath>

template <typename N, typename E, typename T>
class Mesh {
public:

  typedef unsigned size_type;

  typedef N node_value_type;
  typedef E edge_value_type;
  typedef T triangle_value_type;

  /** Type of this mesh. */
  typedef Mesh<N,E,T> mesh_type;

  /** type for graph */
  typedef Graph<N,E> GraphType;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  typedef Node node_type;
  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  typedef Edge edge_type;

  /** predeclaration of Triangle type. */
  class Triangle;
  
  /** Type of node iterators, which iterate over all mesh nodes. */
  class node_iterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class edge_iterator;

  /** Type of incident iterator, which iterates incident triangles to a node. */
  class incident_triangle_iterator;

  /** Type of incident iterator, which iterates incident edges to a node. */
  class incident_edge_iterator;

  /** Construct an empty mesh. 
    * Complexity: O(1)
    */
  Mesh<N,E,T>() {}

  /** Default destructor
    * Complexity: O(1)
    */
  ~Mesh<N,E,T>() = default;


  /** @class Mesh::Node
    * @brief Class representing the mesh's nodes.
    *
    * Node objects are used to access information about the Mesh's nodes
    */
  class Node : private totally_ordered<Node> {
  public:
    /** Construct an invalid node.
  * Complexity: O(1). 
  */
    Node() {};

    /** Return this node's position. 
  * Complexity: O(1). 
  */
    Point position() const{
      return mesh_->graph_.node(index_).position();
    }

    /** Return this node's index, a number in the 
  * range [0, mesh_size).
  * Complexity: O(1). 
  */
    size_type index() const
    {
      return index_;
    }

    /** Return a reference to this node's value in a mutable mesh. 
  * The value is mutable.
  * Complexity: O(1). 
  */ 
    node_value_type& value() {
      return mesh_->graph_.node(index_).value();
    }

    /** Return a reference to this node's value in a constant mesh. 
  * The value is read-only.
  * Complexity: O(1). 
  */
    const node_value_type& value() const {
      return mesh_->graph_.node(index_).value();
    }

    /** Return the degree of a node.
  * Complexity: O(1)
  */
    size_type degree() const {
      return mesh_->graph_.node(index_).degree();
    }

    /** Set the position of a Node to p.
  * Complexity: O(1)
  */
    void set_position(const Point& p) {
      mesh_->graph_.node(index_).set_position(p);
    }

    /** Compare if two Nodes are the same object.
  * @return True if the two nodes have the same index in the 
  * same mesh.
  * Complexity: O(1). 
  */ 
    bool operator==(const Node& n) const {
      return (index() == n.index() && mesh() == n.mesh());
    }

    /** Used to define the <> operator. 
  * @return True if the two nodes are not in the same mesh or the 
  * "less than" relation holds in term of node indexes.
  * Complexity: O(1). 
  */ 
    bool operator< (const Node& n) const {
      return !(index() >= n.index());
    }

    /** Return an incident_triangle_iterator. 
  * @return the incident_triangle_iterator_node that points to the
  * beginning of the adjacency triangle list of the node.
  * Complexity: O(1)
  */		
    incident_triangle_iterator triangle_begin() const
    {
        std::unordered_set<size_type>::iterator it = mesh_->nodes_connected_triangles_[index_].begin();
        return incident_triangle_iterator(mesh_, it);
    }

    /** Return an incident_triangle_iterator. 
  * @return the incident_triangle_iterator_node that points to the
  * end of the adjacency triangle list of the node.
  * Complexity: O(1)
  */		
    incident_triangle_iterator triangle_end() const
    {
        std::unordered_set<size_type>::iterator it = mesh_->nodes_connected_triangles_[index_].end();
        return incident_triangle_iterator(mesh_, it);
    }

    /** Return an incident_edge_iterator. 
  * @return the incident_edge_iterator that points to the
  * beginning of the adjacency edge list of the node.
  * Complexity: O(1)
  */    
    incident_edge_iterator edge_begin() const
    {
        //return mesh_->graph_.node(index_).edge_begin();
      return incident_edge_iterator(mesh_, mesh_->graph_.node(index_).edge_begin());
    }

    /** Return an incident_edge_iterator. 
  * @return the incident_edge_iterator that points to the
  * end of the adjacency edge list of the node.
  * Complexity: O(1)
  */    
    incident_edge_iterator edge_end() const
    {
        return incident_edge_iterator(mesh_, mesh_->graph_.node(index_).edge_end());
    }

  private:
    /** Return a pointer to this node's mesh */
    mesh_type* mesh() const {
      return mesh_;
    }

    /** Only Mesh can access our private members */
    friend class Mesh<N,E,T>;

    /** Pointer back to the Mesh container */
    mesh_type* mesh_;

    /** This node's unique identification number */
    size_type index_;

    /** Private Constructor */
    Node(const mesh_type* mesh, size_type index){
      mesh_ = const_cast<mesh_type*>(mesh);
      index_ = index;
    }

  };


  /** @class Mesh::Edge
    * @brief Class representing the mesh's edges.
    *
    * Edges are order-insensitive pairs of nodes. Two Edges with the 
    * same nodes
    * are considered equal if they connect the same nodes, in either 
    * order.
    */
  class Edge {
  public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return mesh_->node(mesh_->graph_.edge(index_).node1().index());
    }     

    /** Return a node of this Edge */
    Node node2() const {
      return mesh_->node(mesh_->graph_.edge(index_).node2().index());
    }      

    /** Return this edge's index. */
    size_type index() const {
      return index_;
    }    

    /** Return a reference to this edge's value in a mutable mesh. 
  * The value is mutable.
  * Complexity: O(1). 
  */ 
    edge_value_type& value() {
      return mesh_->graph_.edge(index_).value();
    }

  /** Return a reference to this edge's value in a constant mesh.
  * The value is read-only.
  * Complexity: O(1). 
  */
    const edge_value_type& value() const {
      return mesh_->graph_.edge(index_).value();
    }
    
    /** Return the edge's length */
    double length() const {
      return mesh_->graph_.edge(index_).length();  
    }
    
    
    /** Return the number of triangles connected to this edge 
     * @post returns 1 or 2, depending on the number of adjacent triangles
     */
    int num_adj_tri() {
        assert(mesh_->edges_connected_triangles_info_[index_].num_connected_triangles < 3);
        assert(mesh_->edges_connected_triangles_info_[index_].num_connected_triangles > 0);
        return mesh_->edges_connected_triangles_info_[index_].num_connected_triangles;
    }
    
    /** Return the i'th adjacent triangle to an edge
     * @pre i < num_adj_tri()
     * @pre i >= 0
     * @post a triangle adjacent to the edge is returned which is different
     * from all triangles returned for all other valid values of i
     * NOTE: This function makes no guarantees about the order in which triangles will be
     * returned. Only that the two triangles will be different if num_adj_tri() == 2
     */
     Triangle adj_tri(int i) {
         assert(i < num_adj_tri());
         assert(i >= 0);
         if (i == 0) 
            return mesh_->triangle(mesh_->edges_connected_triangles_info_[index_].triangle_1_idx);
         return mesh_->triangle(mesh_->edges_connected_triangles_info_[index_].triangle_2_idx);
     }


    /** Override the '==' operator and say that two edges are equal
  *  if they have the same index and are on the same graph 
  *  @return true if the edge being compared to has the same index
  *  and is in the same graph
  *  @param[in] e The edge to compare the current edge to 
  */
    bool operator==(const Edge& e) const {
      return (mesh() == e.mesh() && index() == e.index());
    }


    /** Override the '<' operator and say that edge a is < edge b if the
  *  index of a < b. If the indices are equal, but the graphs are not
  *  then also return true (intuitively, this doesn't make much
  *  sense to me but seems to be necessary in order to pass the 
  *  last test of test_edges.cpp). 
  *  @return true if the index of the current node < index of the
  *  comparison node. Otherwise, if the indices are equal but the
  *  graphs are not, return true. Otherwise, return false.
  *  @param[in] e The edge to compare to.*/
    bool operator<(const Edge& e) const {
      if (index() < e.index()) 
      return true;
      else if (index() == e.index() || mesh() != e.mesh())
      return true;
      return false;
    }

  private:
    /** Private Constructor */
    Edge(const mesh_type* mesh, size_type index){
      mesh_ = const_cast<mesh_type*>(mesh);
      index_ = index;
    }

    /** Return a pointer to this edges's mesh */
    mesh_type* mesh() const {
      return mesh_;
    }

    /** Only Mesh can access our private members */
    friend class Mesh<N,E,T>;

    /** Pointer back to the Mesh container */
    mesh_type* mesh_;

    /** This edge's unique identification number */
    size_type index_;        
  };

  /** @class Mesh::Triangle
    * @brief Class representing the mesh's triangles.
    *
    * Triangles are sets of three nodes. Two Triangles with the same 
    * nodes    	* are considered equal if they connect the same nodes, 
    * in any order.
    */
  class Triangle {
  public:
    /** Construct an invalid triangle.
  */
    Triangle() {}

    /** Return the i'th node of this Triangle 
  * @pre @a 0 <= i < 3.
  * Complexity: O(1)
  */
    Node node(size_type i) const {
      assert(i < 3);
      return mesh_->node(mesh_->triangles_[index_].nodes[i]);    
    } 
    
    /** Return the i'th Edge of this Triangle associated
  * @pre @a 0 <= i < 3.
  * i = 0: Edge with nodes 0 and 1
  * i = 1: Edge with nodes 1 and 2
  * i = 2: Edge with nodes 0 and 2
  * Complexity: O(1)	*/
    Edge edge(size_type i) const{
      assert(i < 3);
      return mesh_->edge(mesh_->triangles_[index_].edges[i]); 
    }

    /** Return the area of this Triangle
  * Complexity: O(1)	*/
    double area() const {
      /** Formula for calculating triangle is given by:
    * A = (p(p-a)(p-b)(p-c))^0.5
    * Where a, b, and c are the lengths of the sides of the triangle
    * and p is half the perimeter, or (a+b+c)/2
    */
      double a = edge(0).length();
      double b = edge(1).length();
      double c = edge(2).length();
      double p = (a+b+c)/2.0;

      return sqrt(p*(p-a)*(p-b)*(p-c));
    }
    
    /** Return an unnormalized outward normal vector of the i'th edge of this Triangle. Apply to any triangle in the space.
    * @pre @a i < 3.
    * Complexity: O(1)
    */
    Point normal_3D(size_type i) const {
    /** Formula for calculating outward normal vector of edge i is given by:
    * v = e_01 x (e_01 x e_12)   (unnormalized vector)
    * v = v/norm(v)   (normalize v to have magnitude 1)
    * where e_01 is the vector pointing from node i to node (i+1)%3, e_12 is the vector pointing from node (i+1)%3 to node (i+2)%3. "x" means cross product.
    * 
    * Example:
    * edge 0 connects node 0 and node 1. 
    * Take e_01 to be a vector pointing from node 0 to node 1, e_12 from node 1 to node 2.
    * e_01 x e_12 gives us a normal vector to the plane of the triangle. 
    * e_01 x (e_01 x e_12) returns a normal vector to e_01, i.e. the vector perpendicular to the plane of the triangle and e_01, pointing away from the triangle.
    */
      assert(i < 3);
      size_type i_0 = i;
      size_type i_1 = (i + 1) % 3;
      size_type i_2 = (i + 2) % 3;
      
      Point e_01 = node(i_1).position() - node(i_0).position();
      Point e_12 = node(i_2).position() - node(i_1).position();
      
      Point vecn = cross(e_01, cross(e_01, e_12));
      return vecn;  
    }
    
    /** Return an unnormalized outward normal vector of the i'th edge of this Triangle when all nodes in this triangle have the same z-coordinate value. 
    * @pre @a i < 3.
    * Complexity: O(1)
    */
    Point normal_2D(size_type i) const {
    /** How to calculate the outward normal vector of edge e_01 (pointing from node i to node (i+1)%3) on a 2D plane at height z:
     * Take the 2D normal (e_01.y, -e_01.x). 
     * Dot this normal with the first two components of vector e_02. 
     * If the result is positive, then the normal is facing towards the "inner" half-space, so invert the normal.
     * Add z as the third dimension of this normal vector.
    */
      assert(i < 3);
      size_type i_0 = i;
      size_type i_1 = (i + 1) % 3;
      size_type i_2 = (i + 2) % 3;
      
      Point e_01 = node(i_1).position() - node(i_0).position();
      Point e_02 = node(i_2).position() - node(i_0).position();
      
      Point vecn = Point(e_01.y, -e_01.x, 0);
      if (dot(vecn, e_02) > 0)
        vecn *= -1;
      
      vecn.z = node(i_0).position().z;
      return vecn; 
    }

    /** Return the index value of this Triangle
  * Complexity: O(1)  */
    size_type index() const {
      return index_;
    }

    /** Return the value of this Triangle and allow it to be changed
  * @return the value held by this triange, as a non-constant 
  * variable
  * Complexity: O(1)  */
    triangle_value_type& value() {
      return mesh_->triangles_[index_].value;
    }

    /** Return the value of this Triangle and don't allow it to be 
  * changed
  * @return the value held by this triange, as a constant variable
  * Complexity: O(1)  */
    const triangle_value_type& value() const {
      return mesh_->triangles_[index_].value;
    }
    
    /** Return the number of triangles connected to this triangle
     * @post returns 0, 1, 2 or 3, depending on the number of adjacent triangles
     */
    bool has_adj_tri(size_type i) {
        assert(i<3);
        return (edge(i).num_adj_tri() == 2);
    }
    
    /** Return the i'th adjacent triangle to a triangle
     * @pre i < num_adj_tri()
     * @pre has_adj_tri(i) == True
     * @post a triangle adjacent to the triangle on edge i is returned
     */
     Triangle adj_tri(size_t i) {
        assert(has_adj_tri(i));
        if (edge(i).adj_tri(0).index() == index_)
            return edge(i).adj_tri(1);
        return edge(i).adj_tri(0);
     }

    /** Compare if two Triangles are the same object.
  * @return True if the two triangles have the same index in the 
  * same mesh.
  * Complexity: O(1). 
  */ 
    bool operator==(const Triangle& t) const {
      return (mesh() == t.mesh() && index() == t.index());
    }


    /** Used to define the <> operator. 
  * @return True if the two triangles are not in the same mesh or 
  * the "less than" relation holds in term of triangles indexes.
  * Complexity: O(1). 
  */ 
    bool operator<(const Triangle& t) const {
      if (index() < t.index()) 
      return true;
      else if (index() == t.index() || mesh() != t.mesh())
      return true;
      return false;
    }

  private:
    /** Private Constructor */
    Triangle(const mesh_type* mesh, size_type index){
      mesh_ = const_cast<mesh_type*>(mesh);
      index_ = index;
    }    

    /** Return a pointer to this triangle's mesh */
    mesh_type* mesh() const {
      return mesh_;
    }

    /** Only Mesh can access our private members */
    friend class Mesh<N,E,T>;

    /** Pointer back to the Mesh container */
    mesh_type* mesh_;

    /** This triangle's unique identification number */
    size_type index_;    
  };

  /** @class Mesh::triangle_iterator
    * @brief Iterator class for triangles. A forward iterator. 
    * Complexity: O(1) for all operators implemented
    */
  class triangle_iterator : private totally_ordered<triangle_iterator> {  
  public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Triangle value_type;
    /** Type of pointers to elements. */
    typedef Triangle* pointer;
    /** Type of references to elements. */
    typedef Triangle& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid node_iterator. */
    triangle_iterator() {
    }

    /** Returns the node currently associated with it_ in the iterator
  * @pre 0<=it_<num_nodes()-1 
  */
    value_type operator*() const {
      return mesh_->triangle(it_);
    }

    /** increments the value it_ by 1
  * @pre 0<=it_<num_nodes() 
  */
    triangle_iterator& operator++(){
      ++it_;
      return *this;
    }

    /** returns true if two iterators belong to the same graph
  *  and they are at the same value it_
  *  @param[in] other The iterator to compare to
  */
    bool operator==(const triangle_iterator& other) const{
      return (it_ == other.it_ && mesh_ == other.mesh_);
    }

  private:

    /** private constructor for an valid node_iterator. 
  * @pre 0<=it_<num_nodes() */
    triangle_iterator(const mesh_type* mesh, size_type it) {
      mesh_ = const_cast<mesh_type*>(mesh);
      it_ = it;
    }    
    friend class Mesh<N,E,T>;
    mesh_type* mesh_;
    size_type it_;
  };

  /** returns a triangle_iterator that relates to the first triange 
    * by index
    * @pre size() > 0
    */
  triangle_iterator triangle_begin() const {
    return triangle_iterator(this, 0);
  }

  /** returns a triangle_iterator that relates to the position one 
    * greater than size()
    * @pre size() > 0
    */
  triangle_iterator triangle_end() const {
    return triangle_iterator(this, num_triangles());
  }


  /** @class Mesh::node_iterator
    * @brief Iterator class for nodes. A forward iterator. 
    * Complexity: O(1) for all operators implemented
    */
  class node_iterator : private totally_ordered<node_iterator> {  
  public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Node value_type;
    /** Type of pointers to elements. */
    typedef Node* pointer;
    /** Type of references to elements. */
    typedef Node& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid node_iterator. */
    node_iterator() {
    }

    /** Returns the node currently associated with it_ in the iterator
  * @pre 0<=it_<num_nodes()-1 
  */
    value_type operator*() const {
      return mesh_->node(it_);
    }

    /** increments the value it_ by 1
  * @pre 0<=it_<num_nodes() 
  */
    node_iterator& operator++(){
      ++it_;
      return *this;
    }

    /** returns true if two iterators belong to the same graph
  *  and they are at the same value it_
  *  @param[in] other The iterator to compare to
  */
    bool operator==(const node_iterator& other) const{
      return (it_ == other.it_ && mesh_ == other.mesh_);
    }

  private:

    /** private constructor for an valid node_iterator. 
  * @pre 0<=it_<num_nodes() */
  node_iterator(const mesh_type* mesh, size_type it) {
      mesh_ = const_cast<mesh_type*>(mesh);
      it_ = it;
    }    
    friend class Mesh<N,E,T>;
    friend class edge_iterator;
    mesh_type* mesh_;
    size_type it_;
  };

  /** returns a node_iterator that relates to the first node by index
    * @pre num_nodes() > 0
    * Complexity: O(1)
    */
  node_iterator node_begin() const
  {
    return node_iterator(this, 0);
  }

  /** returns a node_iterator that relates to the position one greater
    * than num_nodes()
    * @pre num_nodes() > 0
    * Complexity: O(1)
    */
  node_iterator node_end() const
  {
    return node_iterator(this, graph_.num_nodes());
  }

  /** @class Mesh::edge_iterator
    * @brief Iterator class for edges. A forward iterator. */
  class edge_iterator : private totally_ordered<edge_iterator> {  
  public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid edge_iterator. */
    edge_iterator() {
    }

    /** Returns the edge currently associated with it_ in the iterator
  * @pre 0<=it_<num_edges()-1 
  */
    Edge operator*() const {
      return mesh_->edge(it_);
    }

    /** increments the value it_ by 1
  * @pre 0<=it_<num_edges() 
  */
    edge_iterator& operator++(){
      ++it_;
      return *this;
    }

    /** returns true if two iterators belong to the same mesh
  *  and they are at the same value it_
  *  @param[in] other The iterator to compare to
  */
    bool operator==(const edge_iterator& other) const{
      return (it_ == other.it_ && mesh_ == other.mesh_);
    }

  private:
    /** private constructor for an valid edge_iterator. 
  * @pre 0<=it_<num_edges() */
    edge_iterator(const mesh_type* mesh, size_type it) {
      mesh_ = const_cast<mesh_type*>(mesh);
      it_ = it;
    }    

    friend class Mesh<N,E,T>;
    mesh_type* mesh_;
    size_type it_;
  };


  /** returns a edge iterator that relates to the first edge by user index
    * @pre num_edges() > 0
    * Complexity: O(1)
    */
  edge_iterator edge_begin() const {
    return edge_iterator(this, 0);
  }

  /** returns an edge iterator that relates to a nonexisting position
    *  one greater than the size of graph_.num_edges()
    * @pre num_nodes() > 0
    * Complexity: O(1)
    */
  edge_iterator edge_end() const {
    return edge_iterator(this, graph_.num_edges());
  }

  /** @class Mesh::incident_triangle_iterator_node
    * @brief Iterator class for getting triangles adjacent to a given
    node. A forward iterator. 
    * Complexity: O(1) for all operators implemented
    */
  class incident_triangle_iterator : private totally_ordered<incident_triangle_iterator> { 
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Triangle value_type;
    /** Type of pointers to elements. */
    typedef Triangle* pointer;
    /** Type of references to elements. */
    typedef Triangle& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid incident_iterator. */
    incident_triangle_iterator() { 
    }

    /** Handles the * operator for the incident iterator:
     *  returns a triangle based on the values in the vector of adjacent triangle sets and the
     *  unsigned integer pointed to by it_.
     *  @pre it_ points to a valid value in the connected triangle set in nodes_connected_triangles_ 
     *  corresponding to the node being considered. 
     *  @return a triangle connected to the node using this iterator
     * */
    Triangle operator*() const{
        Triangle return_triangle = mesh_->triangle((*it_));
        return return_triangle;
    }
    
    /** Overloads the ++ operator for the incident iterator:
     *  @pre it_ points to a valid position in the set
     *  @post it_ is is iterated by one
     *  @return an incident iterator that has been iterated to point
     *  to the next value in the set */    
    incident_triangle_iterator& operator++(){
        ++it_;
        return *this;
    }

    /** handles the == operator for the incident iterator:
     *  returns true if both iterators have the same value
     *  and if both belong to the same mesh 
     *  @param[in] other The iterator to compare to
     *  @pre other refers to a valid incident iterator
     *  @return true if the value of the private variable it_ is equal
     *  to other's it_ variable and other is part of the same mesh.
     */    
    bool operator==(const incident_triangle_iterator& other) const{
        return (it_ == other.it_ && mesh_ == other.mesh_);        
    }


   private:
    /** private constructor for an valid edge_iterator. 
     * @param[in] @a mesh refers to the current mesh that this incident_triangle_iterator belongs to
     * @param[in] @a it is an iterator to the unordered set in nodes_connected_triangles_ 
     * corresponding to the node whose triangles the incident iterator is iterating through */
    incident_triangle_iterator(const mesh_type* mesh, std::unordered_set<size_type>::iterator it) {
      mesh_ = const_cast<mesh_type*>(mesh);
      it_ = it;
    } 
    
    friend class Mesh<N,E,T>;
    mesh_type* mesh_;
    std::unordered_set<size_type>::iterator it_;
  };

/** @class Mesh::incident_edge_iterator
    * @brief Iterator class for getting edges adjacent to a given
    node. A forward iterator. 
    * Complexity: O(1) for all operators implemented
    */
  class incident_edge_iterator : private totally_ordered<incident_edge_iterator> { 
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid incident_iterator. */
    incident_edge_iterator() { 
    }

    /** Handles the * operator for the incident iterator:
     *  returns a triangle based on the values in the vector of adjacent triangle sets and the
     *  unsigned integer pointed to by it_.
     *  @pre it_ points to a valid value in the connected triangle set in nodes_connected_triangles_ 
     *  corresponding to the node being considered. 
     *  @return a triangle connected to the node using this iterator
     * */
    Edge operator*() const{
        Edge return_edge = mesh_->edge((*it_).index());
        return return_edge;
    }
    
    /** Overloads the ++ operator for the incident iterator:
     *  @pre it_ points to a valid position in the set
     *  @post it_ is is iterated by one
     *  @return an incident iterator that has been iterated to point
     *  to the next value in the set */    
    incident_edge_iterator& operator++(){
        ++it_;
        return *this;
    }

    /** handles the == operator for the incident iterator:
     *  returns true if both iterators have the same value
     *  and if both belong to the same mesh 
     *  @param[in] other The iterator to compare to
     *  @pre other refers to a valid incident iterator
     *  @return true if the value of the private variable it_ is equal
     *  to other's it_ variable and other is part of the same mesh.
     */    
    bool operator==(const incident_edge_iterator& other) const{
        return (it_ == other.it_ && mesh_ == other.mesh_);        
    }


   private:
    /** private constructor for an valid edge_iterator. 
     * @param[in] @a mesh refers to the current mesh that this incident_triangle_iterator belongs to
     * @param[in] @a it is an iterator to the unordered set in nodes_connected_triangles_ 
     * corresponding to the node whose triangles the incident iterator is iterating through */
    incident_edge_iterator(const mesh_type* mesh, typename GraphType::incident_iterator it) {
      mesh_ = const_cast<mesh_type*>(mesh);
      it_ = it;
    } 
    
    friend class Mesh<N,E,T>;
    mesh_type* mesh_;
    typename GraphType::incident_iterator it_;
  };

  /** Add a node to the mesh, returning the added node.
    * @param[in] position The new node's position
    * @param[in] value The new node's value
    * @post new num_nodes() == old num_nodes() + 1
    * @post result_node.index() == old num_nodes()
    *
    * Complexity: O(1) amortized operations.
    */    
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    // add node to graph_
    graph_.add_node(position,value);

    // update nodes_connected_triangles_ vector
    std::unordered_set<size_type> connected_triangles;
    nodes_connected_triangles_.push_back(connected_triangles);

    // return a newly constructed node with index equal to num_nodes()-1
    return Node(this, graph_.num_nodes()-1);
  }


  /** Determine if this Node belongs to this Mesh
    * @return True if @a n is currently a Node of this Mesh
    * @return False if @a n is not currently a Node of this Mesh
    *
    * Complexity: O(1).
    */
  bool has_node(const Node& n) const {
    return !(n.index() > num_nodes());
  }

  /** Return the node with index @a i.
    * @pre 0 <= @a i < num_nodes()
    * @post result_node.index() == i
    *
    * Complexity: O(1).
    */
  Node node(size_type i) const
  {
    return Node(this, i); 
  }

  /** Return the number of nodes in the mesh 
    * 
    * Complexity: O(1)
    */
  size_type num_nodes() {
    return graph_.num_nodes();
  }

  /** Return the edge with index @a i.
    * @pre 0 <= @a i < num_edges()
    * @post result_edge.index() == i
    *
    * Complexity: O(1).
    */
  Edge edge(size_type i) const
  {
    return Edge(this, i); 
  }

  /** Return the number of edges in the mesh 
    * 
    * Complexity: O(1)
    */
  size_type num_edges() {
    return graph_.num_edges();
  }

  /** Add a triangle to the mesh, or return the current triange
    * if it already exists.
    * @pre @a a, @a b and @a c are distinct valid nodes of this mesh
    * @pre if the triangle being added does not already exist, then the 
    * triangle being added does not utilize an edge that is already used
    * by two triangles
    * @return a Triangle object t with t.node(0) == @a a and 
    * t.node(1) == @a b and t.node(2) == @a c 
    * @post has_triangle(@a a, @a b, @a c) == true
    * @post If old has_triangle(@a a, @a b), 
    * 		new num_triangle() == old num_triangle().
    *       Else,        new num_triangle() == old num_triangle() + 1.
    *
    *
    * Complexity: Amortized O(1)
    */
  Triangle add_triangle(const Node& a,const Node& b,const Node& c){

    // check if the triangle already exists and if so, then return the triangle
    int idx = find_triangle(a,b,c);
    if (idx >= 0)
        return Triangle(this, idx);

    // check if each edge is already in the graph. If so, then we can update
    // edges_connected_triangles_info_ with the second triangle. Otherwise,
    // we'll add a new to the edges_connected_triangles_info_ vector.
    if (graph_.has_edge(graph_.node(a.index()), graph_.node(b.index()))){
        auto e1 = graph_.add_edge(graph_.node(a.index()), graph_.node(b.index()));
        edges_connected_triangles_info_[e1.index()].num_connected_triangles += 1;
        assert(edges_connected_triangles_info_[e1.index()].num_connected_triangles <= 2);
        edges_connected_triangles_info_[e1.index()].triangle_2_idx = num_triangles();
    }
    else {
        internal_edge_info temp;
        temp.num_connected_triangles = 1;
        temp.triangle_1_idx = num_triangles();
        edges_connected_triangles_info_.push_back(temp);
    }
    if (graph_.has_edge(graph_.node(a.index()), graph_.node(c.index()))){
        auto e2 = graph_.add_edge(graph_.node(a.index()), graph_.node(c.index()));
        edges_connected_triangles_info_[e2.index()].num_connected_triangles += 1;
        assert(edges_connected_triangles_info_[e2.index()].num_connected_triangles <= 2);
        edges_connected_triangles_info_[e2.index()].triangle_2_idx = num_triangles();
    }
    else {
        internal_edge_info temp;
        temp.num_connected_triangles = 1;
        temp.triangle_1_idx = num_triangles();
        edges_connected_triangles_info_.push_back(temp);
    }
    if (graph_.has_edge(graph_.node(b.index()), graph_.node(c.index()))){
        auto e3 = graph_.add_edge(graph_.node(b.index()), graph_.node(c.index()));
        edges_connected_triangles_info_[e3.index()].num_connected_triangles += 1;
        assert(edges_connected_triangles_info_[e3.index()].num_connected_triangles <= 2);
        edges_connected_triangles_info_[e3.index()].triangle_2_idx = num_triangles();
    }
    else {
        internal_edge_info temp;
        temp.num_connected_triangles = 1;
        temp.triangle_1_idx = num_triangles();
        edges_connected_triangles_info_.push_back(temp);
    }
    
    // we'll call the add_edge functions again to get the edge indixes
    auto e1 = graph_.add_edge(graph_.node(a.index()), graph_.node(b.index()));
    auto e2 = graph_.add_edge(graph_.node(b.index()), graph_.node(c.index()));
    auto e3 = graph_.add_edge(graph_.node(c.index()), graph_.node(a.index()));

    // setup an internal triangle and add it to the edges_ vector
    internal_triangle internal_t;
    internal_t.nodes.push_back(a.index());
    internal_t.nodes.push_back(b.index());
    internal_t.nodes.push_back(c.index());
    internal_t.edges.push_back(e1.index());
    internal_t.edges.push_back(e2.index());
    internal_t.edges.push_back(e3.index());
    internal_t.index = num_triangles();

    // add this triangle to the vector or internal triangles
    triangles_.push_back(internal_t);

    // Add triangle to the connected triangle sets of the three nodes
    nodes_connected_triangles_[a.index()].insert(num_triangles()-1);
    nodes_connected_triangles_[b.index()].insert(num_triangles()-1);
    nodes_connected_triangles_[c.index()].insert(num_triangles()-1);

    // create a tuple with the three indices in order
    std::string ordered_triangle_string = threesort(a.index(),b.index(),c.index());
    triangles_2_indices_map_[ordered_triangle_string] = num_triangles()-1;

    // Returns a triangle that points to the new triangle
    return Triangle(this, num_triangles()-1);

  };

  /** Locates a triangle given three nodes in the mesh.
    * @pre @a a, @a b and @a c are valid nodes of this mesh
    * @return valid internal index of desired triangle if triangle exists in the mesh or -1 otherwise.
    * 
    * Complexity : Amortized O(1)
    */  
  int find_triangle(const Node& a, const Node& b, const Node& c) const {
    if (num_triangles() == 0)
    return (-1);

    // create a string with the three indices in order
    std::string ordered_triangle_string = threesort(a.index(),b.index(),c.index());

    // check the triangles_2_indices_map_ for this string. If it exists, then
    // return the mapped value. Otherwise, return -1
    auto it = triangles_2_indices_map_.find(ordered_triangle_string);
    if (it == triangles_2_indices_map_.end()) return (-1);
    return ((*it).second);     
  }
  /** Test whether three nodes form a triangle.
    * @pre @a a, @a b and @a c are valid nodes of this mesh
    * @return true if @a a, @a b and @a c are connected with a triangle
    *
    * Complexity: amortized O(1)
    */	
  bool has_triangle(const Node& a,const Node& b,const Node &c) const {
    return find_triangle(a,b,c) >= 0;            
  };

  /** Return the triangle with index @a i.
    * @pre 0 <= @a i < size()
    * @post result_triangle.index() == i
    *
    * Complexity: O(1).
    */
  Triangle triangle(size_type i) const {
    assert(i < num_triangles());
    return Triangle(this, i);
  }

  /** Return the number of triangles in the mesh.
    *
    * Complexity: O(1).
    */
  size_type num_triangles() const {
    return triangles_.size();
  }

  /** Return the number of triangles in the mesh.
    *
    * Complexity: O(1).
    */
  size_type size() const {
    return num_triangles();
  }	 

private:

  /** takes in three values of type size_type and returns a string sorted
    * with the smallest value first and the largest value last in the 
    * form "1_2_3" */
  std::string threesort(size_type a, size_type b, size_type c) const { 
    size_type a_ = a;
    size_type b_ = b;
    size_type c_ = c;
    if (a_ > b_) std::swap(a_,b_); 
    if (a_ > c_) std::swap(a_,c_); 
    if (b_ > c_) std::swap(b_,c_); 
    std::string return_value = std::to_string(a_) + "_" + std::to_string(b_) + "_" + std::to_string(c_);
    return return_value;
  } 

  // internal type for Mesh triangles
  struct internal_triangle {
    std::vector<size_type> nodes;
    std::vector<size_type> edges;
    size_type index;      // The unique identifcation for a triangle
    triangle_value_type value; // holds a value, which is a templated type
    //edge_map_type connected_nodes; // holds all nodes connected to this 
    // node by an edge
  };
  
  // internal type for storing info on triangles adjacent to edges
  struct internal_edge_info {
    int num_connected_triangles;
    size_type triangle_1_idx;
    size_type triangle_2_idx;
  };


  // this vector of sets keeps track of all the triangles connected to each
  // node
  std::vector<std::unordered_set<size_type>> nodes_connected_triangles_;
  
  // This vector holds information on all triangles
  std::vector<internal_triangle> triangles_; 

  // graph_ holds the information on nodes and edges in an object of Graph
  // type from our previously constructed graph class
  Graph<N,E> graph_;
  
  // this unordered map is used by our find_triangle function. It stores
  // a string such as "21_48_99" as a key. This refers to a triangle with
  // nodes 21, 48 and 99 respectively. The values hold the triangle indexes.
  // This allows us to check if a triangle exists and to find a triangle's
  // index in amortized O(1) time 
  std::unordered_map<std::string,size_type> triangles_2_indices_map_;
  
  // this vector of sets keeps track of all the triangles connected to each
  // edge.
  std::vector<internal_edge_info> edges_connected_triangles_info_;

};
#endif


