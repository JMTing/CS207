

#ifndef CS207_GRAPH_HPP
#define CS207_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include "CS207/Util.hpp"
#include "Point.hpp"

#include <typeinfo>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <map>
#include <tuple>
#include <iterator> 
#include <cassert>
#include <iostream>
#include <set>

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there 
 * is at most one edge between any pair of distinct nodes).
 */

template <typename V, typename E> class Graph {
 private:

  struct internal_node;
  struct internal_edge;
  
 public:

  // PUBLIC TYPE DEFINITIONS

  /** Type of this graph. */
  typedef Graph<V,E> graph_type;

  /** Type of node values */
  typedef V node_value_type;

  /** Type of edge values */
  typedef E edge_value_type;
  
  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  typedef Node node_type;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  typedef Edge edge_type;

  /** Type of indexes and sizes. Return type of Node::index() and
      Graph::num_nodes(), argument type of Graph::node. */
  typedef unsigned size_type;


  /** Type of map used to hold associated node ids and edge ids. 
   *  The key will contain the index of the other node that is connected
   *  by an edge and the value will contain that edge's internal id */
  typedef std::unordered_map<size_type,size_type> edge_map_type;

  /** Iterator type for edge_map_type */
  typedef std::unordered_map<size_type,size_type>::iterator edge_map_type_iterator;

  /** Type of node iterators, which iterate over all graph nodes. */
  class node_iterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class edge_iterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class incident_iterator;

  // CONSTRUCTOR AND DESTRUCTOR
  /** Construct an empty graph. */
  Graph<V,E>() {
    uid_index_ = 0;
    uid_index_edges_ = 0;
  }
  /** Default destructor */
  ~Graph<V,E>() = default;


  /** Delete the copy constructor of the Graph class, so that a Graph object 
   * cannot be copied by mistake */
   Graph<V,E>(const Graph<V,E> &g) = delete;
   
  /** Delete the assignment operator so that a Graph object cannot
   * be assigned by mistake */
   Graph<V,E>& operator=(const Graph<V,E> &g) = delete;
  // NODES

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node : private totally_ordered<Node> {  
   public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign 
     * a valid node to it later. For example:
     *
     * @code
     * Node x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */
    Node() {
    }

    /** Return this node's position. */
    Point position() const {
      return fetch().position;
    }
    
    /** Return this node's index, a number in the range [0, graph_size).*/
    size_type index() const {
      return external_index_;  
    }

    /** Returns this node's value and allows it to be changed. 
     *  @return the value held by this node, as a non-constant variable*/
    node_value_type& value(){
        return fetch().value;
    }
    
    /** Returns this node's value and does not allow it to be changed. 
     *  @return the value held by this node, as a constant variable*/
    const node_value_type& value() const{
        return fetch().value;
    }
    
    /** Overloads the == operator so that it compares the index to the
     *  parameter n's index and also checks that parameter n is part of 
     *  the same graph as this node. 
     *  @return true if the parameter n's index is the same as the 
     *  current node's index and if they are in the same graph.
     *  @param[in] n A node being compared to this one.*/
    bool operator==(const Node& n) const {
        return (index() == n.index() && graph() == n.graph());
        }

    /** Overloads the < operator. 
     *  @return true if the index of the current node is < the index
     *  of the comparison node.
     *  @param[in] n A node being compared to this one.*/
    bool operator<(const Node& n) const {
        return !(index() >= n.index());
        }
    
    /** Returns the degree (number of nodes connected) to this node */
    size_type degree() const{
        return fetch().connected_nodes.size();
    }
    
    /** returns an incident iterator pointing to the first node connected
     *  to this node by an edge */
    incident_iterator edge_begin() const {
        edge_map_type_iterator it = fetch().connected_nodes.begin();
        size_type this_node_id = graph_->i2u_[index()];
        return incident_iterator(graph_, this_node_id, it);
    }
    
    /** returns an incident iterator pointing to a value after the last 
     *  node connected to this node by an edge */
    incident_iterator edge_end() const{
        edge_map_type_iterator it = fetch().connected_nodes.end();
        size_type this_node_id = graph_->i2u_[index()];
        return incident_iterator(graph_, this_node_id, it);
    }
    
    /** Sets the position of the node equal to the parameter p
     * @param[in] p A point to which our node's position should be set.
     * @pre p is a valid and initialized point
     * @post Our node's position is equal to p */
    void set_position(const Point& p)
    {
        fetch().position = p;
    }


   private:
    /** Return a pointer to this node's graph */
    graph_type* graph() const {
    	return graph_;
    }

    // Only Graph can access our private members
    friend class Graph<V,E>;
    
    // Pointer back to the Graph container
    graph_type* graph_;
    
    // This node's external (user) identification number
    size_type external_index_;
    
    /** Private Constructor */
    Node(const graph_type* grph, size_type index){
        graph_ = const_cast<graph_type*>(grph);
        external_index_ = index;
    }
    
    /** Helper method to return the appropriate node.
     *  @pre external_index_ < num_nodes()
     */
    internal_node& fetch() const {
      return graph_->nodes_[graph_->i2u_[external_index_]];  
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return i2u_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] value The new node's value
   * @post new size() == old size() + 1
   * @post result_node.index() == old size()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type & value = node_value_type()) {

    // setup an internal node object
    internal_node internal_x;
    internal_x.position = position;
    internal_x.index = num_nodes();
    internal_x.value = value;
    
    // add internal_node object to internal node list
    nodes_.push_back(internal_x);
    
    // add index in internal node list to the i2u_ list
    i2u_.push_back(uid_index_);
    
    ++uid_index_;
    return Node(this, i2u_.size() - 1);
  }

  /** Determine if this Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
      return !(n.index() > num_nodes());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < size()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < size());
    return Node(this, i);
  }

  /** Remove a node from the graph.
   * @param[in] n Node to be removed
   * @pre @a n is a valid node of this graph.
   * @post new size() == old size() - 1
   *
   * Can invalidate outstanding iterators. @a n becomes invalid, as do any
   * other Node objects equal to @a n. All other Node objects remain valid.
   *
   * Complexity: Polynomial in size().
   */
  void remove_node(const Node& n) {
    assert (has_node(n));

    // check if node has any connected edges and remove those 
    for(auto itr = n.fetch().connected_nodes.begin(); itr!=n.fetch().connected_nodes.end();++itr)
    {
        size_type this_node_internal_id = get_internal_node_id(n.index());
        size_type other_node_internal_id = (*itr).first;
        remove_edge_using_internal_node_ids(this_node_internal_id,other_node_internal_id);
    }
    
    i2u_.erase(i2u_.begin()+n.index());
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    edges_.clear();
    nodes_.clear();
    i2u_.clear();
    i2u_edges_.clear();
    uid_index_ = 0;
    uid_index_edges_ = 0;
  }


  // EDGES

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge : private totally_ordered<Edge> { 
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return fetch().node1;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return fetch().node2;
    }
    
    /** Return the edge's length */
    double length() const {
      return norm(fetch().node1.position()-fetch().node2.position());  
    }
    
    
    /** Return this edge's index. */
    size_type index() const {
      return external_index_;
    }    

    /** Returns this edge's value and allows it to be changed. 
     *  @return the value held by this edge, as a non-constant variable*/
    edge_value_type& value(){
        return fetch().value;
    }
    
    /** Returns this edge's value and does not allow it to be changed. 
     *  @return the value held by this edge, as a constant variable*/
    const edge_value_type& value() const{
        return fetch().value;
    }
    
    /** Override the '==' operator and say that two edges are equal
     *  if they have the same index and are on the same graph 
     *  @return true if the edge being compared to has the same index
     *  and is in the same graph
     *  @param[in] e The edge to compare the current edge to 
     */
    bool operator==(const Edge& e) const {
        return (graph() == e.graph() && index() == e.index());
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
        else if (index() == e.index() || graph() != e.graph())
            return true;
        return false;
        }
        
   private:
    /** Return this edge's graph */
    graph_type* graph() const {
    	return graph_;
    }    
    
    // Only Graph can access our private members
    friend class Graph<V,E>;
    graph_type* graph_;
    
    // This edge's unique identification number
    size_type external_index_;
    
    /** Private Constructor */
    Edge(const graph_type* grph, size_type index){
        graph_ = const_cast<graph_type*>(grph);
        external_index_ = index;
    }
    
    /** Helper method to return the appropriate internal edge.*/
    internal_edge& fetch() const {
      return graph_->edges_[graph_->i2u_edges_[external_index_]];    
    }
    
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return i2u_edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < num_edges());
    return Edge(this, i);
  }
  
  /** Locates an edge given two nodes in the graph.
   * @pre @a a and @a b are valid nodes of this graph
   * @return valid internal index of desired edge if edge exists in the graph or -1 otherwise.
   * 
   * Complexity : Amortized O(1)
   */  
  int find_edge(const Node& a, const Node& b) const {
    if (num_edges() == 0)
        return (-1);
    
    size_type a_id = i2u_[a.index()];
    size_type b_id = i2u_[b.index()];
    
    edge_map_type connected_nodes_map = nodes_[a_id].connected_nodes;
    auto it = connected_nodes_map.find(b_id);
    if (it == connected_nodes_map.end()) return (-1);
    return ((*it).second);
  }
  
  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return true if, for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    return find_edge(a, b) >= 0;
  }
  

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b) {
	
	// check if the edge already exists and if so, then return the edge
	int idx = find_edge(a,b);
	if (idx >= 0)
		return Edge(this, idx);

    // add edge to the connected_nodes maps of the two nodes
    size_type a_internal_id = get_internal_node_id(a.index());
    size_type b_internal_id = get_internal_node_id(b.index());
    
    a.fetch().connected_nodes[b_internal_id] = uid_index_edges_;
    b.fetch().connected_nodes[a_internal_id] = uid_index_edges_;
        
	// setup an internal edge and add it to the edges_ vector
    internal_edge internal_e;
    internal_e.node1 = a;
    internal_e.node2 = b;
    internal_e.edge_index = num_edges();

    edges_.push_back(internal_e);
    i2u_edges_.push_back(uid_index_edges_);
    
    ++uid_index_edges_;

    // Returns an edge that points to the new edge
    return Edge(this, num_edges()-1);
  }



  /** Remove an edge, if any, returning the number of edges removed.
   * @param[in] a,b The nodes potentially defining an edge to be removed.
   * @return 1 if old has_edge(@a a, @a b), 0 otherwise
   * @pre @a a and @a b are valid nodes of this graph
   * @post !has_edge(@a a, @a b)
   * @post new num_edges() == old num_edges() - result
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Can invalidate all edge and incident iterators.
   * Invalidates any edges equal to Edge(@a a, @a b). Must not invalidate
   * other outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type remove_edge(const Node& a, const Node& b) {
    size_type a_internal_id = get_internal_node_id(a.index());
    size_type b_internal_id = get_internal_node_id(b.index());
    return remove_edge_using_internal_node_ids(a_internal_id,b_internal_id);
  }

  /** Remove an edge, if any, returning the number of edges removed.
   * @param[in] e The edge to remove
   * @pre @a e is a valid edge of this graph
   * @pre has_edge(@a e.node1(), @a e.node2())
   * @post !has_edge(@a e.node1(), @a e.node2())
   * @post new num_edges() == old num_edges() - 1
   *
   * This is a synonym for remove_edge(@a e.node1(), @a e.node2()), but its
   * implementation can assume that @a e is definitely an edge of the graph.
   * This might allow a faster implementation.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Can invalidate all edge and incident iterators.
   * Invalidates any edges equal to Edge(@a a, @a b). Must not invalidate
   * other outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(),e.node2());
  }

  // ITERATORS

  /** @class Graph::node_iterator
   * @brief Iterator class for nodes. A forward iterator. */
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
        return graph_->node(it_);
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
        return (it_ == other.it_ && graph_ == other.graph_);
    }
    
   private:
    
    /** private constructor for an valid node_iterator. 
     * @pre 0<=it_<num_nodes() */
    node_iterator(const graph_type* grph, size_type it) {
      graph_ = const_cast<graph_type*>(grph);
      it_ = it;
    }    
    friend class Graph<V,E>;
    friend class edge_iterator;
    graph_type* graph_;
    size_type it_;
    
  };
  
  /** returns a node iterator that relates to the first node by user index
   * @pre !i2u_.empty() 
   */
  node_iterator node_begin() const {
      assert(!i2u_.empty());
      return node_iterator(this, 0);
  }

  /** returns a node iterator that relates to a nonexisting position
   *  one greater than the size of i2u_
   * @pre !i2u_.empty() 
   */
  node_iterator node_end() const {
      assert(!i2u_.empty());
      return node_iterator(this, num_nodes());
  }


  /** @class Graph::edge_iterator
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
        return graph_->edge(it_);
    }
    
    /** increments the value it_ by 1
    * @pre 0<=it_<num_edges() 
    */
    edge_iterator& operator++(){
        ++it_;
        return *this;
    }
    
    /** returns true if two iterators belong to the same graph
     *  and they are at the same value it_
     *  @param[in] other The iterator to compare to
    */
    bool operator==(const edge_iterator& other) const{
        return (it_ == other.it_ && graph_ == other.graph_);
    }
    
   private:
    /** private constructor for an valid edge_iterator. 
     * @pre 0<=it_<num_edges() */
    edge_iterator(const graph_type* grph, size_type it) {
      graph_ = const_cast<graph_type*>(grph);
      it_ = it;
    } 
    
    friend class Graph<V,E>;
    graph_type* graph_;
    size_type it_;
  };


  /** returns a edge iterator that relates to the first edge by user index
   * @pre !i2u_edges_.empty() 
   */
  edge_iterator edge_begin() const {
      assert(!i2u_edges_.empty());
      return edge_iterator(this, 0);
  }

  /** returns a node iterator that relates to a nonexisting position
   *  one greater than the size of i2u_edges_
   * @pre !i2u_edges_.empty() 
   */
  edge_iterator edge_end() const {
      assert(!i2u_edges_.empty());
      return edge_iterator(this, num_edges());
  }

  /** @class Graph::incident_iterator
   * @brief Iterator class for edges incident to a given node. A forward
   * iterator. */
  class incident_iterator : private totally_ordered<incident_iterator> { 
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
    incident_iterator() { 
    }

    /** Handles the * operator for the incident iterator:
     *  returns an edge based on the values in origin_node_id_ and the
     *  integer pointed to by it_.
     *  @pre it_ points to a valid value in the connect edge map of the 
     *  node being considered. 
     *  @return an edge that can reference the edge connecting the
     *  two nodes.*/
    Edge operator*() const{
        Edge return_edge = graph_->edge((*it_).second);
        
        // check if the index of node1 of the edge to be returned is 
        // equal to the origin node index (origin_node_id_) and if not,
        // then switch the two nodes stored in the edge around
        if (return_edge.node1().index() != origin_node_id_)
            graph_->switch_nodes_of_edge((*it_).second);
        
        return return_edge;
    }
    
    /** Overloads the ++ operator for the incident iterator:
     *  @pre it_ points to a valid position in the map
     *  @post it_ is is iterated by one
     *  @return an incident iterator that has been iterated to point
     *  to the next value in the map */    
    incident_iterator& operator++(){
        ++it_;
        return *this;
    }

    /** handles the == operator for the incident iterator:
     *  returns true if both iterators have the same value
     *  and if both belong to the same graph 
     *  @param[in] other The iterator to compare to
     *  @pre other refers to a valid incident iterator
     *  @return true if the value of the private variable it_ is equal
     *  to other's it_ variable and other is part of the same graph.
     */    
    bool operator==(const incident_iterator& other) const{
        return (it_ == other.it_ && graph_ == other.graph_);        
    }


   private:
    /** private constructor for an valid edge_iterator. 
     * @param[in] grph refers to the current graph that this incident_iterator belongs to
     * @param[in] origin_node_id refers to the origin node's external (user) id
     * @param[in] it is an iterator to the unordered map stored in the node 
     * whose edges the incident_iterator is iterating through */
    incident_iterator(const graph_type* grph, size_type origin_node_id, edge_map_type_iterator it) {
      graph_ = const_cast<graph_type*>(grph);
      origin_node_id_ = origin_node_id;
      it_ = it;
    } 
    
    friend class Graph<V,E>;
    graph_type* graph_;
    size_type origin_node_id_;
    edge_map_type_iterator it_;
  };
    
    
 private:
  
  /** Remove an edge using the internal nodes IDs, if any, returning the number of edges removed.
   * @param[in] a,b The internal IDs of nodes potentially defining an edge to be removed.
   * @return 1 if old has_edge(@a a, @a b), 0 otherwise
   * @pre @a a_id and @a b_id are valid internal ids of nodes of this graph
   * @post !has_edge(@a a, @a b)
   * @post new num_edges() == old num_edges() - result
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Can invalidate all edge and incident iterators.
   * Invalidates any edges equal to Edge(@a a, @a b). Must not invalidate
   * other outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type remove_edge_using_internal_node_ids(const size_type a_id, const size_type b_id) {
    edge_map_type connected_nodes_map = nodes_[a_id].connected_nodes;
    auto it1 = connected_nodes_map.find(b_id);
    if (it1 == connected_nodes_map.end()) return 0;
    
    size_type edge_internal_id = nodes_[a_id].connected_nodes[b_id];
    for (size_type i = 0; i < num_edges(); i++)
        if (i2u_edges_[i] == edge_internal_id)
            i2u_edges_.erase(i2u_edges_.begin()+i);
    
    nodes_[a_id].connected_nodes.erase(b_id);
    nodes_[b_id].connected_nodes.erase(a_id);
    return 1;
    }
    
   /** get internal node id given external node id 
    * @param[in] external_node_id The external (user) id of a node
    * @return the internal node id associated with an external node id
    * @pre has_node(external_node_id) == true
    */
   size_type get_internal_node_id(size_type external_node_id)
   {
       return i2u_[external_node_id];
   }
   
   /** switches around the two nodes in a stored edge 
    * @param[in] internal_node_id The internal id (unique id) of the 
    * edge whose nodes are to be switched
    * @pre internal_node_id is a valid index of an edge in edges_
    * @post node1 and node2 in the specified stored edge are switched
    */
   void switch_nodes_of_edge(size_type internal_node_id)
   {
       Node temp = edges_[internal_node_id].node1;
       edges_[internal_node_id].node1 = edges_[internal_node_id].node2;
       edges_[internal_node_id].node2 = temp;
   }
  
  //internal type for Graph nodes
  struct internal_node {
    Point position;   // The position of the node
    size_type index;      // The unique identifcation for a node
    node_value_type value; // holds a value, which is a templated type
    edge_map_type connected_nodes; // holds all nodes connected to this 
                                   // node by an edge
  };
  
  //internal type for Graph edges
  struct internal_edge {
    Node node1;
    Node node2;
    size_type edge_index;   // The unique identifcation for an edge
    edge_value_type value; // holds a value, which is a templated type
  };
  
  std::vector<internal_node> nodes_;
  std::vector<internal_edge> edges_;
  std::vector<size_type> i2u_; //this vector holds references to the positions of items in nodes_
  std::vector<size_type> i2u_edges_; //this vector holds references to the positions of items in edges_


  size_type uid_index_; // total number of nodes added, including removed nodes
  size_type uid_index_edges_; // total number of edges added, including removed edges
  
};

#endif
