#ifndef CS207_GRAPH_HPP
#define CS207_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include "CS207/Util.hpp"
#include "Point.hpp"

#include <algorithm>
#include <vector>
#include <cassert>

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {
 private:
  // Internal type for set elements
  struct internal_node;

 public:

  // PUBLIC TYPE DEFINITIONS

  /** Type of this graph. */
  typedef Graph graph_type;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  typedef Node node_type;
  /** Supports the user-specified value for the node*/
  typedef V node_value_type;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  typedef Edge edge_type;
  /** Supports the user-specified value for the edge*/
  typedef E edge_value_type;

  /** Type of indexes and sizes. Return type of Node::index() and
      Graph::num_nodes(), argument type of Graph::node. */
  typedef unsigned size_type;

  /** Type of node iterators, which iterate over all graph nodes. */
  class node_iterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class edge_iterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class incident_iterator;

  // CONSTRUCTOR AND DESTRUCTOR

  /** Construct an empty graph. */
  Graph():edge_counter_(0){
  }
  /** Default destructor */
  ~Graph() = default;

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
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
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
      return graph_->nodes_[uid_].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_->nodes_[uid_].index;
    }

    /** Return this node's uid, a number in the range [0, total number of nodes added). */
    size_type uid() const{
      return uid_;
    }

    /** Return this node's value, 0 by default */
    node_value_type& value(){
      return graph_->nodes_[uid_].value;
    }

    /** Return this node's value, 0 by default */
    const node_value_type& value() const{
      return graph_->nodes_[uid_].value;
    }

    /** compares the 2 nodes by checking the indices
    *   @pre Node n is a valid node in the nodes_ container
    */
    bool operator< (const Node& n) const{
      return index() < n.index();
    }

    /** checks if 2 nodes are equal by comparing the positions and their respective graph
    *   @pre Node n is a valid node in the nodes_ container
    */
    bool operator==(const Node& n) const{
      return uid() == n.uid() && graph_ == n.graph_;
    }

    /** Return the number of edges incident to the node
    *   Time complexity: O(1)
    */
    size_type degree() const{
      return graph_->edges_[uid()].size();
    }

    /** Return incident interator at the beginning by passing in node uid and 0th index   
    *   Time complextity O(1)
    */  
    incident_iterator edge_begin() const{
      return incident_iterator(graph_, uid(), 0);
    }

    /** Return incident interator at the beginning by passing in node uid and last index   
    *   Time complextity O(1)
    */  
    incident_iterator edge_end() const{
      return incident_iterator(graph_, uid(), degree());
    }

    /** Sets the node's position   
    *   Time complextity O(1)
    */  
    void set_position(const Point& p){
      graph_->nodes_[uid_].position = p;
    }

   private:
    // Only Graph can access our private members
    friend class Graph;
    // Use this space declare private data members for Node objects
    // Graph needs a way to construct valid Node objects

    // Pointer back to the graph container  
    graph_type* graph_;
    // This node's unique id number
    size_type uid_;

    /** Private Constructor */
    Node(const graph_type* graph, size_type uid)
      : graph_(const_cast<graph_type*>(graph)), uid_(uid) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return node_i2u_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new size() == old size() + 1
   * @post result_node.index() == old size()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    node_i2u_.push_back(nodes_.size());
    nodes_.push_back({position, value, node_i2u_.size()-1});
    edges_.resize(edges_.size() + 1);
    return Node(this, nodes_.size()-1);
  }

  /** Determine if this Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (n.graph_ == this && n.index() < num_nodes());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < size()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < size());
    return Node(this, node_i2u_[i]);
  }

  /** Remove a node from the graph and edges corresponding to it 
   * @param[in] n Node to be removed
   * @pre @a n is a valid node of this graph.
   * @post new size() == old size() - 1
   * @post index value in nodes container is updated
   * @post edge container is updated
   *
   * Can invalidate outstanding iterators. @a n becomes invalid, as do any
   * other Node objects equal to @a n. All other Node objects remain valid.
   *
   * Complexity: Polynomial in size().-O(number of nodes added + E)
   */
  void remove_node(const Node& n) {
    node_i2u_.erase(node_i2u_.begin() + nodes_[n.uid()].index);
    for (auto it = nodes_.begin() + n.uid(); it != nodes_.end(); ++it)
      --it->index;

    if (n.uid() < edges_.size()){
      for (auto it1 = edges_[n.uid()].begin(); it1 != edges_[n.uid()].end(); ++it1){
        for (auto it2 = edges_[*it1].begin(); it2 != edges_[*it1].end(); ++it2){
          if (*it2 == n.uid()){
            edges_[*it1].erase(it2);
            break;
          }
        }
      }
      edges_[n.uid()].clear();
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    node_i2u_.clear();
    edges_.clear();
    edge_counter_ = 0;
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
      return graph_->node(node1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->node(node2_);
    }

    /** Return this edges's length */
    double length() const{
      return norm(node1().position() - node2().position());
    }

    /** Return this edges's value, (1,0) by default, O(d)*/
    edge_value_type& value(){
      // ensures that edge(a,b) and edge(b,a) points to the same value type by sorting them
      size_type bigger = (node1_ > node2_) ? node1_ : node2_;
      size_type smaller = (node1_ < node2_) ? node1_ : node2_;

      for (auto it = graph_->edges_[bigger].begin(); it != graph_->edges_[bigger].end(); ++it){
        if ((*it).node_uid == smaller)
          return (*it).value;
      }
      return graph_->edges_[node1_][node2_].value;
    }

    /** Return this edges's value, (1,0) by default, O(d)*/
    const edge_value_type& value() const{
      // ensures that edge(a,b) and edge(b,a) points to the same value type by sorting them
      size_type bigger = (node1_ > node2_) ? node1_ : node2_;
      size_type smaller = (node1_ < node2_) ? node1_ : node2_;
      
      for (auto it = graph_->edges_[bigger].begin(); it != graph_->edges_[bigger].end(); ++it){
        if ((*it).node_uid == smaller)
          return (*it).value;
      }
      return graph_->edges_[node1_][node2_].value;
    }

    /** compares 2 edges by comparing its nodes and its respective graph
    * @pre Edge @a e is a valid edge in the graph
    */
    bool operator==(const Edge& e) const{
      return (graph_ == e.graph_) && 
             ((node1() == e.node1() && node2() == e.node2()) || 
              (node2() == e.node1() && node1() == e.node2()));
    }

    /** compares 2 edges by comparing its nodes and respective graph
    * @pre Edge @a e is a valid edge in the graph
    */
    bool operator<(const Edge& e) const{
      return (graph_ != e.graph_) || (node1() < e.node1());
    }

   private:
    // Only Graph can access our private members
    friend class Graph;
    // Use this space declare private data members for Edge objects
    // Graph needs a way to construct valid Edge objects
    
    // Pointer back to the graph container  
    graph_type* graph_;
    // The edge's 2 node indices
    size_type node1_;
    size_type node2_;

    /** Private Constructor with index value*/
    Edge(const graph_type* graph, size_type node1, size_type node2)
        : graph_(const_cast<graph_type*>(graph)), node1_(node1), node2_(node2){
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less, O(1)
   */
  size_type num_edges() const {
    return edge_counter_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less, O(E)
   */
  Edge edge(size_type i) const {
    assert (i < num_edges());
    edge_iterator ei = edge_begin();
    for (size_type j = 0; j < i; j++)
      ++ei;
    return *ei;
  }

  Edge edge(size_type i, size_type j) const {
    return Edge(this, i, j);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return true if, for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less, O(adjacent nodes)
   */
  bool has_edge(const Node& a, const Node& b) const {
    if (a.uid() >= edges_.size() || b.uid() >= edges_.size())
      return false;
  
    for (auto it = edges_[a.uid()].begin(); it != edges_[a.uid()].end(); ++it){
      if ((*it).node_uid == b.uid())
        return true;
    }

    return false;
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
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less, O(adjacent nodes)
   */
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = edge_value_type()){
    if (has_edge(a, b))
      return Edge(this, a.uid(), b.uid());

    // resize the edges_ vector to fit new nodes if necessary
    size_type bigger = (a.uid() > b.uid()) ? a.uid() : b.uid();
    if (bigger >= edges_.size())
      edges_.resize(bigger+1);

    edges_[a.uid()].push_back({b.uid(), value});
    edges_[b.uid()].push_back({a.uid(), value});
    ++edge_counter_;

    return Edge(this, a.uid(), b.uid());
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
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less O(adjacent nodes*2)
   */
  size_type remove_edge(const Node& a, const Node& b) {
    if (a.uid() >= edges_.size() || b.uid() >= edges_.size()){
      return 0;
    }

    for (auto it = edges_[a.uid()].begin(); it != edges_[a.uid()].end(); ++it){
      if (*it == b.uid()){
        edges_[a.uid()].erase(it);
        break;
      }
    }

    for (auto it = edges_[b.uid()].begin(); it != edges_[b.uid()].end(); ++it){
      if (*it == a.uid()){
        edges_[b.uid()].erase(it);
        --edge_counter_;
        return 1;
      }
    }

    return 0;
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
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less, O(adjacent nodes*2)
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(node(e.node1(), e.node2()));
  }

  // ITERATORS

  /** @class Graph::node_iterator
   * @brief Iterator class for nodes. A forward iterator. */
  class node_iterator : private totally_ordered<node_iterator>{
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
    /** Type of the node vector iterator type*/
    typedef std::vector<size_type>::const_iterator iterator_type;

    /** Construct an invalid node_iterator. */
    node_iterator() {
    }

    /** Returns the node the iterator is at */
    Node operator*() const{
      return graph_->node(*it_);
    }

    /** Iterates to next node if not at end and return this node_iterator */ 
    node_iterator& operator++(){
      if (it_ != graph_->node_i2u_.end())
        ++it_;
      return *this;
    }

    /** Checks if iterators are equal */ 
    bool operator==(const node_iterator& ni) const{
      return it_ == ni.it_;
    }

   private:
    friend class Graph;
    friend class edge_iterator;
    // Pointer back to the graph container  
    graph_type* graph_;
    // The iterator to the node_i2u container
    iterator_type it_;

    // private constructor
    node_iterator(const graph_type* graph, iterator_type it)
      : graph_(const_cast<graph_type*>(graph)), it_(it){
    }
  };

  /** Return node interator at the beginning by passing in beginning iterator of container node_i2u 
  *   Time complextity O(1)
  */  
  node_iterator node_begin() const{
    return node_iterator(this, node_i2u_.begin());
  }

  /** Return node interator at the end by passing in end iterator of container node_i2u 
  *   Time complextity O(1)
  */  
  node_iterator node_end() const{
    return node_iterator(this, node_i2u_.end());
  }


  /** @class Graph::edge_iterator
   * @brief Iterator class for edges. A forward iterator. */
  class edge_iterator : private totally_ordered<edge_iterator>{
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

    /** Returns the edge that the iterator is at 
    *   Ensures that edge is called once by putting the paramters in acsending order
    */
    Edge operator*() const{
      return Edge(graph_, graph_->edges_[uid_][index_].node_uid, uid_);
    }

    /** Iterates to next Edge if not at end and return this edge_iterator */
    edge_iterator& operator++(){
      if (uid_ == graph_->edges_.size()-1 && index_ == graph_->edges_[uid_].size())
        return *this;

      do {
        if (index_ != graph_->edges_[uid_].size()-1)
          ++index_;
        else{
          ++uid_;
          index_ = 0;
        }
      } while (uid_ < graph_->edges_[uid_][index_].node_uid);
      return *this;
    }

    /** Checks if edge's 2 node uids are equal */ 
    bool operator==(const edge_iterator& ei) const{
      return uid_ == ei.uid_ && index_ == ei.index_;
    }

   private:
    friend class Graph;
    // Pointer back to the graph container  
    graph_type* graph_;
    // The iterator to the outer vector of the edges_ container, which is indexed by the node's uid
    size_type uid_;
    // The iterator to the inner vector of the edges_ container, where values are the node's uid
    size_type index_;

    // private constructor
    edge_iterator(const graph_type* graph, size_type uid, size_type index)
      : graph_(const_cast<graph_type*>(graph)), uid_(uid), index_(index){
    }
  };

  /** Return edge interator at the beginning by passing in beginning iterator of container edges_
  *   Time complextity O(1)
  */ 
  edge_iterator edge_begin() const{
    return edge_iterator(this, 0, 0);
  }

  /** Return edge interator at the end by passing in end iterator of container edges_
  *   Time complextity O(1)
  */ 
  edge_iterator edge_end() const{
    return edge_iterator(this, edges_.size()-1 ,edges_[size()-1].size()-1);
  }


  /** @class Graph::incident_iterator
   * @brief Iterator class for edges incident to a given node. A forward
   * iterator. */
  class incident_iterator : private totally_ordered<incident_iterator>{
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

    //** Returns the edge that the iterator is at */
    Edge operator*() const{
      return Edge(graph_, uid_, graph_->edges_[uid_][index_].node_uid);
    }

    /** Iterates to next Edge if not at end and return this edge_iterator */
    incident_iterator& operator++(){
      if (index_ < graph_->edges_[uid_].size())
        ++index_;
      return *this;
    }
    
    /** Checks if edge's 2 uid and index are equal */ 
    bool operator==(const incident_iterator& ii) const{
      return uid_ == ii.uid_ && index_ == ii.index_;
    }

   private:
    friend class Graph;
    // Pointer back to the graph container
    graph_type* graph_;
    // uid of incident node
    size_type uid_;
    // index to adjacent node
    size_type index_;

    // private constructor
    incident_iterator(const graph_type* graph, size_type uid, size_type index)
      : graph_(const_cast<graph_type*>(graph)), uid_(uid), index_(index){
    }
  };

 private:

  // Use this space for your Graph class's internals:
  //   helper functions you might need, data members, and so forth.

  // node data structure
  struct internal_node{
    Point position;         // The node's position
    node_value_type value;  // The optional node value
    size_type index;        // The node's index in the i2u container
  };

  // edge data structure
  struct internal_edge{
    size_type node_uid;     // The edge's node2 uid
    edge_value_type value;  // the optional edge value
  };

  // container for nodes indexed by node's uid
  std::vector<internal_node> nodes_;

  // container for nodes indices and value is node's uid
  std::vector<size_type> node_i2u_;

  /* container for edges: adjancy list of nodes in edges
  // each node is listed twice-ie edges[a] has b and edges_[b] has a
  */
  std::vector< std::vector<internal_edge> > edges_;

  // keeps track of the number of edges
  size_type edge_counter_;
};

#endif
