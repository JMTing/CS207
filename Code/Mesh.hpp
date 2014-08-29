/**
 * Jason Ting
 * Joe Pavlisko
 */

#pragma once

/**
 * @file Mesh.hpp
 */

/** @class Mesh
 * @brief A template for 3D triangular meshes.
 *
 * Users can add triangles and retrieve nodes, edges, and triangles.
 */

#include "Graph.hpp"
#include "Point.hpp"
#include <vector>

template <typename N, typename E, typename T>
class Mesh {
  public:

    /** Construct an empty mesh. */
    Mesh() {
    }

    /** Default destructor */
    ~Mesh() = default;

    /** Supports the user-specified value for the triangle*/
    typedef T triangle_value_type;

    /** Type of indices and sizes. */
    typedef int size_type;

    // edge data structure
    struct edge_triangles : E {
        size_type tri1;
        size_type tri2;
        E edge_data;
    };


    /** Synonym for Node (following STL conventions) from graph */
    typedef typename Graph<N,edge_triangles>::node_type Node;
    typedef typename Graph<N,edge_triangles>::node_type node_type;
    typedef typename Graph<N,edge_triangles>::node_iterator node_iterator;
    typedef typename Graph<N,edge_triangles>::edge_iterator edge_iterator;

    /** Predeclaration of triangle type. */
    class Triangle;

    /** Return the number of nodes in the mesh. */
    size_type num_nodes() const {
        return graph_.num_nodes();
    }

    /** Return the number of edges in the mesh. */
    size_type num_edges() const {
        return graph_.num_edges();
    }

    /** Return the number of triangles in the mesh. */
    size_type num_triangles() const {
        return myTriangles_.size();
    }

    /** Adds node to the its graph class */
    Node add_node(Point p) {
   	  node_adjTri.resize(node_adjTri.size()+1);
      return graph_.add_node(p);
    }

  private:
    // triangle data structure
    struct internal_triangle{
      size_type n1_uid;
      size_type n2_uid;
      size_type n3_uid;
      triangle_value_type value;
    };

    // graph to handles nodes and edges
    Graph<N,edge_triangles> graph_;

    // container for triangle indexed by triangle's uid
    std::vector<internal_triangle> myTriangles_;

    // adj list of node's neighboring triangles
    std::vector< std::vector<size_type> > node_adjTri; 

  public:
    class Triangle {
      public:

      	/** Return the triangle's value */
        triangle_value_type& value(){ 
        	return mesh_->myTriangles_[tri_uid_].value;
        }

        /** Return the triangle's value */
        const triangle_value_type value() const{
        	return mesh_->myTriangles_[tri_uid_].value;
        }

        /** Return the triangle's nodes */
        Node node1() {
          return mesh_->graph_.node(mesh_->myTriangles_[tri_uid_].n1_uid);
        }
        Node node2() {
          return mesh_->graph_.node(mesh_->myTriangles_[tri_uid_].n2_uid);
        }
        Node node3() {
          return mesh_->graph_.node(mesh_->myTriangles_[tri_uid_].n3_uid);
        }

        /** Determines if a triangle is valid if uid is a valid index of triangle vector */
        bool is_valid(){
          return tri_uid_ != -1;
        }

        /** Return's the triangle's uid */
        size_type index() {
          return tri_uid_;
        }

        /** Returns the adjacent Triangle given the pair of nodes/edge
          * Returns a triangle with negative uid to indicate its invalid
          * @pre: Nodes form a valid edge of this triangle
          */
        Triangle adj_triangle(Node A, Node B) const {
          auto e = mesh_->graph_.edge(A.uid(), B.uid());
          if (e.value().tri1 == tri_uid_) {
            auto val = e.value().tri2;
            return mesh_->triangle(val);
          }
          return mesh_->triangle(e.value().tri1);
        }

        /** Finds the third node of in a triangle
          * @pre: Nodes form a valid edge of this triangle
          */
        Node opposite(Node& a, Node& b){
        	if ((int(a.index()) == mesh_->myTriangles_[tri_uid_].n1_uid  && 
        		 int(b.index()) == mesh_->myTriangles_[tri_uid_].n2_uid) || 
        		(int(a.index()) == mesh_->myTriangles_[tri_uid_].n2_uid  && 
        		 int(b.index()) == mesh_->myTriangles_[tri_uid_].n1_uid))
        		return mesh_->graph_.node(mesh_->myTriangles_[tri_uid_].n3_uid);

        	else if ((int(a.index()) == mesh_->myTriangles_[tri_uid_].n2_uid  && 
        		      int(b.index()) == mesh_->myTriangles_[tri_uid_].n3_uid) || 
        		     (int(a.index()) == mesh_->myTriangles_[tri_uid_].n3_uid  && 
        		      int(b.index()) == mesh_->myTriangles_[tri_uid_].n2_uid))
        		return mesh_->graph_.node(mesh_->myTriangles_[tri_uid_].n1_uid);

        	else
        		return mesh_->graph_.node(mesh_->myTriangles_[tri_uid_].n2_uid);
        }

        /** Calculates the normal vector of the triangle's edge with node @a and node @b
        * @pre node @a a and node @a b are valid nodes in the triangle
        * normal vector is 2D
        */
    	Point normal(Node a, Node b) {
    		Point e1 = a.position() - b.position();
    		Point e2 = a.position() - opposite(a,b).position();
        Point e3 = Point(e1.y, -e1.x, 0);
        if (dot(e3,e2) < 0)
          return Point(-e3.x, -e3.y, 0);
        return e3;
    	}

        /** Calculate and return the triangle's area */
        float area() {
          return std::abs(node1().position().x*(node2().position().y-node3().position().y) + 
                          node2().position().x*(node3().position().y-node1().position().y) + 
                          node3().position().x*(node1().position().y-node2().position().y))/2.0; 
        }

        /** Calculate and return the triangle's corner */
        Point center(){
          double x = (node1().position().x + node2().position().x + node3().position().x)/3;
          double y = (node1().position().y + node2().position().y + node3().position().y)/3;
          return Point(x,y,0);
        }

      private:
        // Only Mesh can access our private members
        friend class Mesh;

        Mesh* mesh_;
        size_type tri_uid_;

        /** Private Constructor */
        Triangle(Mesh* mesh, size_type uid) : mesh_(const_cast<Mesh*>(mesh)),tri_uid_(uid) {
        }
    };

    /** Return the triangle with index @a i.
    * @pre 0 <= @a i < size()
    * @post result_triangle.uid() == i
    *
    * Complexity: O(1).
    */
    Triangle triangle(size_type i) {
    	assert(i < num_triangles());
    	return Triangle(this, i);
    }

    /** @class Mesh::Triangle_iterator
     * @brief Iterator class for triangles. A forward iterator. */
    class Triangle_iterator {
      public:
        /** Construct an invalid triangle_iterator. */
        Triangle_iterator() {
        }

        // Returns a copy of the Triangle that the iterator is currently pointing to.
        Triangle operator*() const {
          return iterator_parent_mesh->triangle(iterator_index);
        }

        // Iterates to next triangle
        Triangle_iterator& operator++() {
          ++iterator_index;
          return *this;
        }

        // checks it iterators are equal
        bool operator==(const Triangle_iterator& other_iter) const {
          return iterator_index == other_iter.iterator_index && iterator_parent_mesh == other_iter.iterator_parent_mesh;
        }

        // checks it iterators are not equal
        bool operator!=(const Triangle_iterator& other_iter) const {
          return iterator_index != other_iter.iterator_index; //|| !(iterator_parent_mesh == other_iter.iterator_parent_mesh)
        }

       private:
        friend class Mesh;

        Mesh* iterator_parent_mesh;
        size_type iterator_index;

        // private constructor
        Triangle_iterator(const Mesh* mesh, size_type index): iterator_parent_mesh(const_cast<Mesh*>(mesh)), iterator_index(index) {
        }
    };

    /** Return node triangle at the beginning by passing in first index
    *   Time complextity O(1)
    */
    Triangle_iterator triangle_begin() {
      return Triangle_iterator(this, 0);
    }

    /** Return node triangle at the beginning by passing in last index
    *   Time complextity O(1)
    */
    Triangle_iterator triangle_end() {
      return Triangle_iterator(this, myTriangles_.size());
    }

    /** @class Mesh::Node_incident_iterator
     * @brief Iterates through all of the valid triangles adjacent to a given node. */
    class Node_incident_iterator {
      public:
        int get_index() {return index;}

        /** Construct an invalid incident_iterator. */
        Node_incident_iterator() {
        }

        // Returns a copy of the Triangle that the iterator is currently pointing to.
        Triangle operator*() const {
          return parent_mesh->triangle( parent_mesh->node_adjTri[nid][index] );
        }

        // Iterates to next triangle
        Node_incident_iterator& operator++() {
          ++index;
          return *this;
        }

        // checks it iterators are equal
        bool operator==(const Node_incident_iterator& iter2) const {
          return (index == iter2.index) && (parent_mesh->graph_.node(nid) == iter2.parent_mesh->graph_.node(nid));
        }

        // checks it iterators are not equal
        bool operator!=(const Node_incident_iterator& iter2) const {
          return !((index == iter2.index) && (parent_mesh->graph_.node(nid) == iter2.parent_mesh->graph_.node(nid)));
        }

      private:
        friend class Mesh;
        Mesh* parent_mesh;
        size_type index;
        size_type nid;

        // private constructor
        Node_incident_iterator(Mesh* m, size_type it_index, size_type uid) : 
          parent_mesh(const_cast<Mesh*>(m)), index(it_index), nid(uid) {
        }
    };

    /** Return first node_incident_iterator
    *   Time complextity O(1)
    */
    Node_incident_iterator adj_triangle_begin(size_type UID) {
      return Node_incident_iterator(this, 0, UID);
    }

    /** Return last node_incident_iterator (invalid)
    *   Time complextity O(1)
    */
    Node_incident_iterator adj_triangle_end(size_type UID) {
      return Node_incident_iterator(this, node_adjTri[UID].size(), UID);
    }

    /** Add an triangle to the graph and returns the triangle.
      * @pre @a a, @a b, and s@a c are distinct valid nodes of this graph
      * @return an triangle object t with t.node1()==@a a, t.node2()==@a b, and  t.node3()==@a c
      * updates edges in graph
      * updates the adjacent triangle list
      * updates each node's adjacent triangle list
      */
    Triangle add_triangle(const Node& a, const Node& b, const Node& c){
        // updates edge and edge values
    	edge_triangles tris;
    	tris.tri1 = myTriangles_.size();
        tris.tri2 = -1;
    	if (graph_.has_edge(a,b))
    		graph_.edge(a.uid(),b.uid()).value().tri2 = myTriangles_.size();
    	else
    		graph_.add_edge(a,b,tris);

    	if (graph_.has_edge(b,c))
    		graph_.edge(b.uid(),c.uid()).value().tri2 = myTriangles_.size();
    	else
    		graph_.add_edge(b,c,tris);

    	if (graph_.has_edge(c,a))
    		graph_.edge(c.uid(),a.uid()).value().tri2 = myTriangles_.size();
    	else
    		graph_.add_edge(c,a,tris);

    	node_adjTri[a.uid()].push_back(myTriangles_.size());
    	node_adjTri[b.uid()].push_back(myTriangles_.size());
    	node_adjTri[c.uid()].push_back(myTriangles_.size());

        // updates triangle values
    	internal_triangle temp;
    	temp.n1_uid = a.uid();
    	temp.n2_uid = b.uid();
    	temp.n3_uid = c.uid();
        temp.value = triangle_value_type();

    	myTriangles_.push_back(temp);
        return triangle(myTriangles_.size()-1);
    }

  // graph's node iterator begin
  node_iterator node_begin() {
    return graph_.node_begin();
  }

  // graph's node iterator end
  node_iterator node_end() {
    return graph_.node_end();
  }

  // graph's node iterator begin
  edge_iterator edge_begin() {
    return graph_.edge_begin();
  }

  // graph's node iterator end
  edge_iterator edge_end() {
    return graph_.edge_end();
  }

};
