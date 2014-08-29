/**
 * @file shortest_path.cpp
 * Test script for using our templated Graph to determine shortest paths.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"

#include "Graph.hpp"

#include <vector>
#include <fstream>
#include <queue>

/** Comparator that compares the distance from a given point p.
 */
struct MyComparator {
   Point p_;
   MyComparator(const Point& p) : p_(p) {
   };

   template <typename NODE>
   bool operator()(const NODE& node1, const NODE& node2) const {
    double distance1 = norm(node1.position() - p_);
    double distance2 = norm(node2.position() - p_);
    if (distance1 < distance2)
      return true;
    else
      return false;
  }
};


/** Calculate shortest path lengths in @a g from the nearest node to @a point.
 * @param[in,out] g Input graph
 * @param[in] point Point to find the nearest node to.
 * @post Graph has modified node values indicating the minimum path length
 *           to the nearest node to @a point
 * @post Graph nodes that are unreachable to the nearest node to @a point have
 *           the value() -1.
 * @return The maximum path length found.
 *
 * Finds the nearest node to @a point and treats that as the root node for a
 * breadth first search.
 * This sets node's value() to the length of the shortest path to
 * the root node. The root's value() is 0. Nodes unreachable from
 * the root have value() -1.
 */
int shortest_path_lengths(Graph<int>& g, const Point& point) {
  // TYPE DEFINIATIONS
  typedef Graph<int> GraphType;
  typedef GraphType::node_type Node;
  typedef std::pair <Node,int> Pair;

  // apply find closest node
  auto root = std::min_element(g.node_begin(), g.node_end(), MyComparator(point));

  // initiate queue for BFS
  std::queue<Pair> bfs;

  // initialize default distance from point and push node and distance value in queue
  int distance = 1;
  bfs.push(Pair(*root, distance));

  // initialize pair variables for the queue
  Pair node_distance;
  Pair next;

  // Dequeue and add adjacent nodes not yet searched to the queue while updating distance
  while(!bfs.empty()){
    node_distance = bfs.front();
    for (auto it = node_distance.first.edge_begin(); it != node_distance.first.edge_end(); ++it){
      if ((*it).node1() != node_distance.first)
        next.first = (*it).node1();
      else
        next.first = (*it).node2();

      if (next.first.value() == 0){
        distance = node_distance.second;
        next.first.value() = distance;
        next.second = ++distance;
        bfs.push(next);
      }
    }
    bfs.pop();
  }

  // change the value of nodes not found to -1 and update root node value to 0
  for (auto it = g.node_begin(); it != g.node_end(); ++it){
    if ((*it).value() == 0)
      (*it).value() = -1;
  }
  (*root).value() = 0;

  return distance;
}

// calls the color map
struct Color {
  int d_;
  Color(int d) : d_(d) {}
  CS207::Color operator()(const Graph<int>::node_type n) {
    return CS207::Color::make_heat(static_cast<float> (n.value())/d_);
  }
};

int main(int argc, char* argv[])
{
  // Check arguments
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a Graph
  typedef Graph<int> GraphType;
  GraphType graph;
  std::vector<GraphType::node_type> nodes;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interprit each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  while (CS207::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interprit each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CS207::getline_parsed(tets_file, t))
    for (unsigned i = 1; i < t.size(); ++i)
      for (unsigned j = 0; j < i; ++j)
        graph.add_edge(nodes[t[i]], nodes[t[j]]);

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  viewer.launch();

  // Use shortest_path_lengths to set the node values to the path lengths
  // Construct a Color functor and view with the SDLViewer
  auto node_map = viewer.empty_node_map(graph);
  int distance = shortest_path_lengths(graph, {-1,0,1});
  viewer.add_nodes(graph.node_begin(), graph.node_end(), Color(distance), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  return 0;
}
