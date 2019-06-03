#ifndef DIRECTED_GRAPH_H
#define DIRECTED_GRAPH_H

//A large selection of data structures from the standard
//library. You need not feel compelled to use them all,
//but as you can't add any, they're all here just in case.
#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <array>
#include <list>
#include <forward_list>
#include <deque>
#include <map>
#include <cstddef>
#include <string>

//Forward declarations for classes below so they can be used below without worrying too much about the ordering.
template <typename vertex> class vertex_iterator;
template <typename vertex> class neighbour_iterator;
template <typename vertex> class directed_graph;


template <typename vertex>
class directed_graph {

private:

  //You will need to add some data members here
  //to actually represent the graph internally,
  //and keep track of whatever you need to.
	std::vector<vertex> vertices; //Vector containing the vertices in the graph
	std::vector<std::vector<int>> adj_matrix; //Matrix containing all the edges e.g. connecting the vertices
	std::size_t graph_size; //The amount of vertices in the graph
	std::size_t number_edges; //Number of Edges of the Graph, removes the need for a loop to check if true in the num_edges method	

public:


  directed_graph(); //A constructor for directed_graph. The graph should start empty.
  ~directed_graph(); //A destructor. Depending on how you do things, this may
  //not be necessary.
  
  int get_index(const vertex&) const; //Returns the index of the specified vertex in graph_vertices for the adj_matrix use e.g. graph_edges.
  
  bool contains(const vertex&) const; //Returns true if the given vertex is in the graph, false otherwise.

  bool adjacent(const vertex&, const vertex&); //Returns true if the first vertex is adjacent to the second, false otherwise.

  void add_vertex(const vertex&); //Adds the passed in vertex to the graph (with no edges).
  void add_edge(const vertex&, const vertex&); //Adds an edge from the first vertex to the second.

  void remove_vertex(const vertex&); //Removes the given vertex. Should also clear any incident edges.
  void remove_edge(const vertex&, const vertex&); //Removes the edge between the two vertices, if it exists.

  std::size_t in_degree(const vertex&); //Returns number of edges coming in to a vertex.
  std::size_t out_degree(const vertex&); //Returns the number of edges leaving a vertex.
  std::size_t degree(const vertex&); //Returns the degree of the vertex (both in and out edges).
  
  std::size_t num_vertices() const; //Returns the total number of vertices in the graph.
  std::size_t num_edges() const; //Returns the total number of edges in the graph.

  std::vector<vertex> get_vertices() const; //Returns a vector containing all the vertices.
  std::vector<vertex> get_neighbours(const vertex&) const; //Returns a vector containing the neighbours of the given vertex.

  vertex_iterator<vertex> begin(); //Returns a graph_iterator pointing to the start of the vertex set.
  vertex_iterator<vertex> end(); //Returns a graph_iterator pointing to one-past-the-end of the vertex set.

  neighbour_iterator<vertex> nbegin(const vertex&); //Returns a neighbour_iterator pointing to the start of the neighbour set for the given vertex.
  neighbour_iterator<vertex> nend(const vertex&); //Returns a neighbour_iterator pointing to one-past-the-end of the neighbour set for the given vertex.

  std::vector<vertex> depth_first(const vertex&); //Returns the vertices of the graph in the order they are visited in by a depth-first traversal starting at the given vertex.
  std::vector<vertex> breadth_first(const vertex&); //Returns the vertices of the graph in the order they are visisted in by a breadth-first traversal starting at the given vertex.

  directed_graph<vertex> out_tree(const vertex&); //Returns a spanning tree of the graph starting at the given vertex using the out-edges.
  directed_graph<vertex> in_tree(const vertex&); //Returns a spanning tree of the graph starting at the given vertex using the in-edges.

  bool reachable(const vertex&, const vertex&); //Returns true if the second vertex is reachable from the first (can you follow a path of out-edges to get from the first to the second?). Returns false otherwise.

};

//The vertex_iterator class provides an iterator
//over the vertices of the graph.
//This is one of the harder parts, so if you're
//not too comfortable with C++ leave this for last.
//If you are, there are many ways of doing this,
//as long as it passes the tests, it's okay.
//You may want to watch the videos on iterators before starting.
template <typename vertex>
class vertex_iterator {

private:

  //You may need data members here.
  const directed_graph<vertex> * new_graph; //Pointer to iterator graph
  std::size_t current_pos; //Tracks current position of iterator

public:
  vertex_iterator(const vertex_iterator<vertex>&); //Copy Constructor
  vertex_iterator(const directed_graph<vertex>&, std::size_t); //Constructor
  ~vertex_iterator(); //Deconstructor
  vertex_iterator<vertex> operator=(const vertex_iterator<vertex>&); //Assignment Operator
  bool operator==(const vertex_iterator<vertex>&) const; //Equals Operator
  bool operator!=(const vertex_iterator<vertex>&) const; //Does not Equals Operator
  vertex_iterator<vertex> operator++(); //Pre-Increment Operator
  vertex_iterator<vertex> operator++(int); //Post-Increment Operator in int form
  vertex operator*(); //Indirection Operator
  vertex* operator->(); //Arrow Operator
};

//The neighbour_iterator class provides an iterator
//over the neighbours of a given vertex. This is
//probably harder (conceptually) than the graph_iterator.
//Unless you know how iterators work.
template <typename vertex>
class neighbour_iterator {

private:

  const directed_graph<vertex> * new_graph; //Pointer to graph address
  vertex new_vertex; //Vertex of whose neighbours we iterate through
  std::size_t current_pos; //Current Position of Iterator

public:
  neighbour_iterator(const neighbour_iterator<vertex>&); //Copy Constructor
  neighbour_iterator(const directed_graph<vertex>&, const vertex&, std::size_t); //Iterator Constructor
  ~neighbour_iterator(); //Deconstructor
  neighbour_iterator<vertex> operator=(const neighbour_iterator<vertex>&); //Assignment Operator
  bool operator==(const neighbour_iterator<vertex>&) const; //Equals Operator
  bool operator!=(const neighbour_iterator<vertex>&) const; //Does not equals Operator
  neighbour_iterator<vertex> operator++(); //Pre-Increment Operator
  neighbour_iterator<vertex> operator++(int); //Post-Increment Operator			
  vertex operator*(); //Indirection Operator
  vertex* operator->(); //Arrow Operator
};


//Define all your methods down here (or move them up into the header, but be careful you don't double up). If you want to move this into another file, you can, but you should #include the file here.
//Although these are just the same names copied from above, you may find a few more clues in the full
//method headers. Note also that C++ is sensitive to the order you declare and define things in - you
//have to have it available before you use it.

template <typename vertex> directed_graph<vertex>::directed_graph() {
	
	this->graph_size = 0; //Initialises Graph size to 0
	this->number_edges = 0; //Initialises number of edges to 0
}

template <typename vertex> directed_graph<vertex>::~directed_graph() {}
	
template <typename vertex> bool directed_graph<vertex>::contains(const vertex& u) const { 
	
	for (unsigned i = 0; i < vertices.size(); i++) { //Loops through the vertices vector O(n)
		
		if (vertices[i] == u) { //If the vertex exists it will return true else it will return false
			
			return true;
			
		}
		
	}
	
	return false; 
}

template <typename vertex> bool directed_graph<vertex>::adjacent(const vertex& u, const vertex& v) { 
	
	return adj_matrix[get_index(u)][get_index(v)]; //Returns either true or false depending if the edge exists
}
	
template <typename vertex> void directed_graph<vertex>::add_vertex(const vertex& u) {
	
	vertices.push_back(u); //Adds the vertex to the vector
	this->graph_size++; //Increases the size of the graph
	
	adj_matrix.resize(this->graph_size); //Resizes the outer vector
	
	for (unsigned i = 0; i < this->graph_size; i++) { // Loops through the vector and increases the size each inner vector in the matrix O(n)

		adj_matrix[i].resize(this->graph_size);
	}
}

template <typename vertex> void directed_graph<vertex>::add_edge(const vertex& u, const vertex& v) {
	
	if ((get_index(u) >= 0) && (get_index(u) < num_vertices()) && (get_index(v) >= 0) && (get_index(v) < num_vertices()) && (get_index(u) != get_index(v))) { // Checks if the input is valid
		
		adj_matrix[get_index(u)][get_index(v)] = true; //Sets the edge to true creating a new edge
		this->number_edges++; 
	}
}

template <typename vertex> void directed_graph<vertex>::remove_vertex(const vertex& u) {
	
	if (get_index(u) != -1) {
	
		int v = get_index(u); //Gets the position of the vertex for use in the below iterators
	
		for (unsigned j = 0; j < this->graph_size; j++) { //Loops through the edges and erases the column of the matrix for the vertex O(n^2)

			adj_matrix[j].erase(adj_matrix[j].begin() + v);
		}
		
		this->number_edges - degree(u); //Removes the edges associated with the vertex
		adj_matrix.erase(adj_matrix.begin() + v); //Removes the row associated with the vertex
		vertices.erase(vertices.begin() + v); //Removes the vertex using iterator
		this->graph_size--;	
	}
}

template <typename vertex> void directed_graph<vertex>::remove_edge(const vertex& u, const vertex& v) {
	
	if ((get_index(u) >= 0) && (get_index(u) < num_vertices()) && (get_index(v) >= 0) && (get_index(v) < num_vertices()) && (get_index(u) != get_index(v))) { // Checks if the input is valid
		
		adj_matrix[get_index(u)][get_index(v)] = false; //Sets the edge to false thereby removing it
		this->number_edges--;
	}
}

template <typename vertex> std::size_t directed_graph<vertex>::in_degree(const vertex& u) { 
	
	int in_count = 0; //Counter used to count number of incoming edges
	
	for (unsigned i = 0; i < num_vertices(); i++) { //Loops through edges coming in O(n)
		
		if (adj_matrix[i][get_index(u)]) { //If incoming edge is found increments count
			
			in_count++;
		}
	}
	
	return in_count; //returns in_degree
}

template <typename vertex> std::size_t directed_graph<vertex>::out_degree(const vertex& u)  { 
	
	int out_count = 0; //Counter used to count number of outgoing edges
	
	for (unsigned i = 0; i < num_vertices(); i++) { //Loops through edges going out O(n)
		
		if (adj_matrix[get_index(u)][i]) { //If outgoing edge is found increments count
			
			out_count++;
		}
	}
	
	return out_count; //Returns the out_degree
}

template <typename vertex> std::size_t directed_graph<vertex>::degree(const vertex& u)  { 
	
	return in_degree(u) + out_degree(u); //Adds the in_degree and out_degree values together
}

template <typename vertex> std::size_t directed_graph<vertex>::num_vertices() const { 
	
	return this->graph_size; //Returns graph_size
}

template <typename vertex> std::size_t directed_graph<vertex>::num_edges() const { 
		
	return number_edges; //Returns the number of edges
}

template <typename vertex> std::vector<vertex> directed_graph<vertex>::get_vertices() const { 
	
	return vertices; //Returns the vector of the vertices
}

template <typename vertex> std::vector<vertex> directed_graph<vertex>::get_neighbours(const vertex& u) const { 
	
	std::vector<vertex> neighbours; //Creates the vector
	
	for (unsigned i = 0; i < num_vertices(); i++) { //Loops through the connecting outgoing vertices O(n)
		
		if (adj_matrix[get_index(u)][i] == true) { //Checks if i is a neighbour e.g. outgoing edge
			
			neighbours.push_back(vertices[i]); //If the vertex is a neighbour adds to the neigbour vector		
		}
	}

	return neighbours;  //Returns the neighbour vector
}

template <typename vertex> vertex_iterator<vertex> directed_graph<vertex>::begin() { return vertex_iterator<vertex>(*this, 0); } //Beings the iterator of this graph at 0

template <typename vertex> vertex_iterator<vertex> directed_graph<vertex>::end() { return vertex_iterator<vertex>(*this, num_vertices()); } //Ends the iterator at the size of the graph

template <typename vertex> neighbour_iterator<vertex> directed_graph<vertex>::nbegin(const vertex& u) { return neighbour_iterator<vertex>(*this, u, 0); } //Begins the iterator at 0

template <typename vertex> neighbour_iterator<vertex> directed_graph<vertex>::nend(const vertex& u) { return neighbour_iterator<vertex>(*this, u, out_degree(u)); } //Ends the iterator at the number of outgoing edges from the vertex

template <typename vertex> std::vector<vertex> directed_graph<vertex>::depth_first(const vertex& u) { 
	
	std::vector<vertex> ordered_list; //The returned ordered traversal
	std::stack<vertex> to_process; //The container of to process vertices, a stack for depth first, hence the order is different to bfs
	std::vector<bool> visited(graph_size); //Keeps track of the visited vertices
	
	for (unsigned i = 0; i < graph_size; i++){ //Sets visited to false
		
		visited[i] = false;	
	}
	
	to_process.push(u); //Starts the traversal
	while(!to_process.empty()) { //Keeps loop going while there are no more vertices to traverse O(n + m)
		
		vertex v = to_process.top(); //The next vertex is added to process
		to_process.pop(); //Removes processed vertex from vector
		
		if (!visited[get_index(v)]) { //Checks if the vertex has been looped and if it has, it is not checked
			
			visited[get_index(v)] = true; //Sets visited
			ordered_list.push_back(v); //Adds to the order if it has not been visited
			
			for (int i = num_vertices(); i != 0; i--) { //Loops through outgoing edges
				
				if (adj_matrix[get_index(v)][i-1]) { //Checks if it has been visited
					
				to_process.push(vertices[i-1]); //Adds the unvisited nodes to be processed
					
				}
			}
		}
	}
	
	return ordered_list; //returns dfs_order
}

template <typename vertex> std::vector<vertex> directed_graph<vertex>::breadth_first(const vertex& u) { 
	
	std::vector<vertex> ordered_list; //The returned ordered traversal
	std::queue<vertex> to_process; //The container of to process vertices, a queue for breadth first, hence the order is different
	std::vector<bool> visited(graph_size); //Keeps track of the visited vertices
	
	for (unsigned i = 0; i < graph_size; i++){ //Sets visited to false
		
		visited[i] = false;	
	}
	
	to_process.push(u); //Starts the traversal
	while(!to_process.empty()) { //Keeps loop going while there are no more vertices to traverse O(n + m)
		
		vertex v = to_process.front(); //The next vertex is added to process
		to_process.pop(); //Removes processed vertex from vector
		
		if (!visited[get_index(v)]) { //Checks if the vertex has been looped and if it has, it is not checked
			
			visited[get_index(v)] = true; //Sets visited
			ordered_list.push_back(v); //Adds to the order if it has not been visited
			
			for (int i = 0; i < graph_size; ++i) { //Loops through outgoing edges
				
				if (adj_matrix[get_index(v)][i]) { //Checks if it has been visited
					
					to_process.push(vertices[i]); //Adds the unvisited nodes to be processed
					
				}
			}
		}
	}
	
	return ordered_list; //returns bfs_order
}

template <typename vertex> directed_graph<vertex> directed_graph<vertex>::out_tree(const vertex& u) { 
	
	std::vector<bool> visited(graph_size); //Keeps a record of the vertices that have been visited
	directed_graph<vertex> o_tree; //Initalises the new tree to be created
	
	for (unsigned i = 0; i < graph_size; i++){ //Sets visited to false
		
		visited[i] = false;	
	}
	
	std::queue<std::pair<vertex, vertex>> to_process; //Creates the edges to add from the two vertices
	
	to_process.push({u, u}); //Starts the process
	
	while (!to_process.empty()) { //While there are still vertices to process the loop continues O(n^2)
		
		std::pair<vertex, vertex> current = to_process.front(); //Initialises the current pair to be tested
		to_process.pop(); //Removes the current from being processed
		
		if (!visited[current.first]) { //Checks if the vertex has been visited
			
			visited[current.first] = true; //Sets true if visited
			o_tree.add_vertex(current.first); //Adds the vertex to the tree
			o_tree.add_edge(current.second, current.first); //Adds the out-going edge to the tree
			
			for (unsigned i = 0; i < graph_size; i++) { //Loops through the vertex edges
				
				if (adj_matrix[get_index(current.first)][i]) { //Tests if visted and if an edge exists
					
					to_process.push({vertices[i], current.first}); //If the above criteia is met, the new vertex is pushed and added to be processed
				}
			}
		}
	}	
	
	return o_tree;
}

template <typename vertex> directed_graph<vertex> directed_graph<vertex>::in_tree(const vertex& u) {
	
	std::vector<bool> visited(graph_size); //Vector that keeps track of which vertex's have been visited
	directed_graph<vertex> i_tree; //Initialises the tree to be created
	
	for (unsigned i = 0; i < graph_size; i++){ //Initialises visited vector as all false
		
		visited[i] = false;
	}
	
	std::queue<std::pair<vertex, vertex>> to_process; //Vector that keeps track of vertex's to process and unordered pairs
	
	to_process.push({u, u}); //Begins the process by pushing pair in
	
	while (!to_process.empty()) { //Loop of process, continues whilst there are still pairs to process O(n^2)
		
		std::pair<vertex, vertex> current = to_process.front(); //Gets pair that will be processed from the top of to_process queue
		to_process.pop(); //Removes the current pair from the two process
		
		if (!visited[current.first]) { //Check if the vertex has been visited
			
			visited[current.first] = true; //Sets true when visited
			i_tree.add_vertex(current.first); //Adds the vertex to the tree
			i_tree.add_edge(current.first, current.second); //Adds the in-coming edge to the tree
			
			for (unsigned i = 0; i < graph_size; i++) { //Begins the loop to check the incoming edges
				
				if (adj_matrix[i][get_index(current.first)]) { //If there is an incoming edge and it hasn't been visited add the pair to be processed
					
					to_process.push({vertices[i], current.first}); //If above conditions are met, add the pair to be processed
				}
			}
		}
	}
	
	return i_tree; //Return Tree
}

template <typename vertex> bool directed_graph<vertex>::reachable(const vertex& u, const vertex& v) {
	
	directed_graph<vertex> v_tree = in_tree(v); //Makes a in_tree based on the vertex v O(n^2)
	
	return v_tree.contains(u); //If the v tree contains u, v can be reached from u as it has in-edges leading to it
}

template <typename vertex> int directed_graph<vertex>::get_index(const vertex& u) const { //This method returns the index for the vertex called
	
	for (unsigned i = 0; i < this->graph_size; i++) { //Loop that goes through the size of the vector and checks O(n) complexity
		
		if (vertices[i] == u) {
			
			return i; //Returns the index if it exists
		}
	}
	
	return -1;
}

template <typename vertex> vertex_iterator<vertex>::vertex_iterator(const vertex_iterator<vertex>& other) { //The copy constructor for the iterator
	
	new_graph = other.new_graph; //Sets the copy to the new memory address if there is one
	current_pos = other.current_pos; //Sets the copy to the new postition
}

template <typename vertex> vertex_iterator<vertex>::vertex_iterator(const directed_graph<vertex>& graph, std::size_t position) : new_graph(&graph), current_pos(position) {} //The constructor for the iterator

template <typename vertex> vertex_iterator<vertex>::~vertex_iterator() {} //Deconstructor

template <typename vertex> vertex_iterator<vertex> vertex_iterator<vertex>::operator=(const vertex_iterator<vertex>& other) { //Equals operator, left side for position assignment

	if (this != &other) { //Checks if the address of the current iterator != the assigned other
		
		current_pos = other.current_pos;
	}
	
	return *this;
}
	

template <typename vertex> bool vertex_iterator<vertex>::operator==(const vertex_iterator<vertex>& other) const {  //Checks if the values are the same
	
	
	return new_graph == other.new_graph && current_pos == other.current_pos; //Checks for equality of objects and position
}

template <typename vertex> bool vertex_iterator<vertex>::operator!=(const vertex_iterator<vertex>& other) const { //Checks if the values are not the same
	
	return other.current_pos != current_pos || other.new_graph != new_graph; //Checks for inequality of objects and position
}

template <typename vertex> vertex_iterator<vertex> vertex_iterator<vertex>::operator++() { //Increments the position of the iterator
	
	++current_pos; //post-increment
	return *this; 
}

template <typename vertex> vertex_iterator<vertex> vertex_iterator<vertex>::operator++(int) { //Increments the positon of the iterator and represents as an int
	
	auto temp = current_pos; //pre-increment initialisation
	++current_pos;
	
	return temp; 
}

template <typename vertex> vertex vertex_iterator<vertex>::operator*() { //Indirection operator, returns what the iterator is pointing at
	
	std::vector<vertex> vertices = new_graph->get_vertices(); //Returns object of vertices vector
	
	return vertices[current_pos]; 
}

template <typename vertex> vertex* vertex_iterator<vertex>::operator->() { //Arrow operator, returns a pointer to the vertex
	
	std::vector<vertex> vertices = new_graph->get_vertices(); //Returns object of vertices vector
	
	return &vertices[current_pos]; 
}

template <typename vertex> neighbour_iterator<vertex>::neighbour_iterator(const neighbour_iterator<vertex>& other) { //Copy constructor for the iterator
	
	new_graph = other.new_graph; //Sets the graph to the copied graph
	current_pos = other.current_pos; //Sets the position to the new position
	new_vertex = other.new_vertex; //Sets the vertex to the new vertex
}

template <typename vertex> neighbour_iterator<vertex>::neighbour_iterator(const directed_graph<vertex>& graph, const vertex& u, std::size_t position) : new_graph(&graph), new_vertex(u), current_pos(position)	{} //Constructor for the iterator

template <typename vertex> neighbour_iterator<vertex>::~neighbour_iterator() {} //Deconstructor

template <typename vertex> neighbour_iterator<vertex> neighbour_iterator<vertex>::operator=(const neighbour_iterator<vertex>& other) { //Equals operator, left side used for setting position 
	
	if (this != &other) { //Checks if the address of the current iterator != the assigned other
		
		current_pos = other.current_pos;
	}
	
	return *this;; 
}

template <typename vertex> bool neighbour_iterator<vertex>::operator==(const neighbour_iterator<vertex>& other) const { //Determins if the iterator is pointing at the same thing
	
	return new_graph == other.new_graph && current_pos == other.current_pos && new_vertex == other.new_vertex; //Checks for equality of objects and position
}

template <typename vertex> bool neighbour_iterator<vertex>::operator!=(const neighbour_iterator<vertex>& other) const { //Determins if the operator is pointing at different things
	
	return other.current_pos != current_pos || other.new_graph != new_graph || other.new_vertex != new_vertex;  //Checks for inequality of objects and positon
}

template <typename vertex> neighbour_iterator<vertex> neighbour_iterator<vertex>::operator++() { //Increments the position of the iterator
	
	++current_pos; //Pre-increment
	return *this; 
}

template <typename vertex> neighbour_iterator<vertex> neighbour_iterator<vertex>::operator++(int) { //Increments the position of the iterator and returns an int
	
	auto temp = current_pos; //Post-increment
	++current_pos;
	
	return temp; 
}		

template <typename vertex> vertex neighbour_iterator<vertex>::operator*() { //Indirection operator, returns what the iterator is pointing at
	
	std::vector<vertex> vertices = new_graph->get_neighbours(new_vertex); //Returns object of neighbours vector
	
	return vertices[current_pos];  
}

template <typename vertex> vertex* neighbour_iterator<vertex>::operator->() { //Arrow operator, returns a pointer to what the iterator is pointing at
	
	std::vector<vertex> vertices = new_graph->get_neighbours(new_vertex); //Returns object of neighbours vector
	
	return &vertices[current_pos];
}


#endif
