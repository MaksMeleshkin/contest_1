#include <iostream>
#include <vector>
#include <queue>
#include <utility>
#include <algorithm>
#include <cassert>
#include <set>

class Graph {
public:
	typedef size_t Vertex;
	Graph(size_t count, bool is_directed) {
		is_dir = is_directed;
		vertex_count = count;
	}

	virtual std::vector<Vertex> get_neighbors(Vertex v) const = 0;

	size_t get_vertex_count() const {
		return vertex_count;
	}

	bool is_directed() const {
		return is_dir;
	}

	virtual void add_edge(const Vertex &start, const Vertex &finish) = 0;

	size_t get_edge_count() const {
		return edge_count;
	}

protected:
	bool is_dir;
	size_t vertex_count, edge_count = 0;
};


class GraphAdjList: public Graph {
private:
	std::vector<std::vector<Vertex>> adj_list;
public:
	GraphAdjList(size_t count, bool is_directed): Graph(count, is_directed), adj_list(count, std::vector<Vertex>()) {}
	
	virtual void add_edge(const Vertex &start, const Vertex &finish) {
		++edge_count;
		adj_list[start].push_back(finish);
		if (!is_dir) {
			adj_list[finish].push_back(start);
		}
	}
	
	virtual std::vector<Vertex> get_neighbors(Vertex v) const {
		return adj_list[v];
	}
};

class GraphAdjMatrix: public Graph {
private:
	std::vector<std::vector<bool>> adj_matrix;
public:
	GraphAdjMatrix (size_t count, bool is_directed): Graph(count, is_directed), adj_matrix(count, std::vector<bool>(count, false)) {}

	virtual void add_edge(const Vertex &start, const Vertex &finish) {
		++edge_count;
		adj_matrix[start][finish] = true;
		if (!is_dir) {
			adj_matrix[finish][start] = true;
		}
	}

	virtual std::vector<Vertex> get_neighbors(Vertex v) const {
		std::vector<Vertex> result;
		for (Vertex i = 0; i < vertex_count; ++i) {
			if (adj_matrix[v][i]) {
				result.push_back(i);
			}
		}
		return result;
	}
};
////////////////////////////////////////////////////////////////////////////
class WeightedGraph {
typedef size_t Vertex;
private:
	std::vector<std::vector<std::pair<Vertex, int>>> adj_list;
	size_t vertex_count, edge_count = 0;
	bool is_dir;
public:
	WeightedGraph(size_t count, bool is_directed): is_dir(is_directed), vertex_count(count), adj_list(count, std::vector<std::pair<Vertex, int>>()) {} 

	void add_edge(const Vertex &from, const Vertex &to, const int &weight) {
		++edge_count;
		adj_list[from].push_back(std::pair<Vertex, int>(to, weight));
		if (!is_dir) {
			adj_list[to].push_back(std::pair<Vertex, int>(from, weight));
		}
	} 

	std::vector<std::pair<Vertex, int>> get_neighbors(const Vertex v) const{
		return adj_list[v];
	}

	bool is_directed() const {
		return is_dir;
	}

	size_t get_vertex_count() const {
		return vertex_count;
	}
	
	size_t get_edge_count() const {
		return edge_count;
	}
};

class Edge {
private:
	size_t from, to;
	size_t weight;
public:
	Edge(size_t start, size_t finish, size_t w): from(start), to(finish), weight(w) {}
	size_t get_from() const {return from;}
	size_t get_to() const {return to;}
	size_t get_weight() const {return weight;}
};
////////////////////////////////////////////////////////////////////////////
namespace GraphProcessing {
	namespace{
		typedef size_t Vertex;
		enum VertexMark {
			WHITE,
			GRAY,
			BLACK
		};
		
		const int NO_PARENT = -1;
		const int INF = 1000000000;

		void dfs_for_art_points(std::vector<VertexMark> &vertex_mark, const Graph &g, const Vertex v, int parent, std::vector<size_t> &low, std::vector<size_t> &enter, size_t &time, std::set<Vertex> &points_list) {
			vertex_mark[v] = GRAY;
			++time;
			low[v] = time;
			enter[v] = time;
			int children = 0;
			for (Vertex i: g.get_neighbors(v)) {
				if (parent == i) {
					continue;
				}
				if (vertex_mark[i] == GRAY) {
					low[v] = std::min(low[v], enter[i]);
				} else {
					dfs_for_art_points(vertex_mark, g, i, v, low, enter, time, points_list);
					low[v] = std::min(low[v], low[i]);
					if (low[i] >= enter[v] && parent != NO_PARENT) {
						points_list.insert(v);
					} 
					++children;
				}
			}
			if (children > 1 && parent == NO_PARENT) {
				points_list.insert(v);
			}
		}
		
		void dfs_for_cc(std::vector<VertexMark> &vertex_mark, const Graph &g, const Vertex v, std::vector<std::vector<Vertex>> &components) {
			vertex_mark[v] = GRAY;
			(*components.rbegin()).push_back(v);
			for (Vertex i: g.get_neighbors(v)) {
				if (vertex_mark[i] == WHITE) {
					dfs_for_cc(vertex_mark, g, i, components);
				}
			}
			vertex_mark[v] = BLACK;
		}

		void dfs_for_ts(std::vector<VertexMark> &vertex_mark, const Graph &g, const Vertex v, std::vector<Vertex> &order) {
			vertex_mark[v] = GRAY;
			for (Vertex i: g.get_neighbors(v)) {
				if (vertex_mark[i] == WHITE) {
					dfs_for_ts(vertex_mark, g, i, order); 
				}
			}
			order.push_back(v);
		}

		void dfs_for_bridges(std::vector<VertexMark> &vertex_mark, const WeightedGraph &g, const Vertex v, size_t &time, std::vector<size_t> &low, std::vector<size_t> &enter, std::vector<Edge> &bridges, const Vertex parent) {
			++time;
			enter[v] = time;
			low[v] = time;
			vertex_mark[v] = GRAY;
			for (std::pair<Vertex, size_t> i: g.get_neighbors(v)) {
				if (i.first == parent) {
					continue;
				}
				if (vertex_mark[i.first] == GRAY) {
					low[v] = std::min(low[v], enter[i.first]);
				} else {
					dfs_for_bridges(vertex_mark, g, i.first, time, low, enter, bridges, v);
					low[v] = std::min(low[v], low[i.first]);
					if (low[i.first] > enter[v]) {
						int count_of_mul_edges = 0;
						for (auto j: g.get_neighbors(v)) {
							if (i.first == j.first) {
								++count_of_mul_edges;
							}
						}
						if (count_of_mul_edges == 1) { 
							bridges.push_back(Edge(v, i.first, i.second));
						}
					}
				}	
			}
		}

		bool dfs(std::vector<VertexMark> &vertex_mark, const Graph &g, const Vertex v) {
			vertex_mark[v] = GRAY;
			for (Vertex i: g.get_neighbors(v)) {
				if (vertex_mark[i] == WHITE) {
					if (!dfs(vertex_mark, g, i)) {
						return 0;
					} 
				}
				if (vertex_mark[i] == GRAY) {
					return 0;
				}
			}
			vertex_mark[v] = BLACK;
			return 1;
		}

		void dfs_for_scc(std::vector<VertexMark> &vertex_mark, const Graph &g, const Vertex v, std::vector<Vertex> &scc, const size_t &scc_id) {
			vertex_mark[v] = GRAY;
			scc[v] = scc_id;
			for (Vertex i: g.get_neighbors(v)) {
				if (vertex_mark[i] == WHITE) {
					dfs_for_scc(vertex_mark, g, i, scc, scc_id);
				}
			}
		}
		
		class cmp_for_prim{
		public:
			bool operator() (std::pair<Vertex, int> a, std::pair<Vertex, int> b) {
				return a.second > b.second;
			}
		};
	
	}
	std::vector<Vertex> top_sort(const Graph &g) {
		assert(g.is_directed());
		std::vector<VertexMark> vertex_mark(g.get_vertex_count(), WHITE);
		std::vector<Vertex> ts_list;
		for (Vertex i = 0; i < g.get_vertex_count(); ++i) {
			if (vertex_mark[i] == WHITE) {
				dfs_for_ts(vertex_mark, g, i, ts_list);
			}
		}
		return ts_list;
	}

	void bfs(const Graph &g, Vertex v, std::vector<size_t> &distance, std::vector<int> &parent) {
		std::vector<VertexMark> vertex_mark(g.get_vertex_count(), WHITE);
		std::queue<Vertex> vertex_queue;
		vertex_mark[v] = GRAY;
		vertex_queue.push(v);
		while (!vertex_queue.empty()) {
			Vertex u = vertex_queue.front();
			vertex_queue.pop();
			for (Vertex i: g.get_neighbors(u)) {
				if (vertex_mark[i] == WHITE) {
					vertex_mark[i] = GRAY;
					distance[i] = distance[u] + 1;
					parent[i] = u;
					vertex_queue.push(i);
				}
			}
			vertex_mark[u] = BLACK;
		}
	}
	
	std::vector<Vertex> shortest_path(const Graph &g, const Vertex start, const Vertex finish) {
		std::vector<size_t> distance(g.get_vertex_count(), 0);
		std::vector<int> parent(g.get_vertex_count(), NO_PARENT);
		std::vector<Vertex> result;
		bfs(g, start, distance, parent);
		if (start == finish) {
			result.push_back(start);
			return result;
		}
		if (parent[finish] == NO_PARENT) {
			return result;
		}
		size_t current = finish;
		while (start != current) {
			result.push_back(current);
			current = parent[current];
		}
		result.push_back(start);
		return result;
	}	

	std::vector<std::vector<Vertex>> getConnectedComponents(const Graph &g) {
		std::vector<std::vector<Vertex>> components;
		std::vector<VertexMark> vertex_mark(g.get_vertex_count(), WHITE);
		for (Vertex v = 0; v < g.get_vertex_count(); ++v) {
			if (vertex_mark[v] == WHITE) {
				components.push_back(std::vector<Vertex>());
				dfs_for_cc(vertex_mark, g, v, components);
			}
		}
		return components;
	}

	bool check_acyclicity(const Graph &g) {
		if (!g.is_directed()) {
			return false;
		}
		std::vector<VertexMark> vertex_mark(g.get_vertex_count(), WHITE);
		for (Vertex i = 0; i < g.get_vertex_count(); ++i) {
			if (vertex_mark[i] == WHITE) {
				if (!dfs(vertex_mark, g, i)) {
					return 0;
				}
			}
		}
		return 1;
	}
	
	std::vector<Edge> find_bridges(WeightedGraph &g) {
		std::vector<VertexMark> vertex_mark(g.get_vertex_count(), WHITE);
		std::vector<Vertex> low(g.get_vertex_count(), g.get_vertex_count());
		std::vector<Vertex> enter(g.get_vertex_count(), g.get_vertex_count());
		std::vector<Edge> bridges;
		size_t time = 0;
		for (Vertex i = 0; i < g.get_vertex_count(); ++i) {
			if (vertex_mark[i] == WHITE) {
				dfs_for_bridges(vertex_mark, g, i, time, low, enter, bridges, 0);
			}
		}
		return bridges;
	}
	
	std::pair<std::vector<Vertex>, size_t> SCC(const Graph &g) {
		std::vector<VertexMark> vertex_mark(g.get_vertex_count(), WHITE);
		std::vector<Vertex> order = top_sort(g);
		size_t scc_id = 0;
		GraphAdjList g_transposed(g.get_vertex_count(), true);
		for (Vertex i = 0; i < g.get_vertex_count(); ++i) {
			for (Vertex j: g.get_neighbors(i)) {
				g_transposed.add_edge(j, i);
			}
		}
		for (auto i: vertex_mark) {
			i = WHITE;
		}
		std::vector<Vertex> scc(g.get_vertex_count(), 0);
		for (auto i = order.rbegin(); i != order.rend(); ++i) {
			if (vertex_mark[*i] == WHITE) {
				++scc_id;
				dfs_for_scc(vertex_mark, g_transposed, *i, scc, scc_id);
			}
		}
		return std::pair<std::vector<Vertex>, size_t>(scc, scc_id); 
	}

	std::set<Vertex> find_art_points(const Graph &g) {
		std::vector<VertexMark> vertex_mark(g.get_vertex_count(), WHITE);	
		std::vector<Vertex> low(g.get_vertex_count(), g.get_vertex_count());
		std::vector<Vertex> enter(g.get_vertex_count(), g.get_vertex_count());
		std::set<Vertex> art_points;
		size_t time = 0;
		for (Vertex i = 0; i < g.get_vertex_count(); ++i) {
			if (vertex_mark[i] == WHITE) {
				dfs_for_art_points(vertex_mark, g, i, NO_PARENT, low, enter, time, art_points);
			}
		}
		return art_points;
	}

	int find_MST_weight(WeightedGraph &g) {
		assert(!g.is_directed());
		std::vector<bool> used(g.get_vertex_count(), false);
		int mst_weight = 0;
		std::priority_queue<std::pair<Vertex, int>, std::vector<std::pair<Vertex, int>>, cmp_for_prim> q;
		q.push(std::pair<Vertex, int>(0, 0));	
		while (!q.empty()) {
			Vertex v = q.top().first;
			int d = q.top().second;
			q.pop();
			if (used[v]) {
				continue;
			}	
			used[v]	= true;	
			mst_weight += d;
			for (std::pair<Vertex, int> i: g.get_neighbors(v)) {
				if (!used[i.first]) {
					q.push(i);
				}
			}
		}	
		return mst_weight;
	}	
	
}

int main() {
	int n, m;
	std::cin >> n >> m;
	WeightedGraph g(n, false);
	for (int i = 0; i < m; ++i) {
		size_t from, to;
		int weight;
		std::cin >> from >> to >> weight;
		g.add_edge(from - 1, to - 1, weight);
	}
	std::cout << GraphProcessing::find_MST_weight(g);
}
