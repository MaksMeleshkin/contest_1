#include <vector>
#include <queue>
#include <iostream>

std::vector<int> foo(const int n) {
	std::vector<int> res;
	if (n / 1000 != 9) {
		res.push_back(n + 1000);
	} else {
		res.push_back(n);
	}
	if (n % 10 != 1) {
		res.push_back(n - 1);
	} else {
		res.push_back(n);
	}
	res.push_back((n % 1000) * 10 + n / 1000);
	res.push_back((n % 10) * 1000 + n / 10);
	return res;
}

std::vector<int> bfs(int start, int finish) {
	const int NO_PARENT = -1;
	const int MAX_N = 9999;
	std::vector<bool> used(MAX_N, false);	
	std::queue<int> q;
	std::vector<int> parent(MAX_N, NO_PARENT);
	used[start] = true;
	q.push(start);
	while (!q.empty()) {
		int v = q.front();
		q.pop();
		if (v == finish) break;
		for (int i: foo(v)) {	
			if (!used[i]) {
				used[i] = true;
				q.push(i);
				parent[i] = v;
			}
		}
	}
	std::vector<int> path;
	int current = parent[finish];
	path.push_back(finish);
	while (current != start) {
		path.push_back(current);
		current = parent[current];
	}
	path.push_back(start);
	return path;
}

int main() {
	int start, finish;
	std::cin >> start >> finish;
	std::vector<int> path = bfs(start, finish);
	for (auto i = path.rbegin(); i != path.rend(); ++i) {
		std::cout << *i << std::endl;
	}		
}		
