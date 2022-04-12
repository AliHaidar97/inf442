#include <iostream>
#include <fstream>
#include<string>
#include "point.cpp"
#include "cloud.cpp"
#include<queue>
#include<utility>

int main (int argc, char **argv)
{
	
	int dimension = std::stoi(argv[1]);
	std::ifstream is(argv[2]);
	int k = std::stoi(argv[3]);
	int iterations = std::stoi(argv[4]);
	cloud data(dimension, 100000);
	data.load(is);
	for (int i = 0; i < iterations; i++) {
		cloud temp(dimension, 100000);
		for (int j = 0; j < data.get_n(); j++){
			std::priority_queue < std :: pair<double, int >> s;
			for (int v = 0; v < data.get_n(); v++) {
				//if (v == j)continue;
				s.push(std::make_pair(data.get_point(j).dist(data.get_point(v)), v));
				if (s.size() > k) {
					s.pop();
				}
			}
			point p;
			while (s.size()) {
				for (int d = 0; d < dimension; d++) {
					p.coords[d] += data.get_point(s.top().second).coords[d];
				}
				s.pop();
			}

			for (int d = 0; d < dimension; d++) {
				p.coords[d] /= k;
			}
			temp.add_point(p);
		}
		for (int j = 0; j < data.get_n(); j++) {
			for (int d = 0; d < dimension; d++) {
				data.get_point(j).coords[d]=temp.get_point(j).coords[d];
			}
		}
	}
	freopen("out.txt", "w", stdout);
	for (int i = 0; i < data.get_n(); i++) {
		for (int d = 0; d < dimension; d++) {
			std::cout << data.get_point(i).coords[d] << ' ';
		}std::cout << '\n';
	}
	fclose(stdout);
	return 0;
}
