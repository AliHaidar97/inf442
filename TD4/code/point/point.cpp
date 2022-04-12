#include "point.hpp"
#include<iostream>



int point::d = 0;

point::point() {

	coords = new double[d];
	label = 0;

	for (int i = 0; i < d; i++) {
		coords[i] = 0;
	}
}


point::~point() {
	delete[]coords;
}


bool point::set_dim(int _d) {
	if (d > 0) {
		return false;
	}
	d = _d;
	return true;
}

int point::get_dim() {
	return d;
}

void point::print() {

	for (int i = 0; i < d; i++) {
		if (i > 0) {
			std::cout << '\t';
		}
		std::cout << coords[i];
	}
	std::cout << '\n';
}

double point::squared_dist(point& q) {

	double dis = 0;

	for (int i = 0; i < d; i++) {
		dis += (coords[i] - q.coords[i]) * (coords[i] - q.coords[i]);
	}

	return dis;

}