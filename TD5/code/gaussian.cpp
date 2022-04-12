#include <cmath>
#include <iostream>

#include "point.hpp"
#include "cloud.hpp"
#include "gaussian.hpp"

// TODO
double gaussian::volume() {
	double pie = acos(-1);
	return pow(2*pie, point::get_dim() / 2.0);
}

double gaussian::profile(double t) {
	return exp(-t / 2.0);
}

void gaussian::guess_bandwidth() {
	bandwidth = 1.06 * data->standard_deviation() / pow(data->get_n(), 0.2);
}

gaussian::gaussian(cloud* data_, double bandwidth_) :radial(data_, bandwidth_) {}
