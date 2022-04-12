#include <cmath>

#include "point.hpp"
#include "flat.hpp"

// TODO
double flat::volume() {
	double pie = acos(-1);
	double gamma = tgamma(point::get_dim() / 2.0 + 1);
	return pow(pie, point::get_dim() / 2.0) / gamma;
}

double flat::profile(double t) {
	return 1 ? (t <= 1) : 0;
}

flat::flat(cloud* data_, double bandwidth_) :radial(data_, bandwidth_) {}

