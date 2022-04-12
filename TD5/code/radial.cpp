#include <cmath>

#include "point.hpp"
#include "cloud.hpp"
#include "radial.hpp"
#include<iostream>
// TODO
radial::radial(cloud* data_, double bandwidth_):kernel(data_) {

	bandwidth = bandwidth_;
}

double radial::density(point& p) {

	double coef = (volume() * pow(bandwidth, point::get_dim()) * data->get_n());
	double value = 0;

	for (int i = 0; i < data->get_n(); i++) {
		value += profile(pow(p.dist(data->get_point(i))/ bandwidth, 2));
	}

	return  value/coef;
}


