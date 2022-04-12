#pragma once

#include "cloud.hpp"
#include "radial.hpp"

// TODO
class gaussian :public radial {

public:

	double volume();

	double profile(double);

	gaussian(cloud*, double);

	void guess_bandwidth();
};