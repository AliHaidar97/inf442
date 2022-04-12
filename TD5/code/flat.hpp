#pragma once

#include "cloud.hpp"
#include "radial.hpp"

// TODO
class flat :public radial {
	
public:

	double volume();

	double profile(double );

	flat(cloud*, double);
};