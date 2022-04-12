#pragma once

#include "point.hpp"
#include "cloud.hpp"
#include "kernel.hpp"

// TODO
class radial:public kernel {

protected:
	double bandwidth;

public:

	virtual double volume() = 0;

	virtual double profile(double t) = 0;

	radial(cloud*, double);
		
	double density(point&);

};