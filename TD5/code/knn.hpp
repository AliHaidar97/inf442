#pragma once

#include "point.hpp"
#include "cloud.hpp"
#include "kernel.hpp"

// TODO
class knn :public kernel {
private:
	int k,V;
public:
	knn(cloud*,int,double);
	double density(point&);


};