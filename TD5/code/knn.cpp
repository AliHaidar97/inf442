#include "point.hpp"
#include "cloud.hpp"
#include "kernel.hpp"
#include "knn.hpp"

knn::knn(cloud* data_, int k_, double v_) :kernel(data_) {
	k = k_;
	V = v_;
}

double knn::density(point& p) {
	return k / (2 * data->get_n() * V * data->knn(p,k));
}