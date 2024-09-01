#ifndef RIGIDBODY_NORMAL_H
#define RIGIDBODY_NORMAL_H

#include <tuple>

#include "Eigen/Dense"

using Normal=Eigen::Vector2d;
Normal N_0=Normal::Zero();

std::tuple<double, double> SinCos(const Normal& n) {
    double theta = std::atan2(n.y(), n.x()); // Compute angle theta
    return {std::sin(theta), std::cos(theta)}; // Return sin and cos of theta
}


#endif//RIGIDBODY_NORMAL_H
