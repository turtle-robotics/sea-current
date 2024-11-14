#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <numbers>
#include <algorithm>
#include <tuple>
#include <functional>
#include <queue>
#include <limits>
#include <unordered_set>
#include <unordered_map>
#include <optional>
#include <string>

#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

#include <toppra/toppra.hpp>
#include <toppra/geometric_path.hpp>
#include <toppra/geometric_path/piecewise_poly_path.hpp>
#include <toppra/constraint/linear_joint_velocity.hpp>
#include <toppra/constraint/linear_joint_acceleration.hpp>
#include <toppra/algorithm/toppra.hpp>
#include <toppra/parametrizer.hpp>
#include <toppra/parametrizer/spline.hpp>

#include "utilities.hpp"
#include "planning_space.hpp"
#include "planner.hpp"
#include "bezier.hpp"

using toppra::value_type;
using toppra::constraint::LinearJointVelocity;
using vel_lim_func = std::function<std::tuple<toppra::Vector, toppra::Vector>(value_type time)>;

class LinearJointVelocityVarying : public LinearJointVelocity {
    public:
        vel_lim_func calc_lim;
        LinearJointVelocityVarying(int nDof, vel_lim_func calc_lim);

    protected:
        void computeVelocityLimits(value_type time);
};


LinearJointVelocityVarying::LinearJointVelocityVarying(int nDof, vel_lim_func calc_lim) : 
LinearJointVelocity (-1*toppra::Vector::Ones(1), 1*toppra::Vector::Ones(1)) {
    this->calc_lim = calc_lim;
    computeVelocityLimits(0);
}

void LinearJointVelocityVarying::computeVelocityLimits(value_type time) {
    std::tie(m_lower, m_upper) = calc_lim(time);
}