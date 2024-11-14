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
#include "LinearJointVelocityVarying.hpp"

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#ifdef DEBUG
#define SC_ASSERT(cnd, msg)                                                           \
    do {                                                                              \
        static_assert(                                                                \
                !std::is_pointer_v<decltype(cnd)>,                                    \
                "Do not use SC_ASSERT with raw pointers"                              \
                "and instead do SC_ASSERT(cnd != nullptr) "                           \
                "to avoid implicit pointer-to-bool conversion.");                     \
        if (bool(cnd) == false) {                                                     \
            std::cerr << "SC_ASSERT failed: " << #cnd << std::endl;                   \
            std::cerr << "SC_ASSERT message: " << msg << std::endl;                   \
            std::cerr << "SC_ASSERT failed at: ";                                     \
            std::cerr << __func__ << " " << __FILE__ << " " << __LINE__ << std::endl; \
            std::exit(1);                                                             \
        }                                                                             \
    } while(0)
#else
#define SC_ASSERT(cnd, msg)
#endif


namespace turtle::sc {

    using namespace Eigen;
    using namespace std::complex_literals;

    using point_set = std::unordered_set<Vector2f, hash_vector2f>;

    using toppra::value_type;
    using toppra::constraint::LinearJointVelocity;

    using vel_lim_func = std::function<std::tuple<toppra::Vector, toppra::Vector>(value_type time)>;

    template <int N> requires (N >= 1)
    velocity_profile gen_vel_prof(const Vector<value_type, N>& pos_end,
                                  const Vector<value_type, N>& pos_start,
                                  const Vector<value_type, N>& vel_end,
                                  const Vector<value_type, N>& vel_start,
                                  const vel_lim_func& vel_lim,
                                  const Vector<value_type, N>& acc_min,
                                  const Vector<value_type, N>& acc_max,
                                  const float dt=0.02) {
        using namespace toppra;
        using namespace toppra::constraint;

        const int dof = pos_end.rows();

        LinearJointVelocityVarying vel_con(dof, vel_lim);

        toppra::LinearConstraintPtr ljv, lja;
        ljv = std::make_shared<LinearJointVelocityVarying>(vel_con);
        lja = std::make_shared<toppra::constraint::LinearJointAcceleration>(acc_min, acc_max);
        lja->discretizationType(toppra::DiscretizationType::Interpolation);
        toppra::LinearConstraintPtrs constraints{ljv, lja};

        toppra::Vectors positions = {pos_start, pos_end};

        toppra::Vectors velocities = {vel_start, vel_end};

        std::vector<toppra::value_type> steps;
        steps = std::vector<toppra::value_type>{0, 1};

        toppra::PiecewisePolyPath hermite = toppra::PiecewisePolyPath::CubicHermiteSpline(positions, velocities, steps);

        toppra::GeometricPathPtr path = std::make_shared<PiecewisePolyPath>(hermite);

        toppra::algorithm::TOPPRA algo(constraints, path);
        toppra::ReturnCode rc1 = algo.computePathParametrization(0, 0);

        SC_ASSERT(rc1 == toppra::ReturnCode::OK, "");

        toppra::ParametrizationData pd = algo.getParameterizationData();

        toppra::Vector gridpoints = pd.gridpoints;
        toppra::Vector vsquared = pd.parametrization;
        toppra::parametrizer::Spline spp(path, gridpoints, vsquared);

        Eigen::Matrix<toppra::value_type, 1, 2> interval = spp.pathInterval();

        int length = std::ceil((interval(1) - interval(0)) / dt);
        // std::cout << "length " << length << std::endl;
        toppra::Vector times = toppra::Vector::LinSpaced(length, interval(0), interval(1));

        toppra::Vectors path_pos = spp.eval(times, 0);
        toppra::Vectors path_vel = spp.eval(times, 1);
        toppra::Vectors path_acc = spp.eval(times, 2);

        std::vector<VectorXf> pos(dof);
        std::vector<VectorXf> vel(dof);
        std::vector<VectorXf> acc(dof);

        for (int j = 0; j < dof; ++j) {
            pos[j] = VectorXf::Zero(length);
            vel[j] = VectorXf::Zero(length);
            acc[j] = VectorXf::Zero(length);
        }

        // TODO: there is probably a better way to copy this data
        for (int i = 0; i < path_pos.size(); ++i) {
            for (int j = 0; j < dof; ++j) {
                pos[j](i) = path_pos[i](j);
                vel[j](i) = path_vel[i](j);
                acc[j](i) = path_acc[i](j);
            }
        }

        return velocity_profile(pos, vel, acc, times);
    }

    template <typename T>
    std::vector<T> format_vec_vecx(const std::vector<VectorXf>& prof) {
        std::vector<T> prof_ser;
        for (std::size_t i = 0; i < prof.size(); ++i) {
            std::vector<T> tmp = std::vector<T>(prof[i].data(), prof[i].data()+prof[i].size());
            for (std::size_t j = 0; j < tmp.size(); ++j) {
                prof_ser.push_back(tmp[i]);
            }
        }
        return prof_ser;
    }

    // quick and dirty serialization
    json serialize_path_to_json(const bezier_spline& spline, const velocity_profile& vel_prof, const arclength_data& arclens, const std::vector<float>& ang_vel) {
        json j = json::array();

        //j["position"] = format_vec_vecx<float>(vel_prof.pos);
        const std::vector<float> velocities = format_vec_vecx<float>(vel_prof.vel);
        const std::vector<float> accelerations = format_vec_vecx<float>(vel_prof.acc);
        const std::vector<float> times = std::vector<float>(vel_prof.time.data(), vel_prof.time.data() + vel_prof.time.size());

        auto&& pts = spline.pts;
        std::vector<float> pos_x = std::vector<float>(pts.col(0).data(), pts.col(0).data() + pts.rows());
        std::vector<float> pos_y = std::vector<float>(pts.col(1).data(), pts.col(1).data() + pts.rows());

        // const std::vector<std::vector<float>> segments = format_vec_vecx<float>(arclens.segments);
        // const std::vector<std::vector<float>> positions = format_vec_vecx<float>(arclens.positions);
        // j["arclength"] = { {"arclength", arclens.arclength}, {"segments", segments}, {"positions", positions} };

        for (std::size_t i = 0; i < pts.rows(); ++i) {
            json j2;
            j2["time"] = times[i];
            j2["velocity"] = velocities[i];
            j2["acceleration"] = accelerations[i];
            j2["angularVelocity"] = ang_vel[i];
            // j2["curvature"]

            // j2["pose"]["rotation"]["radians"] 
            j2["pose"]["translation"]["x"] = pos_x[i];
            j2["pose"]["translation"]["y"] = pos_y[i];

            j2["holonomicRotation"] = 0.0;
            j2["holonomicAngularVelocity"] = 0.0;
            j.push_back(j2);
        }

        return j;
    }
}
