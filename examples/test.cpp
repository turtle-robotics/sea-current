#include "../sea_current.hpp"

#include "../third-party/matplotlib-cpp/matplotlibcpp.h"

#include <Eigen/Dense>

#include <toppra/toppra.hpp>
#include <toppra/geometric_path.hpp>
#include <toppra/geometric_path/piecewise_poly_path.hpp>
#include <toppra/constraint/linear_joint_velocity.hpp>
#include <toppra/constraint/linear_joint_acceleration.hpp>
#include <toppra/algorithm/toppra.hpp>
#include <toppra/parametrizer.hpp>
#include <toppra/parametrizer/spline.hpp>

#include <chrono>
#include <iostream>

using namespace Eigen;
using namespace turtle::sc;

namespace plt = matplotlibcpp;

using namespace toppra;
using namespace toppra::constraint;

class TestVLJV : public LinearJointVelocity {
    public:
    TestVLJV(int nDof) : LinearJointVelocity (-1*toppra::Vector::Ones(1), 1*toppra::Vector::Ones(1)) {
        computeVelocityLimits(0);
    }
    protected:
    void computeVelocityLimits(value_type time) {
        value_type slow = 4;
        value_type fast = 6;
        if (time > 0.5) {
            m_lower(0,0) = -slow;
            m_upper(0,0) = slow;
        } else {
            m_lower(0,0) = -fast;
            m_upper(0,0) = fast;
        }
        // std::cout << (m_lower.array() > m_upper.array()).any() << std::endl;
    }
};

int main() {
    dbg_assert(true, "");

    std::vector<Vector2f> ctrl_pts;

    ctrl_pts.push_back(Vector2f(-1.0, 0.0));
    // ctrl_pts.push_back(Vector2f(-0.6, 0.5));
    ctrl_pts.push_back(Vector2f(-0.5, 0.5));
    ctrl_pts.push_back(Vector2f(0.5, -0.5));
    ctrl_pts.push_back(Vector2f(1.0, 0.0));

    auto start = std::chrono::high_resolution_clock::now();
    bezier_spline bs = bezier_spline::bezier_curve(ctrl_pts, 0.0001);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << (duration.count() / 1000.0)<< std::endl;

    std::vector<float> x(bs.n_pts());
    std::vector<float> y(bs.n_pts());
    for (int i = 0; i < bs.n_pts(); ++i) {
        x[i] = bs.pts(i, 0);
        y[i] = bs.pts(i, 1);
    }

    std::vector<float> x2(bs.ctrl_pts[0].size());
    std::vector<float> y2(bs.ctrl_pts[0].size());
    for (int i = 0; i < bs.ctrl_pts[0].size(); ++i) {
        x2[i] = bs.ctrl_pts[0][i].x();
        y2[i] = bs.ctrl_pts[0][i].y();
    }

    bezier_spline deriv = bs.hodograph();

    // float n_arclen = 0;
    // for (int i = 0; i < deriv.n_pts()-1; ++i) {
    //     Vector2f a(deriv.pts(i, 0), deriv.pts(i, 1));
    //     Vector2f b(deriv.pts(i+1, 0), deriv.pts(i+1, 1));
    //     n_arclen += std::sqrt(std::pow(a(1)-b(1), 2) + std::pow(a(0)-b(0), 2));
    // }
    // std::cout << n_arclen << std::endl;

    start = std::chrono::high_resolution_clock::now();
    float arclen = bs.arclength(0.01).arclength;
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    // std::cout << arclen << std::endl;
    std::cout << (duration.count() / 1000.0)<< std::endl;


    VectorXf arcs = bs.arclength().segments[0];
    for (float arc : arcs) {
        // std::cout << arc << std::endl;
    }

    constexpr int degree = 10;
    start = std::chrono::high_resolution_clock::now();

    cheb_poly b = chebfit(bs.pts.col(0), bs.pts.col(1), degree);
    VectorXf y_hat = chebeval(bs.pts.col(0), b, degree);

    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << (duration.count() / 1000.0)<< std::endl;


    std::vector<float> y_hat_p(y_hat.data(), y_hat.data() + y_hat.rows() * y_hat.cols());

    // bsr = bs.resample();

    // plt::subplot(4,1,1);
    // plt::plot(x, y, "");
    // plt::plot(x2, y2, "ro");
    // plt::plot(x, y_hat_p, "");

    start = std::chrono::high_resolution_clock::now();
    TestVLJV test(1);

    toppra::LinearConstraintPtr ljv, lja;
    ljv = std::make_shared<TestVLJV>(test);  //[[-1, 1], [-0.5, 0.5]]
    lja = std::make_shared<toppra::constraint::LinearJointAcceleration>(-20 * toppra::Vector::Ones(1), 20 * toppra::Vector::Ones(1));  //[[-0.05, 0.2], [-0.1, 0.3]]
    lja->discretizationType(toppra::DiscretizationType::Interpolation);
    toppra::LinearConstraintPtrs constraints{ljv, lja};

    constexpr double arcvtest = 10;

    toppra::Vector position0{1}, position1{1};
    position0 << 0.0;
    position1 << arclen;

    toppra::Vectors positions = {position0, position1};

    toppra::Vector velocity0{1}, velocity1{1};
    velocity0 << 0.0;
    velocity1 << 0.0;

    toppra::Vectors velocities = {velocity0, velocity1};

    std::vector<toppra::value_type> steps;
    steps = std::vector<toppra::value_type>{0, 1};

    toppra::PiecewisePolyPath hermite = toppra::PiecewisePolyPath::CubicHermiteSpline(positions, velocities, steps);


    // toppra::Matrix arcwrap(1,1);
    // arcwrap(0,0) = arcvtest;
    // toppra::Matrices coefficients;
    // // coefficients.push_back(toppra::Matrix::Zero(1,1));
    // // coefficients.push_back(toppra::Matrix::Zero(1,1));
    // coefficients.push_back(arcwrap);
    // coefficients.push_back(toppra::Matrix::Zero(1,1));

    // std::cout << coefficients.size() << std::endl;

    // std::vector<toppra::value_type> breakpoints = {0, arcvtest/2, arcvtest};

    // PiecewisePolyPath path(coefficients, breakpoints);

    // toppra::LinearConstraintPtrs constraints{test, lja};
    toppra::GeometricPathPtr path = std::make_shared<PiecewisePolyPath>(hermite);

    toppra::algorithm::TOPPRA algo(constraints, path);
    toppra::ReturnCode rc1 = algo.computePathParametrization(0, 0);

    dbg_assert(rc1 == toppra::ReturnCode::OK, "");

    toppra::ParametrizationData pd = algo.getParameterizationData();

    toppra::Vector gridpoints = pd.gridpoints;  // Grid-points used for solving the discretized problem.
    toppra::Vector vsquared = pd.parametrization;  // Output parametrization (squared path velocity)
    toppra::parametrizer::Spline spp(path, gridpoints, vsquared);

    Eigen::Matrix<toppra::value_type, 1, 2> interval2;
    interval2 = spp.pathInterval();

    const int length2 = 100;
    toppra::Vector times2 = toppra::Vector::LinSpaced(length2, interval2(0), interval2(1));

    toppra::Vectors path_pos2;
    path_pos2 = spp.eval(times2, 0);  // TODO this function call fails
    toppra::Vectors path_vel2;
    path_vel2 = spp.eval(times2, 1);  // TODO this function call fails
    toppra::Vectors path_acc2;
    path_acc2 = spp.eval(times2, 2);  // TODO this function call fails

    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << (duration.count() / 1000.0)<< std::endl;

    std::cout << path_pos2.size() << std::endl;

    VectorXf pos_plot(length2);
    VectorXf vel_plot(length2);
    VectorXf acc_plot(length2);
    for (int i = 0; i < length2; ++i) {
        for (auto val : path_pos2[i]) {
            pos_plot(i) = val;
            // std::cout << val << std::endl;
        }

        for (auto val : path_vel2[i]) {
            vel_plot(i) = val;
            // std::cout << val << std::endl;
        }

        for (auto val : path_acc2[i]) {
            acc_plot(i) = val;
        }
    }
    // pos_plot(0) = 0;
    arclength_data ad = bs.arclength();
    bezier_spline re = bs.resample(pos_plot, ad, true);

    std::vector<float> x3(re.n_pts());
    std::vector<float> y3(re.n_pts());
    for (int i = 0; i < re.n_pts(); ++i) {
        x3[i] = re.pts(i, 0);
        y3[i] = re.pts(i, 1);
        // std::cout << x3[i] << " " << y3[i] << std::endl;
    }

    plt::plot(x3, y3, "");


    // // plt::subplot(4,1,2);
    // plt::plot(std::vector<double>(times2.data(), times2.data()+times2.size()), std::vector<double>(pos_plot.data(), pos_plot.data()+pos_plot.size()), "");

    // // plt::subplot(4,1,3);
    // plt::plot(std::vector<double>(times2.data(), times2.data()+times2.size()), std::vector<double>(vel_plot.data(), vel_plot.data()+vel_plot.size()), "");

    // // plt::subplot(4,1,4);
    // plt::plot(std::vector<double>(times2.data(), times2.data()+times2.size()), std::vector<double>(acc_plot.data(), acc_plot.data()+acc_plot.size()), "");


    plt::show();
}
