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
    std::cout << arclen << std::endl;
    std::cout << (duration.count() / 1000.0)<< std::endl;


    VectorXf arcs = bs.arclength().segments[0];
    for (float arc : arcs) {
        // std::cout << arc << std::endl;
    }

    constexpr int degree = 10;
    start = std::chrono::high_resolution_clock::now();

    chebpoly b = chebfit(bs.pts.col(0), bs.pts.col(1), degree);
    VectorXf y_hat = chebeval(bs.pts.col(0), b, degree);

    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << (duration.count() / 1000.0)<< std::endl;


    std::vector<float> y_hat_p(y_hat.data(), y_hat.data() + y_hat.rows() * y_hat.cols());


    // plt::subplot(4,1,1);
    // plt::plot(x, y, "");
    // plt::plot(x2, y2, "ro");
    // plt::plot(x, y_hat_p, "");


    Eigen::Vector<value_type, 1> pos_end{arclen};
    Eigen::Vector<value_type, 1> pos_start{0};

    Eigen::Vector<value_type, 1> vel_end{0};
    Eigen::Vector<value_type, 1> vel_start{0};

    Eigen::Vector<value_type, 1> acc_min{-40};
    Eigen::Vector<value_type, 1> acc_max{40};

    auto vel_lim = [](toppra::value_type time) {
        toppra::Vector lower{1};
        toppra::Vector upper{1};

        value_type slow = 4;
        value_type fast = 6;
        if (time > 0.5) {
            lower(0,0) = -slow;
            upper(0,0) = slow;
        } else {
            lower(0,0) = -fast;
            upper(0,0) = fast;
        }

        return std::make_tuple(lower, upper);
    };

    start = std::chrono::high_resolution_clock::now();

    velocity_profile prof = gen_vel_prof<1>(pos_end, pos_start, vel_end, vel_start, vel_lim, acc_min, acc_max);

    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << (duration.count() / 1000.0)<< std::endl;

    VectorXf pos_plot = prof.pos[0];
    VectorXf vel_plot = prof.vel[0];
    VectorXf acc_plot = prof.acc[0];

    // // plt::subplot(4,1,2);
    // plt::plot(std::vector<double>(prof.time.data(), prof.time.data()+prof.time.size()), std::vector<float>(pos_plot.data(), pos_plot.data()+pos_plot.size()), "");

    // // plt::subplot(4,1,3);
    // plt::plot(std::vector<double>(prof.time.data(), prof.time.data()+prof.time.size()), std::vector<float>(vel_plot.data(), vel_plot.data()+vel_plot.size()), "");

    // // plt::subplot(4,1,4);
    // plt::plot(std::vector<double>(prof.time.data(), prof.time.data()+prof.time.size()), std::vector<float>(acc_plot.data(), acc_plot.data()+acc_plot.size()), "");


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


    plt::show();
}
