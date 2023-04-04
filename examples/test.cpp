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

// inline Vector2f calc_start_tangent(const Vector2f& W_0, const Vector2f& W_1, const float theta);


// inline float tangent_magnitude(const Vector2f& W_0, const Vector2f& W_1, const Vector2f& W_2) {
//     return 0.5 * std::min(pt_dist(W_0, W_1), pt_dist(W_1, W_2));
// }

// inline Vector2f calc_tangent(const Vector2f& W_0, const Vector2f& W_1, const Vector2f& W_2) {
//     const Vector2f u = W_0 - W_1;
//     const Vector2f v = W_2 - W_1;
//     // std::cout << "u " << u.x() << " " << u.y() << std::endl;
//     // std::cout << "v " << v.x() << " " << v.y() << std::endl;
//     float u_dot_v = u.dot(v);
//     const float denom = pt_dist(u) * pt_dist(v);
//     const float theta = std::acos((u_dot_v) / denom);
//     // std::cout << (u_dot_v) / denom << std::endl;
//     // std::cout << "theta " << theta << std::endl;

//     return tangent_magnitude(W_0, W_1, W_2) * Vector2f(std::cos(theta), std::sin(theta));

//     // return calc_start_tangent(W_1, W_2, theta);
// }

// inline Vector2f calc_start_tangent(const Vector2f& W_0, const Vector2f& W_1, const float theta) {
//     return tangent_magnitude(W_0, W_1, W_0) * Vector2f(std::cos(theta), std::sin(theta));
// }

// inline Vector2f calc_end_tangent(const Vector2f& W_1, const Vector2f W_2) {
//     return tangent_magnitude(W_1, W_2, W_1) * (W_2 - W_1).normalized();
// }

int main() {

    std::vector<Vector2f> ctrl_pts;
    std::vector<Vector2f> ctrl_pts2;

    const Vector2f W_0(0, 0);
    const Vector2f W_1(0.5, 0.5);
    const Vector2f W_2(1, 0);

    const Vector2f T_0 = calc_start_tangent(W_0, W_1, 0);
    const Vector2f T_1 = calc_tangent(W_0, W_1, W_2);
    const Vector2f T_2 = calc_end_tangent(W_1, W_2);

    constexpr float k = 0.2;
    std::cout << "T_0 " << T_0.x() << " " << T_0.y() << std::endl;
    std::cout << "T_1 " << T_1.x() << " " << T_1.y() << std::endl;
    std::cout << "T_2 " << T_2.x() << " " << T_2.y() << std::endl;

    ctrl_pts.push_back(W_0);
    ctrl_pts.push_back(W_0 + (k) * T_0);
    ctrl_pts.push_back(W_1 - (k) * T_1);
    ctrl_pts.push_back(W_1);

    ctrl_pts2.push_back(W_1);
    ctrl_pts2.push_back(W_1 + (k) * T_1);
    ctrl_pts2.push_back(W_2 - (k) * T_2);
    ctrl_pts2.push_back(W_2);

    // ctrl_pts.push_back(Vector2f(-1.0, 0.0));
    // ctrl_pts.push_back(Vector2f(-0.5, 0.5));
    // ctrl_pts.push_back(Vector2f(0.5, -0.5));
    // ctrl_pts.push_back(Vector2f(1.0, 0.0));

    for (auto& ctrl_pt : ctrl_pts) {
        std::cout << "ctrl_pt: " << ctrl_pt.x() << " " << ctrl_pt.y() << std::endl;
    }

    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "die?" << std::endl;
    bezier_spline bs = bezier_spline::bezier_curve(ctrl_pts, 0.0001);
    bezier_spline bs2 = bezier_spline::bezier_curve(ctrl_pts2, 0.0001);
    std::cout << "die?" << std::endl;
    bs = bezier_spline::join_splines({bs, bs2});
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
    // std::cout << "here";
    std::cout << "arclength: " << arclen << std::endl;
    std::cout << "arclength2: " << bs2.arclength(0.01).arclength << std::endl;
    // std::cout << "here";
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
    std::cout << "die?" << std::endl;
    bezier_spline re = bs.resample(pos_plot, ad, true);
    std::cout << "die?" << std::endl;

    std::vector<float> x3(re.n_pts());
    std::vector<float> y3(re.n_pts());
    for (int i = 0; i < re.n_pts(); ++i) {
        x3[i] = re.pts(i, 0);
        y3[i] = re.pts(i, 1);
        // std::cout << x3[i] << " " << y3[i] << std::endl;
    }

    // plt::plot(x3, y3, "");

    start = std::chrono::high_resolution_clock::now();

    const bounding_rect br = {1, -1, 1, -1};
    // const bounding_rect br = {1, 0, 1, 0};
    planning_space space(br);

    obstacle ob({Vector2f(-0.5, 0), Vector2f(1, 0), Vector2f(1, 1), Vector2f(0, 1)});
    obstacle ob2({Vector2f(0, -0.5), Vector2f(1, 0), Vector2f(1, 1), Vector2f(0, 1)});
    obstacle ob3({Vector2f(-0.6, 0.148), Vector2f(-1, 0.148), Vector2f(-1, 0), Vector2f(-0.6, 0)});
    space.obstacles.push_back(ob);
    space.obstacles.push_back(ob2);
    space.obstacles.push_back(ob3);

    point_set pts_set = space.sample_free(2560);
    std::vector<Vector2f> pts(pts_set.begin(), pts_set.end());

    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << (duration.count() / 1000.0)<< std::endl;

    // for (auto& pt : pts) {
    //     std::cout << "pt: " << pt.x() << " " << pt.y() << std::endl;
    // }

    std::cout << "c: " << ob.contains(Vector2f(0.5, 0.5)) << std::endl;

    std::vector<float> x4(pts.size());
    std::vector<float> y4(pts.size());
    for (int i = 0; i < pts.size(); ++i) {
        // std::cout << "x y: " << pts[i].x() << " " << pts[i].y() << std::endl;
        x4[i] = pts[i].x();
        y4[i] = pts[i].y();
        // std::cout << x3[i] << " " << y3[i] << std::endl;
    }

    start = std::chrono::high_resolution_clock::now();
    std::vector<Vector2f> path = space.fast_marching_trees(Vector2f(-0.5, 1), Vector2f(1, -1), 200, 1).value();
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << (duration.count() / 1000.0)<< std::endl;

    std::cout << path.size() << std::endl;
    std::vector<float> x5;
    std::vector<float> y5;
    x5.reserve(path.size());
    y5.reserve(path.size());
    for (const auto& p : path) {
        x5.push_back(p.x());
        y5.push_back(p.y());
        // std::cout << p.x() << " " << p.y() << std::endl;
    }

    bezier_spline pad = bezier_spline::from_path(path, space);
    std::vector<float> x6(pad.n_pts());
    std::vector<float> y6(pad.n_pts());
    for (int i = 0; i < pad.n_pts(); ++i) {
        x6[i] = pad.pts(i, 0);
        y6[i] = pad.pts(i, 1);
        // std::cout << "x: " << x6[i] << " y: " << y6[i] << std::endl;
    }

    std::vector<float> x7;
    std::vector<float> y7;
    for (int i = 0; i < pad.ctrl_pts[0].size(); ++i) {
        // x7[i] = pad.ctrl_pts[0][i].x();
        // y7[i] = pad.ctrl_pts[0][i].y();
        for (int j = 0; j < pad.ctrl_pts.size() && j < 1; ++j) {
            x7.push_back(pad.ctrl_pts[j][i].x());
            y7.push_back(pad.ctrl_pts[j][i].y());
        }
    }

    plt::plot(x4, y4, "o");

    plt::plot(x5, y5);
    plt::plot(x6, y6);
    plt::plot(x7, y7, "o");

    plt::show();
}
