#include "../sea_current.hpp"

#include "../third-party/matplotlib-cpp/matplotlibcpp.h"

#include <Eigen/Dense>

#include <chrono>
#include <iostream>

using namespace Eigen;
using namespace turtle::sc;

namespace plt = matplotlibcpp;

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
    plt::plot(x, y, "");
    plt::plot(x2, y2, "ro");
    plt::plot(x, y_hat_p, "");
    plt::show();
}
