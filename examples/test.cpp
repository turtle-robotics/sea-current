#include "../sea_current.hpp"

#include "../third-party/matplotlib-cpp/matplotlibcpp.h"

#include <Eigen/Dense>

#include <chrono>
#include <iostream>

using namespace Eigen;
using namespace turtle::sc;

namespace plt = matplotlibcpp;

int main() {
    turtle::sc::dbg_assert(true, "");

    std::vector<Vector2f> ctrl_pts;

    ctrl_pts.push_back(Vector2f(-1.0, 0.0));
    ctrl_pts.push_back(Vector2f(-0.5, 0.5));
    ctrl_pts.push_back(Vector2f(0.5, -0.5));
    ctrl_pts.push_back(Vector2f(1.0, 0.0));

    auto start = std::chrono::high_resolution_clock::now();
    bezier_spline bs = bezier_spline::bezier_curve(ctrl_pts, 0.01);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << (duration.count() / 1000.0)<< std::endl;

    std::vector<float> x(bs.num_pts());
    std::vector<float> y(bs.num_pts());
    for (int i = 0; i < bs.num_pts(); ++i) {
        x[i] = bs.pts(i, 0);
        y[i] = bs.pts(i, 1);
    }

    std::vector<float> x2(bs.num_pts());
    std::vector<float> y2(bs.num_pts());
    for (int i = 0; i < bs.ctrl_pts[0].size(); ++i) {
        x2[i] = bs.ctrl_pts[0][i].x();
        y2[i] = bs.ctrl_pts[0][i].y();
    }


    plt::plot(x, y, "");
    plt::plot(x2, y2, "ro");
    plt::show();
}
