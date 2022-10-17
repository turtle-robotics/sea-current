#include "../sea_current.hpp"

#include <Eigen/Dense>

using namespace Eigen;
using namespace turtle::sc;

int main() {
    turtle::sc::dbg_assert(true, "");

    std::vector<Vector2f> ctrl_pts;

    ctrl_pts.push_back(Vector2f(-1.0, 0.0));
    ctrl_pts.push_back(Vector2f(-0.5, -0.5));
    ctrl_pts.push_back(Vector2f(0.5, 0.5));
    ctrl_pts.push_back(Vector2f(1.0, 0.0));

    bezier_spline bs = bezier_spline::bezier_curve(ctrl_pts, 0.01);
}
