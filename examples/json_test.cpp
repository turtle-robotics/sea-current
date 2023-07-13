#include "../sea_current.hpp"

#include <Eigen/Dense>


using namespace Eigen;
using namespace turtle::sc;

int main() {
    const bounding_rect br = {1, -1, 1, -1};
    // const bounding_rect br = {1, 0, 1, 0};
    planning_space space(br);


    std::vector<Vector2f> path = {Vector2f(0, 0), Vector2f(1, 0)};
    bezier_spline pad = bezier_spline::from_path(path, space);

    arclength_data ad = pad.arclength();
    Eigen::Vector<value_type, 1> pos_end{ad.arclength};
    Eigen::Vector<value_type, 1> pos_start{0};

    Eigen::Vector<value_type, 1> vel_end{0};
    Eigen::Vector<value_type, 1> vel_start{0};

    Eigen::Vector<value_type, 1> acc_min{-40};
    Eigen::Vector<value_type, 1> acc_max{40};

    auto vel_lim = [](toppra::value_type time) {
        toppra::Vector lower{1};
        toppra::Vector upper{1};

        lower(0, 0) = 0;
        upper(0, 0) = 5;

        return std::make_tuple(lower, upper);
    };

    velocity_profile prof = gen_vel_prof<1>(pos_end, pos_start, vel_end, vel_start, vel_lim, acc_min, acc_max);

    std::cout << prof.pos[0][1] << std::endl;
    bezier_spline re = pad.resample(prof.pos[0], ad, true);

    serialize_path_to_file("test.json", re, prof);
}
