#include "../sea_current.hpp"

#include <Eigen/Dense>

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <cmath>


using namespace Eigen;
using namespace turtle::sc;

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: json_util [input file] [output file]" << std::endl;
        return -1;
    }

    std::ifstream istrm(std::string(argv[1]));

    if (!istrm.is_open()) {
        std::cerr << "failed to open file " << std::endl;
        return -1;
    }

    std::string fline;
    std::getline(istrm, fline);

    std::stringstream ss(fline);

    float acc_min_val;
    float acc_max_val;
    float vel_min_val;
    float vel_max_val;
    ss >> acc_min_val;
    ss >> acc_max_val;
    ss >> vel_min_val;
    ss >> vel_max_val;

    std::vector<Vector2f> path;
    float max_x = 0;
    float max_y = 0;
    for (std::string line; std::getline(istrm, line); ) {
        float x, y;
        ss >> x;
        ss >> y;
        if (std::abs(x) > max_x) max_x = std::abs(x);
        if (std::abs(y) > max_y) max_y = std::abs(y);
        path.push_back(Vector2f(x, y));
    }


    const bounding_rect br = {max_x, -max_x, max_y, -max_y};
    // const bounding_rect br = {1, 0, 1, 0};
    planning_space space(br);


    // std::vector<Vector2f> path = {Vector2f(0, 0), Vector2f(1, 1), Vector2f(0, 2)};
    bezier_spline pad = bezier_spline::from_path(path, space);

    const arclength_data ad = pad.arclength();
    const Eigen::Vector<value_type, 1> pos_end{ad.arclength};
    const Eigen::Vector<value_type, 1> pos_start{0};

    const Eigen::Vector<value_type, 1> vel_end{0};
    const Eigen::Vector<value_type, 1> vel_start{0};

    const Eigen::Vector<value_type, 1> acc_min{acc_min_val};
    const Eigen::Vector<value_type, 1> acc_max{acc_max_val};

    auto vel_lim = [&](toppra::value_type time) {
        toppra::Vector lower{1};
        toppra::Vector upper{1};

        lower(0, 0) = vel_min_val;
        upper(0, 0) = vel_max_val;

        return std::make_tuple(lower, upper);
    };

    velocity_profile prof = gen_vel_prof<1>(pos_end, pos_start, vel_end, vel_start, vel_lim, acc_min, acc_max);

    std::cout << prof.pos[0][1] << std::endl;
    bezier_spline re = pad.resample(prof.pos[0], ad, true);

    const std::vector<float> ang_vel = re.angular_velocity(prof);

    // const std::vector<float> curv = re.curvature();

    // std::cout << "------------" << std::endl;
    // for (std::size_t i = 0; i < curv.size(); ++i) {
    //     std::cout << curv[i] << std::endl;
    // }
    // std::cout << "------------" << std::endl;

    serialize_path_to_file("test6.json", re, prof, ad, ang_vel);
}
