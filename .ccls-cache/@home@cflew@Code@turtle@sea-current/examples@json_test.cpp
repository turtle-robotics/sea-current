#include "../sea_current.hpp"

#include <Eigen/Dense>

#include <fstream>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

using namespace Eigen;
using namespace turtle::sc;

int main() {
    const bounding_rect br = {1, -1, 1, -1};
    // const bounding_rect br = {1, 0, 1, 0};
    planning_space space(br);


    // std::vector<Vector2f> path = {Vector2f(0, 0), Vector2f(1, 1), Vector2f(0, 2)};
    std::vector<Vector2f> path = {
    Vector2f(0, 0),
    Vector2f(10, 0),
    Vector2f(10, 10),
    Vector2f(10, 20),
    Vector2f(30, 50),
    Vector2f(0, 0)
    };
    bezier_spline pad = bezier_spline::from_path(path, space);

    const arclength_data ad = pad.arclength();
    const Eigen::Vector<value_type, 1> pos_end{ad.arclength};
    const Eigen::Vector<value_type, 1> pos_start{0};

    const Eigen::Vector<value_type, 1> vel_end{0};
    const Eigen::Vector<value_type, 1> vel_start{0};

    const Eigen::Vector<value_type, 1> acc_min{-0.25};
    const Eigen::Vector<value_type, 1> acc_max{0.25};

    auto vel_lim = [](toppra::value_type time) {
        toppra::Vector lower{1};
        toppra::Vector upper{1};

        lower(0, 0) = 0;
        upper(0, 0) = 1;

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

    std::ofstream o("test6.json");
    json j = serialize_path_to_json(re, prof, ad, ang_vel);
    o << std::setw(4) << j << std::endl;
}
