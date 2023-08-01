#include "../sea_current.hpp"

#include <Eigen/Dense>

#include <zmq.hpp>

#include <string>
#include <sstream>
#include <cstring>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

using namespace Eigen;
using namespace turtle::sc;

int main() {
    zmq::context_t ctx;
    zmq::socket_t socket(ctx, zmq::socket_type::rep);

    // socket.bind("tcp://localhost:8080");
    socket.bind("ipc:///tmp/test");


    while (true) {
        zmq::message_t request;
        const bool res = socket.recv(&request);
        std::string req = std::string(static_cast<char*>(request.data()), request.size());

        std::stringstream ss(req);

        float acc_min_val;
        float acc_max_val;
        float vel_min_val;
        float vel_max_val;
        ss >> acc_min_val;
        ss >> acc_max_val;
        ss >> vel_min_val;
        ss >> vel_max_val;

        // std::cout << acc_min_val << std::endl;
        // std::cout << acc_max_val << std::endl;
        // std::cout << vel_min_val << std::endl;
        // std::cout << vel_max_val << std::endl;

        std::vector<Vector2f> path;
        float max_x = 0;
        float max_y = 0;
        while (ss.tellp() == std::streampos(0)) {
            float x, y;
            ss >> x;
            ss >> y;
            // std::cout << "x y " << x << " " << y << std::endl;
            if (std::abs(x) > max_x) max_x = std::abs(x);
            if (std::abs(y) > max_y) max_y = std::abs(y);
            path.push_back(Vector2f(x, y));
        }

        path.pop_back();

        const bounding_rect br = {max_x, -max_x, max_y, -max_y};
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

        // std::cout << prof.pos[0][1] << std::endl;
        bezier_spline re = pad.resample(prof.pos[0], ad, true);

        const std::vector<float> ang_vel = re.angular_velocity(prof);

        json j = serialize_path_to_json(re, prof, ad, ang_vel);
        std::stringstream js;
        js << std::setw(4) << j << std::endl;
        std::string rep = js.str();
        zmq::message_t reply(rep);
        socket.send(reply, zmq::send_flags::none);
    }
    return 0;
}
