#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <numbers>
#include <algorithm>
#include <tuple>
#include <functional>
#include <queue>
#include <limits>
#include <unordered_set>
#include <unordered_map>
#include <optional>
#include <string>

#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

#include <toppra/toppra.hpp>
#include <toppra/geometric_path.hpp>
#include <toppra/geometric_path/piecewise_poly_path.hpp>
#include <toppra/constraint/linear_joint_velocity.hpp>
#include <toppra/constraint/linear_joint_acceleration.hpp>
#include <toppra/algorithm/toppra.hpp>
#include <toppra/parametrizer.hpp>
#include <toppra/parametrizer/spline.hpp>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#ifdef DEBUG
#define SC_ASSERT(cnd, msg)                                                           \
    do {                                                                              \
        static_assert(                                                                \
                !std::is_pointer_v<decltype(cnd)>,                                    \
                "Do not use SC_ASSERT with raw pointers"                              \
                "and instead do SC_ASSERT(cnd != nullptr) "                           \
                "to avoid implicit pointer-to-bool conversion.");                     \
        if (bool(cnd) == false) {                                                     \
            std::cerr << "SC_ASSERT failed: " << #cnd << std::endl;                   \
            std::cerr << "SC_ASSERT message: " << msg << std::endl;                   \
            std::cerr << "SC_ASSERT failed at: ";                                     \
            std::cerr << __func__ << " " << __FILE__ << " " << __LINE__ << std::endl; \
            std::exit(1);                                                             \
        }                                                                             \
    } while(0)
#else
#define SC_ASSERT(cnd, msg)
#endif


namespace turtle::sc {

    using namespace Eigen;

    using namespace std::complex_literals;

    struct LIDAR_point{ // struct to handle a LIDAR point
        double x;
        double y;

        Point(double xCoord, double yCoord) : x(xCoord), y(yCoord) {}
    }

    std::vector<Point> points;  // creates variable called 'points' of type LIDAR_point

    struct bounding_rect {  // struct creates bounding rectangle that is initially of size
        float x_max;
        float x_min;
        float y_max;
        float y_min;

        inline bool contains(const Vector2f& point) {   // contains method checks is 'point' is inside bounding rectangle (input: ,output: bool)
            return point.x() <= x_max &&
                   point.x() >= x_min &&
                   point.y() <= y_max &&
                   point.y() >= y_min;
        }

        inline void enclose_point(const Vector2f& pt) {    // enclose_point method will check if the point is inside bounding rectangle, if not then the bounding rectangle changes size to accommodate
            if (pt.x() > x_max) x_max = pt.x();
            if (pt.x() < x_min) x_min = pt.x();
            if (pt.y() > y_max) y_max = pt.y();
            if (pt.y() < y_min) y_min = pt.y();
        }

        bounding_rect(const float x_max, const float x_min, const float y_max, const float y_min) : x_max(x_max), x_min(x_min), y_max(y_max), y_min(y_min) {}
    };

    inline float pt_dist(const Vector2f& u, const Vector2f& v={0,0}) {  // pt_dist function returns the distance between points 'u' and 'v'
        return std::sqrt(std::pow(v.x() - u.x(), 2) + std::pow(v.y() - u.y(), 2));
    }

    struct halton_state {   // halton_state struct will set 'f' and 'i' to values it is sent
        int f = 0;
        int i = 0;
        halton_state(int f, int i) : f(f), i(i) {}
        halton_state() : f(0), i(0) {}
    };

    std::vector<float> halton(const int b, const int n, halton_state& state) {  // halton function creates a vector 'nums' that contain a halton seqence
        std::vector<float> nums(n);

        int f = 0;
        int i = 1;

        if (state.i != 0 && state.f != 0) {
            i = state.i;
            f = state.f;
        }

        for (int j = 0; j < n; ++j) {

            const int x = i - f;
            if (x == 1) {
                f = 1;
                i *= b;
            } else {
                int y = std::floor(i / b);
                while (x <= y) {
                    y = std::floor(y / b);
                }
                f = (b + 1) * y - x;
            }

            nums[j] = static_cast<float>(f) / i;
        }

        state.f = f;
        state.i = i;

        return nums;
    }



    inline float cross2d(Vector2f u, Vector2f v) {      // cross2d function returns the magnitude of the cross product of 'u' and 'v'
        return u.x()*v.y() - u.y()*v.x();
    }

    std::tuple<bool, Vector2f> intersects(std::tuple<Vector2f, Vector2f> l, std::tuple<Vector2f, Vector2f> k) {
        const float a = cross2d(std::get<0>(k)-std::get<0>(l), std::get<1>(l)-std::get<0>(l));
        const float b = cross2d(std::get<1>(l)-std::get<0>(l), std::get<1>(k)-std::get<0>(k));

        if (a == 0 && b == 0) {
            // float ax = std::get<0>(l).x();
            // float bx = std::get<1>(l).x();
            // float cx = std::get<0>(k).x();
            // float dx = std::get<1>(k).x();

            // using std::max;
            // using std::min;

            // return (max(ax, bx) >= min(cx, dx) && min(ax, bx) <= min(cx, dx))
            //     || (max(cx, dx) >= min(ax, bx) && min(cx, dx) <= min(ax, bx))
            //     || (min(ax, bx) <= max(cx, dx) && max(ax, bx) >= max(cx, dx))
            //     || (min(cx, dx) <= max(ax, bx) && max(cx, dx) >= max(ax, bx));

            // TODO: consider a different way of handling colinear segments
            return {false, Vector2f(0, 0)};
        } else if (b == 0 && a != 0) {
            return {false, Vector2f(0, 0)};
        } else if (b != 0) {
            const float u = a / b;

            const float c = cross2d(std::get<0>(k)-std::get<0>(l), std::get<1>(k)-std::get<0>(k));
            // float d = cross2d(std::get<1>(l)-std::get<0>(l), std::get<1>(k)-std::get<0>(k));
            const float d = b;

            const float t = c / d;

            if (0 <= u && u <= 1 && 0 <= t && t <= 1) {
                return {true, t*(std::get<1>(l)-std::get<0>(l))+std::get<0>(l)};
            }
        }
        return {false, Vector2f(0, 0)};
    }
    struct obstacle {
        std::vector<std::tuple<Vector2f, Vector2f>> lines;
        // create function that turns points into lines
        std::vector<Vector2f> vertices;
        std::bool closed_obstacle;    // Creates a bool to store whether or not the obstacle is open or closed -CB
        bounding_rect bound_rect = {0, 0, 0, 0};

        bool contains(const Vector2f& point) {
            if (!bound_rect.contains(point)) return false;

            const Vector2f outside_pt(bound_rect.x_max+1, bound_rect.y_max+1);

            int count = 0;

            std::vector<bool> corners(vertices.size());
            for (auto& line : lines) {
                auto [ints, int_pt] = intersects({outside_pt, point}, line);
                if (ints) {
                    bool counted = false;

                    for (int i = 0; i < vertices.size(); ++i) {
                        using std::abs;
                        if (abs(int_pt.x() - vertices[i].x()) <= 0.00001 &&
                            abs(int_pt.y() - vertices[i].y()) <= 0.00001) {

                            if (corners[i]) {
                                counted = true;
                                continue;
                            }


                            corners[i] = true;
                            count++;
                            counted = true;
                            break;
                        }
                    }

                    if (!counted) count++;
                }
            }


            return count % 2 == 1;
        }

        obstacle(std::vector<Vector2f> vertices) : vertices(vertices) {
            bound_rect = {vertices[0].x(), vertices[0].x(), vertices[0].y(), vertices[0].y()};

            if (vertices.size() % 2 == 0) {
                vertices.push_back(vertices[0]);
            }

            lines.resize(vertices.size()-1);

            for (int i = 0; i < vertices.size() - 1; ++i) {
                lines[i] = {vertices[i], vertices[i+1]};

                bound_rect.enclose_point(vertices[i]);
            }
        }

        obstacle(std::vector<Vector2f> vertices, std::vector<std::tuple<int, int>> edges) : vertices(vertices) {
            bound_rect = {vertices[0].x(), vertices[0].x(), vertices[0].y(), vertices[0].y()};
            lines.reserve(edges.size());
            for (const auto& edge : edges) {
                SC_ASSERT(std::get<0>(edge) < vertices.size(), "edge must reference indices within the vertex list");
                SC_ASSERT(std::get<1>(edge) < vertices.size(), "edge must reference indices within the vertex list");

                lines.push_back({vertices[std::get<0>(edge)], vertices[std::get<1>(edge)]});

                bound_rect.enclose_point(vertices[std::get<0>(edge)]);
                bound_rect.enclose_point(vertices[std::get<1>(edge)]);
            }
        }
    };

    struct hash_vector2f {  // hash_vector2f is a struct that converts a float 2d vector to a single hash value
        size_t operator()(const Vector2f v) const {
            float fa = v.x();
            float fb = v.y();
            const int32_t a = reinterpret_cast<int32_t&>(fa);
            const int32_t b = reinterpret_cast<int32_t&>(fb);
            return std::hash<int32_t>()(a) ^ std::hash<int32_t>()(b);
        }
    };
    using point_set = std::unordered_set<Vector2f, hash_vector2f>;

    class planning_space {
        public:
            point_set sample_free(const int n);
            float cost(const Vector2f a, const Vector2f b) const;
            point_set near(const Vector2f b, const point_set& nodes, const float dist) const;
            std::optional<std::vector<Vector2f>> fast_marching_trees(const Vector2f& x_init, const Vector2f& x_goal, const int n, const float rn);
            void register_free_space_constraint(std::function<bool(Vector2f)>);

            planning_space(const bounding_rect& br);

            std::vector<obstacle> obstacles;
            bounding_rect bound_rect;
            halton_state x_state;
            halton_state y_state;
    };

    struct arclength_data {
        float arclength;
        std::vector<VectorXf> segments;
        std::vector<VectorXf> positions;
    };

    struct chebpoly {
        VectorXf coeffs;
        float xmin;
        float xmax;

        chebpoly(VectorXf coeffs, const float xmin, const float xmax) : coeffs(coeffs), xmin(xmin), xmax(xmax) {}
    };

    chebpoly chebfit(const VectorXf& x, const VectorXf& y, const int degree);
    VectorXf chebeval(const VectorXf& x, const chebpoly& b, const int degree);

    inline Vector2f calc_start_tangent(const Vector2f& W_0, const Vector2f& W_1, const float theta);


    inline float tangent_magnitude(const Vector2f& W_0, const Vector2f& W_1, const Vector2f& W_2) {
        return 0.5 * std::min(pt_dist(W_0, W_1), pt_dist(W_1, W_2));
    }

    inline Vector2f calc_tangent(const Vector2f& W_0, const Vector2f& W_1, const Vector2f& W_2) {
        // std::cout << "calc_tangent call ------------------" << std::endl;
        // std::cout << "W_0 " << W_0.x() << " " << W_0.y() << std::endl;
        // std::cout << "W_1 " << W_1.x() << " " << W_1.y() << std::endl;
        // std::cout << "W_2 " << W_2.x() << " " << W_2.y() << std::endl;
        const Vector2f u = W_0 - W_1;
        const Vector2f v = W_2 - W_1;
        // std::cout << "u " << u.x() << " " << u.y() << std::endl;
        // std::cout << "v " << v.x() << " " << v.y() << std::endl;
        float u_dot_v = u.dot(v);
        const float denom = pt_dist(u) * pt_dist(v);
        const float theta = std::acos((u_dot_v) / denom) / 2;

        // const float offset = std::atan2(-u.y(), -u.x());
        const float offset = std::atan2(u.y(), u.x());
        const float test_offset = std::atan2(v.y(), v.x());

        int mult = 1;
        if (test_offset - offset < 0) mult = -1;

        // std::cout << "theta " << (theta * 180 / 3.1415926535) << std::endl;
        // std::cout << "offset " << (offset * 180 / 3.1415926535) << std::endl;

        // std::cout << "offset " << offset << std::endl;

        Vector2f l90(std::sin(offset + (mult * theta)), -std::cos(offset + (mult * theta)));
        l90 = l90.normalized();
        // Vector2f r90 = -l90;

        // std::cout << "calc_tangent end ????????????" << std::endl;
        int mult2 = -1;
        if (pt_dist(W_1+l90, W_2) < pt_dist(W_1+(-l90), W_2)) {
            // return tangent_magnitude(W_0, W_1, W_2) * l90;
            mult2 = 1;
        }
        return tangent_magnitude(W_0, W_1, W_2) * (mult2 * l90);

        // // std::cout << (u_dot_v) / denom << std::endl;
        // std::cout << "theta " << (180 * theta / 3.1415926535) << std::endl;

        // return tangent_magnitude(W_0, W_1, W_2) * Vector2f(std::cos(theta), std::sin(theta));

        // // return calc_start_tangent(W_1, W_2, theta);

    }

    inline Vector2f calc_start_tangent(const Vector2f& W_0, const Vector2f& W_1, const float theta) {
        // Vector2f tmp(std::cos(theta), std::sin(theta));
        // std::cout << "wtfff " << tmp.x() << " " << tmp.y() << std::endl;
        // std::cout << tangent_magnitude(W_0, W_1, W_0) << std::endl;
        // Vector2f tmp2 = tangent_magnitude(W_0, W_1, W_0) * Vector2f(std::cos(theta), std::sin(theta));
        // std::cout << "wtfff2 " << tmp2.x() << " " << tmp2.y() << std::endl;
        return tangent_magnitude(W_0, W_1, W_0) * Vector2f(std::cos(theta), std::sin(theta));
    }

    inline Vector2f calc_end_tangent(const Vector2f& W_1, const Vector2f W_2) {
        return tangent_magnitude(W_1, W_2, W_1) * (W_2 - W_1).normalized();
    }

    struct velocity_profile {
        std::vector<VectorXf> pos;
        std::vector<VectorXf> vel;
        std::vector<VectorXf> acc;
        toppra::Vector time;

        velocity_profile(std::vector<VectorXf> pos, std::vector<VectorXf> vel, std::vector<VectorXf> acc, toppra::Vector time) : pos(pos), vel(vel), acc(acc), time(time) {}
    };

    class bezier_spline {
        public:
            std::vector<std::vector<Vector2f>> ctrl_pts;
            Matrix<float, Dynamic, 2> pts;
            std::vector<VectorXf> positions; // for curvature/hodograph
            std::vector<Matrix<std::complex<float>, Dynamic, 2>> Q_cache;

            bezier_spline() = default;
            bezier_spline(const std::vector<std::vector<Vector2f>>& ctrl_pts, const Matrix<float, Dynamic, 2>& pts, const std::vector<VectorXf>& positions);
            bezier_spline(const std::vector<std::vector<Vector2f>>& ctrl_pts, const std::vector<VectorXf>& positions); // uniform positions for each curve

            static bezier_spline join_splines(const std::vector<bezier_spline>& splines);

            static Vector2f shrink_tangent(const Vector2f& T, const Vector2f& W, const float k, const planning_space& ps);
            static bezier_spline from_path(const std::vector<Vector2f>& path, const planning_space& ps, float start_angle);

            static bezier_spline bezier_curve(const std::vector<Vector2f>& ctrl_pts, const std::vector<float>& positions);
            static bezier_spline bezier_curve(const std::vector<Vector2f>& ctrl_pts, const VectorXf& positions);
            static bezier_spline bezier_curve(const std::vector<Vector2f>& ctrl_pts, const VectorXf& positions, Matrix<std::complex<float>, Dynamic, 2> Q);
            static bezier_spline bezier_curve(std::vector<Vector2f>& ctrl_pts, float precision);

            inline int n_pts() const;
            inline int n_segments() const;
            inline int degree() const;

            arclength_data arclength(const float precision) const;
            std::vector<float> curvature() const;
            bezier_spline hodograph() const;
            std::vector<float> angular_velocity(const velocity_profile& vel_prof) const;
            std::vector<float> angular_velocity2(const velocity_profile& vel_prof) const;

            bezier_spline resample(VectorXf& profile_pos, arclength_data ad, bool nudge_positions) const;

            // helper function for bezier_curve
            static inline std::vector<std::complex<float>> omega_table(int degree);

        private:
            template <typename T, typename Y>
            static inline T coerce(const T& num, const Y& low, const Y& high) {
                if (num < low) return low;
                else if (num > high) return high;
                return num;
            }
    };

    class planner {
        public:
            planning_space ps;

            Vector2f pick_next_goal_point() const;
            // void update_planning_space();
            for (const obstacle& obs : obsacles) {   // Loop that will loop through all of the obstacles -CB
                if (obs.closed_obstacle) {  // check to see if obstacle is closed   -CB
                    continue;   // moves to the next iteration (if obstacle is closed)  -CB
                }
                // look for end point to move towards   -CB
                // pick closest end point to robots location    -CB
                // determine "look-around" point that is in line with the robot and the endpont, and is past the end point  -CB
                // be sure that "look-around" point is within the area scanned by lidar (prevents blind collisions) -CB
                // move  -CB
                // once moved scan  -CB
                // create a line between "close" scanned points (close should be some value slightly larger than how close the lidar will place points on an object -CB   
                // if scanned points are not "close" enough add them to stack for potential objects -CB
                // repet this process until scanned end point is "close" enough to a previously stored end point   -CB
                // update the object to closed (obs.closed_obstacle = true;), ready to move to next iteration of loop   -CB
            }
    };

    Vector2f planner::pick_next_goal_point() const {

    }


    bezier_spline::bezier_spline(const std::vector<std::vector<Vector2f>>& ctrl_pts, const Matrix<float, Dynamic, 2>& pts, const std::vector<VectorXf>& positions) : ctrl_pts(ctrl_pts), pts(pts), positions(positions) {}

    inline int bezier_spline::n_pts() const {
        return pts.rows();
    }

    inline int bezier_spline::n_segments() const {
        return ctrl_pts.size();
    }

    inline int bezier_spline::degree() const {
        return ctrl_pts[0].size()-1;
    }

    // it's ok that this is static as a new matrix would have to be allocated either way
    bezier_spline bezier_spline::join_splines(const std::vector<bezier_spline>& splines) {
        std::vector<std::vector<Vector2f>> ctrl_pts;
        std::vector<VectorXf> positions;
        std::vector<Matrix<std::complex<float>, Dynamic, 2>> Qs;

        int n_pts = 0;
        for (bezier_spline spline : splines) {
            n_pts += spline.n_pts();
            for (std::size_t i = 0; i < spline.n_segments(); ++i) {
                ctrl_pts.push_back(spline.ctrl_pts[i]);
                positions.push_back(spline.positions[i]);
                Qs.push_back(spline.Q_cache[i]);
            }
        }

        Matrix<float, Dynamic, 2> joined_pts = Matrix<float, Dynamic, 2>::Zero(n_pts, 2);

        int n = 0;
        for (bezier_spline spline : splines) {
            joined_pts.block(n, 0, spline.n_pts(), 2) = spline.pts;
            n += spline.n_pts();
        }

        bezier_spline bs = bezier_spline(ctrl_pts, joined_pts, positions);
        bs.Q_cache = Qs;
        return bs;
    }

    inline Vector2f bezier_spline::shrink_tangent(const Vector2f& T, const Vector2f& W, const float k, const planning_space& ps) {
        Vector2f ret_pt = k * T;
        // std::cout << "ret_pt " << ret_pt.x() << " " << ret_pt.y() << std::endl;
        // std::cout << "W " << W.x() << " " << W.y() << std::endl;
        for (auto& obstacle : ps.obstacles) {
            for (auto& oline : obstacle.lines) {
                auto [ints, int_pt] = intersects({W+ret_pt, W}, oline);
                if (ints) {
                    ret_pt = int_pt - W;
                }
                auto [ints2, int_pt2] = intersects({W-ret_pt, W}, oline);
                if (ints2) {
                    // std::cout << "ints2 " << int_pt2.x() << " " << int_pt2.y() << std::endl;
                    // ret_pt = -(int_pt2 - W);
                    ret_pt = W - int_pt2;
                }
                // std::cout << "ret_pt " << ret_pt.x() << " " << ret_pt.y() << std::endl;
            }
        }

        return ret_pt;
    }

    bezier_spline bezier_spline::from_path(const std::vector<Vector2f>& path, const planning_space& ps, float start_angle=NAN) {
        SC_ASSERT(path.size() >= 2, "Not enough points for a path");

        // std::vector<bezier_spline> parts;
        std::vector<Vector2f> tangents;
        tangents.reserve(path.size());

        if (std::isnan(start_angle)) {
            const Vector2f tmp_diff = path[1] - path[0];
            start_angle = std::atan2(tmp_diff.y(), tmp_diff.x());
            // std::cout << "start angle " << (start_angle * 180 / 3.1415926) << std::endl;
        }

        constexpr float k = 1;

        Vector2f T_0 = calc_start_tangent(path[0], path[1], start_angle);
        T_0 = shrink_tangent(T_0, path[0], k, ps);
        // std::cout << "T_0 " << T_0.x() << " " << T_0.y() << std::endl;
        tangents.push_back(T_0);

        Vector2f T_e = calc_end_tangent(path[path.size() - 2], path[path.size() - 1]);
        T_e = shrink_tangent(T_e, path[path.size() - 1], k, ps);
        // std::cout << "T_e " << T_e.x() << " " << T_e.y() << std::endl;
        for (std::size_t i = 1; i < path.size() - 1; ++i) {
            // std::vector<Vector2f> ctrl_pts;
            // std::vector<Vector2f> ctrl_pts2;

            const Vector2f W_0 = path[i-1];
            const Vector2f W_1 = path[i];
            const Vector2f W_2 = path[i+1];

            // const Vector2f T_0 = calc_start_tangent(W_0, W_1, start_angle);
            Vector2f T_1 = calc_tangent(W_0, W_1, W_2);
            // std::cout << "T_1 " << T_1.x() << " " << T_1.y() << std::endl;
            // const Vector2f T_2 = calc_end_tangent(W_1, W_2);


            // T_0 = shrink_tangent(T_0, W_0, k, ps);
            T_1 = shrink_tangent(T_1, W_1, k, ps);
            tangents.push_back(T_1);
            // T_2 = shrink_tangent(T_2, W_2, k, ps);


            // ctrl_pts.push_back(W_0);
            // ctrl_pts.push_back(W_0 + T_0);
            // ctrl_pts.push_back(W_1 - T_1);
            // ctrl_pts.push_back(W_1);

            // ctrl_pts2.push_back(W_1);
            // ctrl_pts2.push_back(W_1 + T_1);
            // ctrl_pts2.push_back(W_2 - T_2);
            // ctrl_pts2.push_back(W_2);

            // bezier_spline bs = bezier_spline::bezier_curve(ctrl_pts, 0.0001);
            // bezier_spline bs2 = bezier_spline::bezier_curve(ctrl_pts2, 0.0001);
            // parts.push_back(bs);
            // parts.push_back(bs2);
        }
        tangents.push_back(T_e);

        std::vector<bezier_spline> parts;
        parts.reserve(path.size());
        for (std::size_t i = 0; i < path.size() - 1; ++i) {
            const Vector2f W_0 = path[i];
            const Vector2f W_1 = path[i+1];

            std::vector<Vector2f> ctrl_pts;

            ctrl_pts.push_back(W_0);
            ctrl_pts.push_back(W_0 + tangents[i]);
            ctrl_pts.push_back(W_1 - tangents[i+1]);
            ctrl_pts.push_back(W_1);
            // std::cout << "T[0] " << tangents[i].x() << " " << tangents[i].y() << std::endl;
            // std::cout << "T[1] " << tangents[i+1].x() << " " << tangents[i+1].y() << std::endl;
            Vector2f tmp = W_1 - tangents[i+1];
            // std::cout << tmp.x() << " " << tmp.y() << std::endl;
            // std::cout << "W_0 " << W_0.x() << " " << T_e.y() << std::endl;
            // std::cout << "W_1 " << ctrl_pts[1].x() << " " << ctrl_pts[1].y() << std::endl;
            // std::cout << "W_2 " << ctrl_pts[2].x() << " " << ctrl_pts[2].y() << std::endl;
            // std::cout << "W_3 " << W_1.x() << " " << W_1.y() << std::endl;
            bezier_spline bs = bezier_spline::bezier_curve(ctrl_pts, 0.0001);
            parts.push_back(bs);
        }
        return join_splines(parts);
    }

    inline bezier_spline bezier_spline::bezier_curve(const std::vector<Vector2f>& ctrl_pts, const std::vector<float>& positions) {
        const VectorXf mapped = VectorXf::Map(&positions[0], positions.size());
        return bezier_curve(ctrl_pts, mapped);
    }

    bezier_spline bezier_spline::bezier_curve(const std::vector<Vector2f>& ctrl_pts, const VectorXf& positions) {
        SC_ASSERT(ctrl_pts.size() >= 2, "ctrl_pts must have at least 2 points");
        SC_ASSERT(positions.minCoeff() >= 0, "positions must be between [0, 1]");
        SC_ASSERT(positions.minCoeff() <= 1, "positions must be between [0, 1]");

        const int degree = ctrl_pts.size() - 1;


        // TODO: figure out if we can avoid constructing a new FFT object
        FFT<float> fft;

        Matrix<std::complex<float>, Dynamic, 2> U = Matrix<float, Dynamic, 2>::Zero(ctrl_pts.size(), 2);

        for (std::size_t i = 0; i <= degree; ++i) {
            U(i, 0) = ctrl_pts[i].x();
            U(i, 1) = ctrl_pts[i].y();
        }

        Matrix<std::complex<float>, Dynamic, 2> Q = Matrix<float, Dynamic, 2>::Zero(ctrl_pts.size(), 2);

        Eigen::VectorXcf tmp_fft(ctrl_pts.size());
        fft.inv(tmp_fft, U.col(0));
        Q.col(0) = tmp_fft;

        fft.inv(tmp_fft, U.col(1));
        Q.col(1) = tmp_fft;

        return bezier_curve(ctrl_pts, positions, Q);
    }

    bezier_spline bezier_spline::bezier_curve(const std::vector<Vector2f>& ctrl_pts, const VectorXf& positions, Matrix<std::complex<float>, Dynamic, 2> Q) {
        SC_ASSERT(ctrl_pts.size() >= 2, "ctrl_pts must have at least 2 points");
        SC_ASSERT(positions.minCoeff() >= 0, "positions must be between [0, 1]");
        SC_ASSERT(positions.minCoeff() <= 1, "positions must be between [0, 1]");

        const int degree = ctrl_pts.size() - 1;

        Matrix<float, Dynamic, 2> B;
        B.setZero(positions.rows(), 2);

        const std::vector<std::complex<float>> omegas = bezier_spline::omega_table(degree);

        // TODO: apply the optimizations described in section 3 of
        // "Efficient computation of Bezier curves from their Bernstein-Fourier representation"

        for (int i = 0; i < positions.rows(); ++i) {
            const float s = positions(i);
            for (int k = 0; k <= degree; ++k) {
                std::complex<float> tmp = std::pow((1.0f+0if) + s*(omegas[k] - (1.0f+0if)), degree);
                B(i, 0) += (Q(k, 0) * tmp).real();
                B(i, 1) += (Q(k, 1) * tmp).real();
            }
        }

        bezier_spline bs = bezier_spline({ctrl_pts}, B, {positions});
        bs.Q_cache = {Q};
        return bs;
    }

    bezier_spline bezier_spline::bezier_curve(std::vector<Vector2f>& ctrl_pts, const float precision) {
        SC_ASSERT(precision < 1 && precision > 0, "spline percision must be in (0, 1)");

        const int n = static_cast<int>(std::round(1.0f / precision));

        SC_ASSERT(std::abs(n - (1.0f / precision)) < 0.0001,
                    "1/percision must be (close) to an integer, for arbitrary position values use the other bezier_curve method");

        std::vector<float> positions(n+1);
        for (std::size_t i = 0; i <= n; ++i) {
            positions[i] = std::max(0.0f, std::min(i * precision, 1.0f));
        }
        return bezier_curve(ctrl_pts, positions);
    }


    // std::tuple<float, std::vector<std::vector<float>>> bezier_spline::arclength(const float precision=0.01) const {
    arclength_data bezier_spline::arclength(const float precision=0.01) const {
        constexpr std::array<double, 32> weights = {
            0.0965400885147278005667648300635757947368606312355700687323182099577497758679466512968173871061464644599963197828969869820251559172455698832434930732077927850876632725829187045819145660710266452161095406358159608874152584850413283587913891015545638518881205600825069096855488296437485836866,
            0.0965400885147278005667648300635757947368606312355700687323182099577497758679466512968173871061464644599963197828969869820251559172455698832434930732077927850876632725829187045819145660710266452161095406358159608874152584850413283587913891015545638518881205600825069096855488296437485836866,
            0.0956387200792748594190820022041311005948905081620055509529898509437067444366006256133614167190847508238474888230077112990752876436158047205555474265705582078453283640212465537132165041268773645168746774530146140911679782502276289938840330631903789120176765314495900053061764438990021439069,
            0.0956387200792748594190820022041311005948905081620055509529898509437067444366006256133614167190847508238474888230077112990752876436158047205555474265705582078453283640212465537132165041268773645168746774530146140911679782502276289938840330631903789120176765314495900053061764438990021439069,
            0.0938443990808045656391802376681172600361000757462364500506275696355695118623098075097804207682530277555307864917078828352419853248607668520631751470962234105835015158485760721979732297206950719908744248285672032436598213262204039212897239890934116841559005147755270269705682414708355646603,
            0.0938443990808045656391802376681172600361000757462364500506275696355695118623098075097804207682530277555307864917078828352419853248607668520631751470962234105835015158485760721979732297206950719908744248285672032436598213262204039212897239890934116841559005147755270269705682414708355646603,
            0.0911738786957638847128685771116370625448614132753900053231278739777031520613017513597426417145878622654027367650308019870251963114683369110451524174258161390823876554910693202594383388549640738095422966058367070348943662290656339592299608558384147559830707904449930677260444604329157917977,
            0.0911738786957638847128685771116370625448614132753900053231278739777031520613017513597426417145878622654027367650308019870251963114683369110451524174258161390823876554910693202594383388549640738095422966058367070348943662290656339592299608558384147559830707904449930677260444604329157917977,
            0.0876520930044038111427714627518022875484497217017572223192228034747061150211380239263021665771581379364685191248848158059408000065275041643745927401342920150588893827207354226012701872322225514682178439577327346929209121046816487338309068375228210705166692551938339727096609740531893725675,
            0.0876520930044038111427714627518022875484497217017572223192228034747061150211380239263021665771581379364685191248848158059408000065275041643745927401342920150588893827207354226012701872322225514682178439577327346929209121046816487338309068375228210705166692551938339727096609740531893725675,
            0.0833119242269467552221990746043486115387468839428344598401864047287594069244380966536255650452315042012372905572506028852130723585016898197140339352228963465326746426938359210160503509807644396182380868089959855742801355208471205261406307895519604387550841954817025499019984032594036141439,
            0.0833119242269467552221990746043486115387468839428344598401864047287594069244380966536255650452315042012372905572506028852130723585016898197140339352228963465326746426938359210160503509807644396182380868089959855742801355208471205261406307895519604387550841954817025499019984032594036141439,
            0.078193895787070306471740918828306671039786798482159190307481553869493700115196435401943819761440851294456424770323467367505109006517482028994114252939401250416132320553639542341400437522236191275346323130525969269563653003188829786549728825182082678498917784036375053244425839341945385297,
            0.078193895787070306471740918828306671039786798482159190307481553869493700115196435401943819761440851294456424770323467367505109006517482028994114252939401250416132320553639542341400437522236191275346323130525969269563653003188829786549728825182082678498917784036375053244425839341945385297,
            0.0723457941088485062253993564784877916043369833018248707397632823511765345816800402874475958591657429073027694582930574378890633404841054620298756279975430795706338162404545590689277985270140590721779502609564199074051863640176937117952488466002340085264819537808079947788437998042296495822,
            0.0723457941088485062253993564784877916043369833018248707397632823511765345816800402874475958591657429073027694582930574378890633404841054620298756279975430795706338162404545590689277985270140590721779502609564199074051863640176937117952488466002340085264819537808079947788437998042296495822,
            0.0658222227763618468376500637069387728775364473732465153710916696852412442018627316280044447764609054151761388378861151807154113495715653711918644796313239555117970398473141615070299152284100887258072240524028885129828725430021172354299810423059697133688823072212214503334259555369485963074,
            0.0658222227763618468376500637069387728775364473732465153710916696852412442018627316280044447764609054151761388378861151807154113495715653711918644796313239555117970398473141615070299152284100887258072240524028885129828725430021172354299810423059697133688823072212214503334259555369485963074,
            0.0586840934785355471452836373001708867501204674575467587150032786132877518019090643743123653437052116901895704813134467814193905269714480573030647540887991405215103758723074481312705449946311993670933802369300463315125015975216910705047901943865293781921122370996257470349807212516159332678,
            0.0586840934785355471452836373001708867501204674575467587150032786132877518019090643743123653437052116901895704813134467814193905269714480573030647540887991405215103758723074481312705449946311993670933802369300463315125015975216910705047901943865293781921122370996257470349807212516159332678,
            0.0509980592623761761961632446895216952601847767397628437069071236525030510385137821267442193868358292147899714519363571211100873456269865150186456681043804358654826791768545393024953758025593924464295555854744882720755747096079325496814455853004350452095212995888025282619932613606999567133,
            0.0509980592623761761961632446895216952601847767397628437069071236525030510385137821267442193868358292147899714519363571211100873456269865150186456681043804358654826791768545393024953758025593924464295555854744882720755747096079325496814455853004350452095212995888025282619932613606999567133,
            0.0428358980222266806568786466061255284928108575989407395620219408911043916962572261359138025961596979511472539467367407419206021900868371610612953162236233351132214438513203223655531564777278515080476421262443325932320214191168239648611793958596884827086182431203349730049744697408543115307,
            0.0428358980222266806568786466061255284928108575989407395620219408911043916962572261359138025961596979511472539467367407419206021900868371610612953162236233351132214438513203223655531564777278515080476421262443325932320214191168239648611793958596884827086182431203349730049744697408543115307,
            0.0342738629130214331026877322523727069948402029116274337814057454192310522168984446294442724624445760666244242305266023810860790282088335398182296698622433517061843276344829146573593201201081743714879684153735672789104567624853712011151505225193933019375481618760594889854480408562043658635,
            0.0342738629130214331026877322523727069948402029116274337814057454192310522168984446294442724624445760666244242305266023810860790282088335398182296698622433517061843276344829146573593201201081743714879684153735672789104567624853712011151505225193933019375481618760594889854480408562043658635,
            0.0253920653092620594557525897892240292875540475469487209362512822192154788532376645960457016338988332029324531233401833547954942765653767672102838323550828207273795044402516181251040411735351747299230615776597356956641506445501689924551185923348003766988424170511157069264716719906995309826,
            0.0253920653092620594557525897892240292875540475469487209362512822192154788532376645960457016338988332029324531233401833547954942765653767672102838323550828207273795044402516181251040411735351747299230615776597356956641506445501689924551185923348003766988424170511157069264716719906995309826,
            0.0162743947309056706051705622063866181795429637952095664295931749613369651752917857651844425586692833071042366002861684552859449530958901379260437604156888337987656773068694383447504913457771896770689760342192010638946676879735404121702279005140285599424477022083127753774756520463311689155,
            0.0162743947309056706051705622063866181795429637952095664295931749613369651752917857651844425586692833071042366002861684552859449530958901379260437604156888337987656773068694383447504913457771896770689760342192010638946676879735404121702279005140285599424477022083127753774756520463311689155,
            0.0070186100094700966004070637388531825133772207289396032320082356192151241454178686953297376907573215077936155545790593837513204206518026084505878987243348925784479817181234617862457418214505322067610482902501455504204433524520665822704844582452877416001060465891907497519632353148380799619,
            0.0070186100094700966004070637388531825133772207289396032320082356192151241454178686953297376907573215077936155545790593837513204206518026084505878987243348925784479817181234617862457418214505322067610482902501455504204433524520665822704844582452877416001060465891907497519632353148380799619};

        constexpr std::array<double, 32> abscissa = {
            -0.0483076656877383162348125704405021636908472517308488971677937345463685926042778777794060365911173780988289503411375793689757446357461295741679964108035347980667582792392651327368009453047606446744575790523465655622949909588624860214137051585425884056992683442137333250625173849291299678673,
            0.0483076656877383162348125704405021636908472517308488971677937345463685926042778777794060365911173780988289503411375793689757446357461295741679964108035347980667582792392651327368009453047606446744575790523465655622949909588624860214137051585425884056992683442137333250625173849291299678673,
            -0.1444719615827964934851863735988106522038459913156355521379528938242184438164519731102406769974924713989580220758441301598578946580142268413547299935841673092513202403499286272686350814272974392746706128556678811982653393383080797337231702069432462445053984587997153683967433095128570624414,
            0.1444719615827964934851863735988106522038459913156355521379528938242184438164519731102406769974924713989580220758441301598578946580142268413547299935841673092513202403499286272686350814272974392746706128556678811982653393383080797337231702069432462445053984587997153683967433095128570624414,
            -0.2392873622521370745446032091655015206088554219602530155470960995597029133039943915553593695844147813728958071901224632260145752503694970545640339873418480550362677768010887468668377893757173424222709744116861683634989914911762187599464033126988486345234374380695224452457957624756811128321,
            0.2392873622521370745446032091655015206088554219602530155470960995597029133039943915553593695844147813728958071901224632260145752503694970545640339873418480550362677768010887468668377893757173424222709744116861683634989914911762187599464033126988486345234374380695224452457957624756811128321,
            -0.3318686022821276497799168057301879961957751368050598360182296306285376829657438169809731852312743263005943551508559377834274303920771100489026913715847854727626540340157368609696698131829681988642689780208633461925468064919389286805624602715005948661328152252049795463242055567997437182143,
            0.3318686022821276497799168057301879961957751368050598360182296306285376829657438169809731852312743263005943551508559377834274303920771100489026913715847854727626540340157368609696698131829681988642689780208633461925468064919389286805624602715005948661328152252049795463242055567997437182143,
            -0.4213512761306353453641194361724264783358772886324433305416613404557190462549837315607633055675740638739884093394574651160978879545562247406839036854173715776910866941643197988581928900702286425821151586000969947406313405310082646561917980302543820974679501841964453794193724645925031841919,
            0.4213512761306353453641194361724264783358772886324433305416613404557190462549837315607633055675740638739884093394574651160978879545562247406839036854173715776910866941643197988581928900702286425821151586000969947406313405310082646561917980302543820974679501841964453794193724645925031841919,
            -0.5068999089322293900237474743778212301802836995994354639743662809707712640478764442266190213124522047999876916596854537447047905434649918210338296049592120273725464263651562560829050004258268002241145951271730860506703690843719936432852920782304931272053564539127514959875734718036950073563,
            0.5068999089322293900237474743778212301802836995994354639743662809707712640478764442266190213124522047999876916596854537447047905434649918210338296049592120273725464263651562560829050004258268002241145951271730860506703690843719936432852920782304931272053564539127514959875734718036950073563,
            -0.5877157572407623290407454764018268584509401154544205727031788473129228586684474311408145102018661764979429510790747919023774933113319119601088669936958908618326367715806216053155906936017362413244183150445492317940727345571648726363597097311647731726438279098059670236086983675374932643925,
            0.5877157572407623290407454764018268584509401154544205727031788473129228586684474311408145102018661764979429510790747919023774933113319119601088669936958908618326367715806216053155906936017362413244183150445492317940727345571648726363597097311647731726438279098059670236086983675374932643925,
            -0.6630442669302152009751151686632383689770222859605053010170834964924461749232229404368981536611965356686820332804126742949900731319113817214392193185613161549689934301410316417342588149871686184296988807305719690974644891055567340650986465615021143958920599684258616066247948224049997371166,
            0.6630442669302152009751151686632383689770222859605053010170834964924461749232229404368981536611965356686820332804126742949900731319113817214392193185613161549689934301410316417342588149871686184296988807305719690974644891055567340650986465615021143958920599684258616066247948224049997371166,
            -0.732182118740289680387426665091267146630270483506629100821139573270385253587797727611292298988652560055905228466313310601075333829094630570926240639601009902567982815376254840388565733846030450161774620971196087756484387383432502715118096615117242484073636640563609696801484680439912327302,
            0.732182118740289680387426665091267146630270483506629100821139573270385253587797727611292298988652560055905228466313310601075333829094630570926240639601009902567982815376254840388565733846030450161774620971196087756484387383432502715118096615117242484073636640563609696801484680439912327302,
            -0.7944837959679424069630972989704289020954794016388354532507582449720593922816426654241878967890821228397041480126630294067578180914548706957761322921470535094589673860419616615738928385807346185892317514562489971543238450942224396667500582904031225063621511429185567036727089257387570529468,
            0.7944837959679424069630972989704289020954794016388354532507582449720593922816426654241878967890821228397041480126630294067578180914548706957761322921470535094589673860419616615738928385807346185892317514562489971543238450942224396667500582904031225063621511429185567036727089257387570529468,
            -0.849367613732569970133693004967742538954886793049759233100219598613724656141562558741881463752754991143937635778596582088915769685796612254240615386941355933272723068952531445772190363422003834495043219316062885999846179078139659341918527603834809670576387535564876596379488780285979062125,
            0.849367613732569970133693004967742538954886793049759233100219598613724656141562558741881463752754991143937635778596582088915769685796612254240615386941355933272723068952531445772190363422003834495043219316062885999846179078139659341918527603834809670576387535564876596379488780285979062125,
            -0.8963211557660521239653072437192122684789964967957595765636154129650249794910409173494503783167666654202705333374285522819507600044591355080910768854012859468015827508424619812224062460791781333400979810176198916239783226706506012473250929962326307746466256167673927887144428859779028909399,
            0.8963211557660521239653072437192122684789964967957595765636154129650249794910409173494503783167666654202705333374285522819507600044591355080910768854012859468015827508424619812224062460791781333400979810176198916239783226706506012473250929962326307746466256167673927887144428859779028909399,
            -0.9349060759377396891709191348354093255286714322828372184584037398118161947182932855418880831417927728359606280450921427988850058691931014887248988124656348299653052688344696135840215712191162135178273756415771123010111796122671724143565383396162107206772781551029308751511942924942333859805,
            0.9349060759377396891709191348354093255286714322828372184584037398118161947182932855418880831417927728359606280450921427988850058691931014887248988124656348299653052688344696135840215712191162135178273756415771123010111796122671724143565383396162107206772781551029308751511942924942333859805,
            -0.9647622555875064307738119281182749603888952204430187193220113218370995254867038008243801877562227002840740910741483519987441236283464394249183812395373150090695515823078220949436846111682404866338388944248976976566275875721000356873959697266702651250019105084704924793016185368873243713355,
            0.9647622555875064307738119281182749603888952204430187193220113218370995254867038008243801877562227002840740910741483519987441236283464394249183812395373150090695515823078220949436846111682404866338388944248976976566275875721000356873959697266702651250019105084704924793016185368873243713355,
            -0.9856115115452683354001750446309019786323957143358063182107821705820305847193755946663846485510970266115353839862364606643634021712823093784875255943834038377710426488328772047833289470320023596895438028281274741367781028592272459887917924171204666683239464005128153533797603112851826904814,
            0.9856115115452683354001750446309019786323957143358063182107821705820305847193755946663846485510970266115353839862364606643634021712823093784875255943834038377710426488328772047833289470320023596895438028281274741367781028592272459887917924171204666683239464005128153533797603112851826904814,
            -0.9972638618494815635449811286650407271385376637294611593011185457862359083917418520130456693085426416474280482200936551645510686196373231416035137741332968299789863385253514914078766236061488136738023162574655835389902337937054326098485227311719825229066712510246574949376367552421728646398,
            0.9972638618494815635449811286650407271385376637294611593011185457862359083917418520130456693085426416474280482200936551645510686196373231416035137741332968299789863385253514914078766236061488136738023162574655835389902337937054326098485227311719825229066712510246574949376367552421728646398};

        static_assert(weights.size() == abscissa.size());

        float total_arclen = 0;
        std::vector<VectorXf> arclens(n_segments());
        std::vector<VectorXf> arclen_positions(n_segments());
        for (int s = 0; s < n_segments(); ++s) {
            SC_ASSERT(precision < 1 && precision > 0, "spline percision must be in (0, 1)");

            const int n = static_cast<int>(std::round(1.0f / precision));

            SC_ASSERT(std::abs(n - (1.0f / precision)) < 0.0001,
                        "1/percision must be (close) to an integer, for arbitrary position values use the other bezier_curve method");
            
            std::vector<float> pos(n+1);
            for (std::size_t i = 0; i <= n; ++i) {
                pos[i] = std::min(i * precision, 1.0f);
            }

            VectorXf deriv_pos = VectorXf::Zero(weights.size() * (pos.size()-1));
            for (int k = 0; k < pos.size()-1; ++k) {
                const float b = pos[k+1];
                const float a = pos[k];

                std::vector<float> ab_pos(weights.size());
                for (int i = 0; i < weights.size(); ++i) {
                    deriv_pos((k*weights.size())+i) = ((b-a)/2)*abscissa[i] + ((b+a)/2);
                }
            }

            bezier_spline deriv = (bezier_spline({ctrl_pts[s]}, Matrix<float, Dynamic, 2>(), {deriv_pos})).hodograph();

            VectorXf arclen_seg = VectorXf::Zero(pos.size());
            VectorXf arclen_seg_pos = VectorXf::Zero(pos.size());
            for (int k = 0; k < pos.size()-1; ++k) {
                const float b = pos[k+1];
                const float a = pos[k];

                for (int i = 0; i < weights.size(); ++i) {
                    const int ind = (k*weights.size())+i;
                    arclen_seg(k+1) += weights[i] * std::sqrt((deriv.pts(ind, 0) * deriv.pts(ind, 0)) + (deriv.pts(ind, 1) * deriv.pts(ind, 1)));
                }
                arclen_seg[k+1] *= ((b-a)/2.0f);
                arclen_seg_pos[k+1] = b;
                total_arclen += arclen_seg[k+1];
            }

            float sum = 0;
            for (float& arc : arclen_seg) {
                sum += arc;
                arc = sum;
            }
            arclens[s] = arclen_seg;
            arclen_positions[s] = arclen_seg_pos;
        }
        arclength_data ad;
        ad.arclength = total_arclen;
        ad.segments = arclens;
        ad.positions = arclen_positions;

        return ad;
    }

    bezier_spline bezier_spline::resample(VectorXf& profile_pos, arclength_data ad, bool nudge_positions=false) const {

        // This method of nudging only works well for isolated cases of weird values
        // TODO: consider a more robust way to handle small numerical errors
        if (nudge_positions) {
            profile_pos(0) = 0;
            profile_pos(profile_pos.rows()-1) = ad.arclength;
            for (int i = 1; i < profile_pos.rows()-1; ++i) {
                if (profile_pos(i) < profile_pos(i-1) || profile_pos(i) > profile_pos(i+1)) {
                    profile_pos(i) = (profile_pos(i-1) + profile_pos(i+1)) / 2;
                }
                if (profile_pos(i) < 0) profile_pos(i) = 0;
                if (profile_pos(i) > ad.arclength) profile_pos(i) = ad.arclength;
            }
        }

        SC_ASSERT(profile_pos.rows() > 0, "The vector of positions to be sampled must not be empty");
        // std::cout << "pp mc " << profile_pos.maxCoeff() << std::endl;
        // std::cout << "pp mc2 " << profile_pos(profile_pos.rows() - 2) << std::endl;
        // std::cout << "arklen " << ad.arclength << std::endl;
        SC_ASSERT(profile_pos.maxCoeff() <= ad.arclength, "The profile can not go beyond the arclength of the spline");
        SC_ASSERT(profile_pos.minCoeff() >= 0, "The profile can not go beyond the arclength of the spline");

        std::vector<bezier_spline> curves(positions.size());

        // std::cout << "ppr " << profile_pos.rows() << std::endl;
        int j = 0;
        float offset = 0;
        std::size_t block_size_sum = 0;
        for (int i = 0; i < positions.size(); ++i) {
            const VectorXf seg = ad.segments[i];
            const float last = seg(seg.rows()-1);
            // std::cout << "first " << seg(0) << std::endl;
            // std::cout << "last " << last << std::endl;
            // std::cout << "offset " << offset << std::endl;
            int start = j;
            for (; j < profile_pos.rows() && (profile_pos(j) - offset <= last); ++j) {} // cursed
            // std::cout << "ppj - off " << (profile_pos(j) - offset) << std::endl;
            // std::cout << "diff " << (j - start) << std::endl;
            if (i + 1 == positions.size() && i == 0) j = profile_pos.rows();
            else if (i + 1 == positions.size()) j = profile_pos.rows() - 1;
            j -= 1;
            // std::cout << "ppj - off2 " << (profile_pos(j) - offset) << std::endl;
            // std::cout << "j " << j << std::endl;
            offset = profile_pos(j);

            const VectorXf block = profile_pos.block(start, 0, ((j+1)-start), 1).array() - profile_pos(start);
            block_size_sum += block.rows();
            std::cout << "block size " << block.rows() << std::endl;
            std::cout << "block end " << block(block.rows() - 1) << std::endl;

            // const int degree = std::max(seg.rows()/2, 2L);
            const int degree = std::min(10L, seg.rows()); // TODO: figure out a better heuristic for polynomial degree


            // polynomial from segment arclength to [0,1]
            const chebpoly poly = chebfit(seg, ad.positions[i], degree);
            VectorXf positions_fixed = chebeval(block, poly, degree);

            // std::cout << "positions_fixed " << positions_fixed(positions_fixed.rows() - 1) << std::endl;

            // std::cout << "degree " << degree << std::endl;
            // std::cout << "ad.positions[0] " << ad.positions[i][0] << std::endl;
            // std::cout << "ad.positions[-1] " << ad.positions[i][ad.positions[i].rows() - 1] << std::endl;
            // std::cout << "ad.positions[max] " << ad.positions[i].maxCoeff() << std::endl;

            // std::cout << positions_fixed.minCoeff() << std::endl;
            // std::cout << positions_fixed.rows() << std::endl;
            for (std::size_t i = 0; i < positions_fixed.rows(); ++i) {
                if (positions_fixed(i) < 0) positions_fixed(i) = 0;
                if (positions_fixed(i) > 1) positions_fixed(i) = 1;
                // std::cout << positions_fixed(i) << std::endl;
            }
            curves[i] = bezier_curve(ctrl_pts[i], positions_fixed, Q_cache[i]);

            // std::cout << "curve size " << curves[i].n_pts() << std::endl;
            // std::cout << "curve end " << curves[i].arclength().arclength << std::endl;
            // std::cout << "curve end pos " << curves[i].positions[0](curves[i].positions[0].rows() - 1) << std::endl;
            // std::cout << "-------------" << std::endl;
        }

        // float sum = 0;
        // for (auto& curve : curves) {
        //     sum += curve.arclength().arclength;
        //     std::cout << curve.arclength().arclength << std::endl;
        // }
        // std::cout << "sum " << sum << std::endl;

        if (block_size_sum > profile_pos.rows()) {
            // remove pts and positions
            const int diff = (block_size_sum - profile_pos.rows());
            auto&& last_curve = curves[curves.size() - 1];
            // std::cout << "die here> " << std::endl;
            last_curve.pts = last_curve.pts.block(0, 0, last_curve.pts.rows() - diff, 2);
            auto&& last_curve_last_pos = last_curve.positions[last_curve.positions.size() - 1];
            // std::cout << "die here> " << std::endl;
            last_curve_last_pos = last_curve_last_pos.block(0, 0, last_curve_last_pos.rows() - diff, 1);
            // std::cout << "die here> " << std::endl;
        }

        bezier_spline fixed_spline = join_splines(curves);
        std::cout << "fixed spline " << fixed_spline.n_pts() << " " << profile_pos.rows() << std::endl;
        std::cout << block_size_sum << std::endl;
        std::cout << profile_pos(profile_pos.rows() - 1) << std::endl;
        std::cout << fixed_spline.arclength().arclength << std::endl;
        SC_ASSERT(fixed_spline.n_pts() == profile_pos.rows(), "fixed_spline.n_pts() == profile_pos.rows()");
        return fixed_spline;
    }

    // std::vector<float> num_diff(std::vector<float> f) {
    //     SC_ASSERT(f.size() >= 3, "");

    //     std::vector<float> d(f.size());
    //     for (std::size_t i = 1; i < f.size() - 1; ++i) {

    //         d[i] =
    //     }
    // }

    std::vector<float> bezier_spline::curvature() const {
        bezier_spline d = hodograph();
        bezier_spline dd = d.hodograph();

        std::vector<float> res(pts.rows());
        for (std::size_t i = 0; i < pts.rows(); ++i) {
            const Vector3f d_pt(d.pts(i, 0), d.pts(i, 1), 0);
            const Vector3f dd_pt(dd.pts(i, 0), dd.pts(i, 1), 0);

            // std::cout << "point 1: " << d_pt.x() << " " << d_pt.y() << std::endl;
            // std::cout << "point 2: " << dd_pt.x() << " " << dd_pt.y() << std::endl;

            const float k = (d_pt.x() * dd_pt.y() - d_pt.y() * dd_pt.x()) / std::pow(d_pt.x() * d_pt.x() + d_pt.y() * d_pt.y(), 1.5);
            res[i] = k;

            // const Vector3f pt3 = d_pt.cross(dd_pt);
            // std::cout << "point 3: " << pt3.x() << " " << pt3.y() << " " << pt3.z() << std::endl;

            // res[i] = d_pt.cross(dd_pt).norm() / std::pow(d_pt.norm(), 3);
        }

        return res;
    }

    bezier_spline bezier_spline::hodograph() const {
        std::vector<bezier_spline> curves(ctrl_pts.size());

        for (std::size_t i = 0; i < ctrl_pts.size(); ++i) {
            std::vector<Vector2f> d_cps(ctrl_pts[i].size()-1);
            for (std::size_t j = 0; j < ctrl_pts[i].size()-1; ++j) {
                d_cps[j] = degree() * (ctrl_pts[i][j+1] - ctrl_pts[i][j]);
            }
            curves[i] = bezier_curve(d_cps, positions[i]);
        }

        return join_splines(curves);
    }

    std::vector<float> bezier_spline::angular_velocity(const velocity_profile& vel_prof) const {

        const std::vector<float> curv = curvature();
        auto&& vels = vel_prof.vel[0];

        SC_ASSERT(curv.size() == vel_prof.vel[0].size(), "curvature and velocity vectors must be the same size");

        std::vector<float> angular_velocities(curv.size());
        for (std::size_t i = 0; i < angular_velocities.size(); ++i) {
            angular_velocities[i] = vels[i] * curv[i];
        }
        return angular_velocities;
    }

    std::vector<float> bezier_spline::angular_velocity2(const velocity_profile& vel_prof) const {
        bezier_spline d = hodograph();

        SC_ASSERT(pts.rows() == vel_prof.vel[0].size(), "");

        std::vector<float> res(pts.rows());

        std::vector<float> rads(pts.rows());
        for (std::size_t i = 0; i < pts.rows(); ++i) {
            Vector3f d_pt(d.pts(i, 0), d.pts(i, 1), 0);
            d_pt.normalize();
            const float rad = std::atan2(d_pt.y(), d_pt.x());
            std::cout << "rad " << rad << std::endl;
            rads[i] = rad;
        }

        res[0] = 0;
        res[res.size() - 1] = 0;

        auto&& times = vel_prof.time;
        for (std::size_t i = 1; i < res.size() - 1; ++i) {
            res[i] = (rads[i+1] - rads[i]) / (times[i+1] - times[i]);
        }

        return res;
    }

    inline std::vector<std::complex<float>> bezier_spline::omega_table(const int degree) {
        using std::complex;

        std::vector<std::complex<float>> omegas(degree+1);
        omegas[0] = 1.0f+0if;
        for (std::size_t i = 1; i <= degree; ++i) {
            omegas[i] = pow(exp((complex<float>(std::numbers::pi) * -2if) / complex<float>(degree+1)), i);
        }

        return omegas;
    }


    chebpoly chebfit(const VectorXf& x, const VectorXf& y, const int degree) {
        SC_ASSERT(degree >= 1, "degree must be a positive integer");
        SC_ASSERT(x.rows() == y.rows(), "x and y must have the same number of rows");

        const int n = degree;
        const int m = x.rows();
        const float xmax = x.maxCoeff();
        const float xmin = x.minCoeff();

        SC_ASSERT(std::abs(xmax - xmin) > 0.00001, "Error: vector x should not have all equal values");
        SC_ASSERT(degree >= 1, "degree must be >= 1");

        const VectorXf x_norm = ((2*x).array() - (xmax + xmin)) / (xmax - xmin);

        MatrixXf T = MatrixXf::Zero(m, n);
        T.col(0) = VectorXf::Ones(m);
        if (n >= 1) T.col(1) = x_norm;

        for (int j = 2; j < n; ++j) {
            T.col(j) = (2*x_norm).array() * T.col(j-1).array() - T.col(j-2).array();
        }

        // ColPivHouseholderQR<MatrixXf> T_Qr = T.colPivHouseholderQr();
        // SC_ASSERT(T_Qr.rank() == degree, "");

        HouseholderQR<MatrixXf> T_Qr = T.householderQr();
        SC_ASSERT(T.colPivHouseholderQr().rank() == degree, "T.colPivHouseholderQr().rank() == degree");

        return chebpoly(T_Qr.solve(y), xmin, xmax);
    }

    VectorXf chebeval(const VectorXf& x, const chebpoly& b, const int degree) {
        SC_ASSERT(degree >= 1, "degree must be a positive integer");

        const int n = degree;
        const int m = x.rows();
        const float xmax = b.xmax;
        const float xmin = b.xmin;

        SC_ASSERT(std::abs(xmax - xmin) > 0.00001, "Error: vector x should not have all equal values");
        SC_ASSERT(degree >= 1, "degree must be >= 1");

        const VectorXf x_norm = ((2*x).array() - (xmax + xmin)) / (xmax - xmin);

        VectorXf y = VectorXf::Zero(m);

        MatrixXf T = MatrixXf::Zero(m, n);
        T.col(0) = VectorXf::Ones(m);
        y += b.coeffs(0) * T.col(0);

        if (n >= 1) {
            T.col(1) = x_norm;
            y += b.coeffs(1) * T.col(1);
        }

        for (int j = 2; j < n; ++j) {
            T.col(j) = (2*x_norm).array() * T.col(j-1).array() - T.col(j-2).array();
            y += b.coeffs(j) * T.col(j);
        }

        return y;
    }

    using toppra::value_type;
    using toppra::constraint::LinearJointVelocity;

    using vel_lim_func = std::function<std::tuple<toppra::Vector, toppra::Vector>(value_type time)>;

    class LinearJointVelocityVarying : public LinearJointVelocity {
        public:
        vel_lim_func calc_lim;
        LinearJointVelocityVarying(int nDof, vel_lim_func calc_lim) : LinearJointVelocity (-1*toppra::Vector::Ones(1), 1*toppra::Vector::Ones(1)) {
            this->calc_lim = calc_lim;
            computeVelocityLimits(0);
        }
        protected:
        void computeVelocityLimits(value_type time) {
            std::tie(m_lower, m_upper) = calc_lim(time);
        }
    };


    template <int N> requires (N >= 1)
    velocity_profile gen_vel_prof(const Vector<value_type, N>& pos_end,
                                  const Vector<value_type, N>& pos_start,
                                  const Vector<value_type, N>& vel_end,
                                  const Vector<value_type, N>& vel_start,
                                  const vel_lim_func& vel_lim,
                                  const Vector<value_type, N>& acc_min,
                                  const Vector<value_type, N>& acc_max,
                                  const float dt=0.02) {
        using namespace toppra;
        using namespace toppra::constraint;

        const int dof = pos_end.rows();

        LinearJointVelocityVarying vel_con(dof, vel_lim);

        toppra::LinearConstraintPtr ljv, lja;
        ljv = std::make_shared<LinearJointVelocityVarying>(vel_con);
        lja = std::make_shared<toppra::constraint::LinearJointAcceleration>(acc_min, acc_max);
        lja->discretizationType(toppra::DiscretizationType::Interpolation);
        toppra::LinearConstraintPtrs constraints{ljv, lja};

        toppra::Vectors positions = {pos_start, pos_end};

        toppra::Vectors velocities = {vel_start, vel_end};

        std::vector<toppra::value_type> steps;
        steps = std::vector<toppra::value_type>{0, 1};

        toppra::PiecewisePolyPath hermite = toppra::PiecewisePolyPath::CubicHermiteSpline(positions, velocities, steps);

        toppra::GeometricPathPtr path = std::make_shared<PiecewisePolyPath>(hermite);

        toppra::algorithm::TOPPRA algo(constraints, path);
        toppra::ReturnCode rc1 = algo.computePathParametrization(0, 0);

        SC_ASSERT(rc1 == toppra::ReturnCode::OK, "");

        toppra::ParametrizationData pd = algo.getParameterizationData();

        toppra::Vector gridpoints = pd.gridpoints;
        toppra::Vector vsquared = pd.parametrization;
        toppra::parametrizer::Spline spp(path, gridpoints, vsquared);

        Eigen::Matrix<toppra::value_type, 1, 2> interval = spp.pathInterval();

        int length = std::ceil((interval(1) - interval(0)) / dt);
        // std::cout << "length " << length << std::endl;
        toppra::Vector times = toppra::Vector::LinSpaced(length, interval(0), interval(1));

        toppra::Vectors path_pos = spp.eval(times, 0);
        toppra::Vectors path_vel = spp.eval(times, 1);
        toppra::Vectors path_acc = spp.eval(times, 2);

        std::vector<VectorXf> pos(dof);
        std::vector<VectorXf> vel(dof);
        std::vector<VectorXf> acc(dof);

        for (int j = 0; j < dof; ++j) {
            pos[j] = VectorXf::Zero(length);
            vel[j] = VectorXf::Zero(length);
            acc[j] = VectorXf::Zero(length);
        }

        // TODO: there is probably a better way to copy this data
        for (int i = 0; i < path_pos.size(); ++i) {
            for (int j = 0; j < dof; ++j) {
                pos[j](i) = path_pos[i](j);
                vel[j](i) = path_vel[i](j);
                acc[j](i) = path_acc[i](j);
            }
        }

        return velocity_profile(pos, vel, acc, times);
    }






    planning_space::planning_space(const bounding_rect& br) : bound_rect(br) {}

    point_set planning_space::sample_free(const int n) {
        // technically this should be a set, but the halton sequence is guaranteed to not repeat
        // so we can avoid element checks for a set
        point_set pts = {Vector2f(0,0)};
        pts.reserve(n);

        while (n > pts.size()) {
            // TODO: consider using different bases, or expose the bases to the user
            std::vector<float> x_test = halton(2, n - pts.size(), x_state);
            std::vector<float> y_test = halton(3, n - pts.size(), y_state);

            for (int i = 0; i < x_test.size(); ++i) {
                const Vector2f test((bound_rect.x_max-bound_rect.x_min)*(x_test[i])+bound_rect.x_min, (bound_rect.y_max-bound_rect.y_min)*(y_test[i])+bound_rect.y_min);

                bool add = true;
                for (auto& obstacle : obstacles) {
                    if (obstacle.contains(test)) add = false;
                }
                if (add) pts.insert(test);
            }
        }

        return pts;
    }

    inline float planning_space::cost(const Vector2f a, const Vector2f b) const {
        const std::tuple<Vector2f, Vector2f> line = {a, b};
        for (auto& obstacle : obstacles) {
            for (auto& oline : obstacle.lines) {
                if (std::get<0>(intersects(line, oline))) {
                    return std::numeric_limits<float>::max();
                }
            }
        }

        return pt_dist(a, b);
    }

    inline point_set planning_space::near(const Vector2f b, const point_set& nodes, const float dist) const {
        point_set nodes_out;
        nodes_out.reserve(nodes.size());
        for (auto& a : nodes) {
            if (pt_dist(a, b) <= std::pow(dist, 2) && a != b) {
                nodes_out.insert(a);
            }
        }
        return nodes_out;
    }

    std::optional<std::vector<Vector2f>> planning_space::fast_marching_trees(const Vector2f& x_init, const Vector2f& x_goal, const int n, const float rn) {
        // TODO; replace std::optional return type with std::variant + some error type
        point_set V_closed;
        point_set V_open = {x_init};
        point_set V_unvisited = sample_free(n);
        V_unvisited.insert(x_goal);

        Vector2f z = x_init;

        std::unordered_map<Vector2f, std::optional<float>, hash_vector2f> cost_map;
        cost_map.emplace(x_init, std::optional<float>{0});

        std::unordered_map<Vector2f, Vector2f, hash_vector2f> parent_map;
        parent_map.emplace(x_init, x_init);


        const float inf = std::numeric_limits<float>::max();

        while (z != x_goal) {
            point_set V_open_new;
            point_set X_near = near(z, V_unvisited, rn);
            for (const auto& x : X_near) {
                const point_set Y_near = near(x, V_open, rn);
                if (Y_near.size() == 0) continue;
                Vector2f y_min = *Y_near.begin();
                for (const auto& y : Y_near) {
                    const float cost_y_min = cost_map[y_min].value_or(inf) + cost(x, y_min);
                    const float cost_y = cost_map[y].value_or(inf) + cost(x, y);
                    if (cost_y < cost_y_min) y_min = y;
                }

                if (cost(x, y_min) != inf) {
                    parent_map.insert_or_assign(x, y_min);
                    V_open_new.insert(x);
                    V_unvisited.erase(x);
                    cost_map.insert_or_assign(x, cost_map[y_min].value_or(inf) + cost(x, y_min));
                }

            }

            V_open.merge(V_open_new);
            V_open.erase(z);
            V_closed.insert(z);

            if (V_open.size() == 0) {
                return std::nullopt;
            }

            z = *V_open.begin();
            for (const auto& y : V_open) {
                const float cost_z = cost_map[z].value_or(inf);
                const float cost_y = cost_map[y].value_or(inf);
                if (cost_y < cost_z) {
                    z = y;
                }
            }
        }

        std::vector<Vector2f> path;

        Vector2f p = z;
        while (p != x_init) {
            path.push_back(p);
            p = parent_map[p];
        }
        path.push_back(x_init);
        std::reverse(path.begin(), path.end());
        return std::optional<std::vector<Vector2f>>{path};
    }


    template <typename T>
    std::vector<T> format_vec_vecx(const std::vector<VectorXf>& prof) {
        std::vector<T> prof_ser;
        for (std::size_t i = 0; i < prof.size(); ++i) {
            std::vector<T> tmp = std::vector<T>(prof[i].data(), prof[i].data()+prof[i].size());
            for (std::size_t j = 0; j < tmp.size(); ++j) {
                prof_ser.push_back(tmp[i]);
            }
        }
        return prof_ser;
    }

    // quick and dirty serialization
    json serialize_path_to_json(const bezier_spline& spline, const velocity_profile& vel_prof, const arclength_data& arclens, const std::vector<float>& ang_vel) {
        json j = json::array();

        //j["position"] = format_vec_vecx<float>(vel_prof.pos);
        const std::vector<float> velocities = format_vec_vecx<float>(vel_prof.vel);
        const std::vector<float> accelerations = format_vec_vecx<float>(vel_prof.acc);
        const std::vector<float> times = std::vector<float>(vel_prof.time.data(), vel_prof.time.data() + vel_prof.time.size());

        auto&& pts = spline.pts;
        std::vector<float> pos_x = std::vector<float>(pts.col(0).data(), pts.col(0).data() + pts.rows());
        std::vector<float> pos_y = std::vector<float>(pts.col(1).data(), pts.col(1).data() + pts.rows());

        // const std::vector<std::vector<float>> segments = format_vec_vecx<float>(arclens.segments);
        // const std::vector<std::vector<float>> positions = format_vec_vecx<float>(arclens.positions);
        // j["arclength"] = { {"arclength", arclens.arclength}, {"segments", segments}, {"positions", positions} };

        for (std::size_t i = 0; i < pts.rows(); ++i) {
            json j2;
            j2["time"] = times[i];
            j2["velocity"] = velocities[i];
            j2["acceleration"] = accelerations[i];
            j2["angularVelocity"] = ang_vel[i];
            // j2["curvature"]

            // j2["pose"]["rotation"]["radians"] 
            j2["pose"]["translation"]["x"] = pos_x[i];
            j2["pose"]["translation"]["y"] = pos_y[i];

            j2["holonomicRotation"] = 0.0;
            j2["holonomicAngularVelocity"] = 0.0;
            j.push_back(j2);
        }

        return j;
    }
}
