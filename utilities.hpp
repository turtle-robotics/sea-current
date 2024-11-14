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

#include "planning_space.hpp"
#include "planner.hpp"
#include "bezier.hpp"
#include "LinearJointVelocityVarying.hpp"

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

struct bounding_rect {
        float x_max;
        float x_min;
        float y_max;
        float y_min;

        inline bool contains(const Vector2f& point);
        inline void enclose_point(const Vector2f& pt);

        bounding_rect(
            const float x_max, 
            const float x_min, 
            const float y_max, 
            const float y_min
        ) : 
            x_max(x_max), 
            x_min(x_min), 
            y_max(y_max), 
            y_min(y_min) 
        {}
};

struct halton_state {
    int f = 0;
    int i = 0;
    halton_state(int f, int i) : f(f), i(i) {}
    halton_state() : f(0), i(0) {}
};

struct obstacle {
    std::vector<std::tuple<Vector2f, Vector2f>> lines;
    std::vector<Vector2f> vertices;
    bounding_rect bound_rect = {0, 0, 0, 0};

    bool contains(const Vector2f& point);

    obstacle(std::vector<Vector2f> vertices) : vertices(vertices);
    obstacle(std::vector<Vector2f> vertices, std::vector<std::tuple<int, int>> edges) : vertices(vertices); //check
};

struct hash_vector2f {
    size_t operator()(const Vector2f v) const;
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

    chebpoly(
        VectorXf coeffs, 
        const float xmin, 
        const float xmax
    ) : 
        coeffs(coeffs), 
        xmin(xmin), 
        xmax(xmax) 
    {}
};

struct velocity_profile {
    std::vector<VectorXf> pos;
    std::vector<VectorXf> vel;
    std::vector<VectorXf> acc;
    toppra::Vector time;

    velocity_profile(
        std::vector<VectorXf> pos, 
        std::vector<VectorXf> vel, 
        std::vector<VectorXf> acc, 
        toppra::Vector time
    ) : 
        pos(pos), 
        vel(vel), 
        acc(acc), 
        time(time) 
    {}
};

inline Vector2f calc_start_tangent(const Vector2f& W_0, const Vector2f& W_1, const float theta);

chebpoly chebfit(const VectorXf& x, const VectorXf& y, const int degree);
VectorXf chebeval(const VectorXf& x, const chebpoly& b, const int degree);


inline bool bounding_rect::contains(const Vector2f& point) {
    return point.x() <= x_max &&
           point.x() >= x_min &&
           point.y() <= y_max &&
           point.y() >= y_min;
}

// grows the bounding rect to enclose a point
inline void bounding_rect::enclose_point(const Vector2f& pt) {
    if (pt.x() > x_max) x_max = pt.x();
    if (pt.x() < x_min) x_min = pt.x();
    if (pt.y() > y_max) y_max = pt.y();
    if (pt.y() < y_min) y_min = pt.y();
}

bool obstacle::contains(const Vector2f& point) {
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

obstacle::obstacle(std::vector<Vector2f> vertices) : vertices(vertices) {
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

obstacle::obstacle(std::vector<Vector2f> vertices, std::vector<std::tuple<int, int>> edges) : vertices(vertices) {
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

size_t hash_vector2f::operator()(const Vector2f v) const {
    float fa = v.x();
    float fb = v.y();
    const int32_t a = reinterpret_cast<int32_t&>(fa);
    const int32_t b = reinterpret_cast<int32_t&>(fb);
    return std::hash<int32_t>()(a) ^ std::hash<int32_t>()(b);
}

std::vector<float> halton(const int b, const int n, halton_state& state) {
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

inline float pt_dist(const Vector2f& u, const Vector2f& v={0,0}) {
    return std::sqrt(std::pow(v.x() - u.x(), 2) + std::pow(v.y() - u.y(), 2));
}

inline float cross2d(Vector2f u, Vector2f v) {
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

// distance from a point to a line segment
float dist_pt_line(const Vector2f& p1, const Vector2f& p2, const Vector2f a) {
    const float x0 = a.x();
    const float y0 = a.y();
    const float x1 = p1.x();
    const float y1 = p1.y();
    const float x2 = p2.x();
    const float y2 = p2.y();
    const float dist = std::abs((x2 - x1) * (y1 - y0) - (x1 - x0) * (y2 - y1)) / (p2 - p1).norm();
    return dist;
}

inline float tangent_magnitude(const Vector2f& W_0, const Vector2f& W_1, const Vector2f& W_2) {
    return 0.5 * std::min(pt_dist(W_0, W_1), pt_dist(W_1, W_2));
}

// inline Vector2f calc_start_tangent(const Vector2f& W_0, const Vector2f& W_1, const float theta);

inline Vector2f calc_start_tangent(const Vector2f& W_0, const Vector2f& W_1, const float theta) {
    // Vector2f tmp(std::cos(theta), std::sin(theta));
    // std::cout << "wtfff " << tmp.x() << " " << tmp.y() << std::endl;
    // std::cout << tangent_magnitude(W_0, W_1, W_0) << std::endl;
    // Vector2f tmp2 = tangent_magnitude(W_0, W_1, W_0) * Vector2f(std::cos(theta), std::sin(theta));
    // std::cout << "wtfff2 " << tmp2.x() << " " << tmp2.y() << std::endl;
    return tangent_magnitude(W_0, W_1, W_0) * Vector2f(std::cos(theta), std::sin(theta));
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

// chebpoly chebfit(const VectorXf& x, const VectorXf& y, const int degree);

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

// VectorXf chebeval(const VectorXf& x, const chebpoly& b, const int degree);

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

inline Vector2f calc_end_tangent(const Vector2f& W_1, const Vector2f W_2) {
    return tangent_magnitude(W_1, W_2, W_1) * (W_2 - W_1).normalized();
}