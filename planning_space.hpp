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

#include "utilities.hpp"
#include "planner.hpp"
#include "bezier.hpp"
#include "LinearJointVelocityVarying.hpp"

using point_set = std::unordered_set<Vector2f, hash_vector2f>;


/*class planning_space {
    public:
        point_set sample_free(const int n);
        float cost(const Vector2f a, const Vector2f b) const;
        point_set near(const Vector2f b, const point_set& nodes, const float dist) const;
        void register_free_space_constraint(std::function<bool(Vector2f)>);

        std::optional<std::vector<Vector2f>> fast_marching_trees(
            const Vector2f& x_init, 
            const Vector2f& x_goal, 
            const int n, 
            const float rn
        );

        planning_space(const bounding_rect& br);

        std::vector<obstacle> obstacles;
        bounding_rect bound_rect;
        halton_state x_state;
        halton_state y_state;
};*/

class planning_space {
    public:
        point_set sample_free(const int n);

        bool is_free(const Vector2f& p);
        bool is_free_space_allocated(const Vector2f& p);
        std::tuple<bool, obstacle> is_obstacle(const Vector2f& p);

        bool is_same_obstacle_fuzzy(const Vector2f& a, const Vector2f& b, const float epsilon);

        float cost(const Vector2f a, const Vector2f b) const;
        point_set near(const Vector2f b, const point_set& nodes, const float dist) const;
        std::variant<std::vector<Vector2f>> fast_marching_trees(const Vector2f& x_init, const Vector2f& x_goal, const int n, const float rn);

        planning_space(const bounding_rect& br);

        std::vector<obstacle> obstacles;
        std::vector<std::function<bool(Vector2f)>> free_space_allocations;
        bounding_rect bound_rect;
        halton_state x_state;
        halton_state y_state;
};


using point_set = std::unordered_set<Vector2f, hash_vector2f>;

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

std::variant<std::vector<Vector2f>> planning_space::fast_marching_trees(const Vector2f& x_init, const Vector2f& x_goal, const int n, const float rn) {
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

bool planning_space::is_free_space_allocated(const Vector2f& p) {
    for (auto& free_space_allocation : free_space_allocations) {
        if (free_space_allocation(p)) return true;
    }
    return false;
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

point_set planning_space::sample_free(const int n) {
    // technically this should be a set, but the halton sequence is guaranteed to not repeat
    // so we can avoid element checks for a set
    point_set pts = {Vector2f(0,0)};
    pts.reserve(n);

    while (n > pts.size()) {
        // TODO: consider using different bases, or expose the bases to the user
        //done
        std::vector<float> x_test = halton(2, n - pts.size(), x_state);
        std::vector<float> y_test = halton(3, n - pts.size(), y_state);
        std::cout<<"X_test base is 2 and y_test base is 3"<<std::endl;
        for (int i = 0; i < x_test.size(); ++i) {
            const Vector2f test((bound_rect.x_max-bound_rect.x_min)*(x_test[i])+bound_rect.x_min, (bound_rect.y_max-bound_rect.y_min)*(y_test[i])+bound_rect.y_min);
            std::cout << "testing point " << test.x() << " " << test.y() << std::endl;
            if (is_free(test)) pts.insert(test);
        }
    }

    return pts;
}