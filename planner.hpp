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
#include "planning_space.hpp"
#include "bezier.hpp"
#include "LinearJointVelocityVarying.hpp"

class planner {
    public:
        planning_space ps;
        std::stack<Vector2f> goal_point_cache;
        int64_t past_id;

        Vector2f pick_next_goal_point(const Vector2f& start, const int n, const float search_radius);
        // void update_planning_space();
};


using point_set = std::unordered_set<Vector2f, hash_vector2f>;

Vector2f planner::pick_next_goal_point(const Vector2f& start, const int n, const float search_radius) {
    // const point_set free_pts = sample_free(n);

    // trace circle around start
    bool is_free_past = false;
    const float epsilon = 0.001;
    std::vector<Vector2f> test_pts;
    for (std::size_t i = 0; i < n; i++) {
        const float angle = 2 * std::numbers::pi * (i * 1.0f/n);
        Vector2f test(search_radius * std::cos(angle), search_radius * std::sin(angle));
        bool isf = ps.is_free(test);
        // Since our local area of known space will be circular, any change from free to non-free along the circle will have to be due to an obstacle
        // TODO: fix condition in the loop below
        for (std::size_t j = 0; j < n / 10; j++) {
            Vector2f test_uk((search_radius + epsilon) * std::cos(angle), (search_radius + epsilon) * std::sin(angle));
            bool is_uk = !ps.is_free_space_allocated(test_uk);
            if (isf != is_free_past && i != 0 && is_uk) {
                test_pts.push_back(test);
            }
        }
        // if (isf != is_free_past && i != 0) {
        //     test_pts.push_back(test);
        // }
        is_free_past = isf;
    }

    if (test_pts.empty()) {
        if (goal_point_cache.empty()) {
            // failed
        }
        // pull point from cache
        test_pts.push_back(goal_point_cache.top());
        goal_point_cache.pop();
    }

    // grab a point and send rest to cache
    const Vector2f next_point = test_pts.back();
    // TODO: grab the point on the same obstacle
    test_pts.pop_back();
    return next_point;

    // pick possible goal points and the use heuristic to select one
    // for (auto& obstacle : ps.obstacles) {
    //     // trace each obstacle to find points that border free space
    //     for (auto& line : obstacle) {
    //         const Vector2f a = std::get<0>(line);
    //         const Vector2f b = std::get<1>(line);
    //         const float dist = (b - a).norm();
    //         const float n = 100; // TODO: think abt this

    //     }
    // }
}