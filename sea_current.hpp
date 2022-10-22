#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <numbers>
#include <algorithm>

#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>


namespace turtle::sc {

    using namespace Eigen;

    using namespace std::complex_literals;

    inline void dbg_assert(bool cnd, const char* msg) {
        #ifdef DEBUG
        if (!cnd) {
            std::cerr << "assertion failed: " << msg << std::endl;
            std::exit(1);
        }
        #endif
    }

    class bezier_spline {
        public:
            std::vector<std::vector<Vector2f>> ctrl_pts;
            Matrix<float, Dynamic, 2> pts;

            bezier_spline(std::vector<Vector2f>& ctrl_pts, Matrix<float, Dynamic, 2>& pts);
            bezier_spline(std::vector<std::vector<Vector2f>>& ctrl_pts, Matrix<float, Dynamic, 2>& pts);

            static bezier_spline bezier_curve(std::vector<Vector2f>& ctrl_pts, std::vector<float>& positions);
            static bezier_spline bezier_curve(std::vector<Vector2f>& ctrl_pts, float precision);

            static std::vector<Vector2f> transform_ctrl_pts(std::vector<Vector2f>& ctrl_pts);
            static std::vector<std::complex<float>> omega_table(int degree);
    };

    bezier_spline::bezier_spline(std::vector<std::vector<Vector2f>>& ctrl_pts, Matrix<float, Dynamic, 2>& pts) : ctrl_pts(ctrl_pts), pts(pts) {}

    bezier_spline bezier_spline::bezier_curve(std::vector<Vector2f>& ctrl_pts, std::vector<float>& positions) {
        const int degree = ctrl_pts.size() - 1;

        ctrl_pts = transform_ctrl_pts(ctrl_pts);

        std::vector<std::complex<float>> omegas = omega_table(degree);

        FFT<float> fft;

        Matrix<std::complex<float>, Dynamic, 2> U;
        U.setZero(ctrl_pts.size(), 2);

        for (std::size_t i = 0; i <= degree; ++i) {
            U(i, 0) = ctrl_pts[i].x();
            U(i, 1) = ctrl_pts[i].y();
        }
        Matrix<std::complex<float>, Dynamic, 2> Q;
        Q.setZero(ctrl_pts.size(), 2);

        Eigen::VectorXcf tmp_fft(ctrl_pts.size());
        fft.inv(tmp_fft, U.col(0));
        Q.col(0) = tmp_fft;

        fft.inv(tmp_fft, U.col(1));
        Q.col(1) = tmp_fft;

        Matrix<float, Dynamic, 2> B;
        B.setZero(positions.size(), 2);

        for (int i = 0; i < positions.size(); ++i) {
            float s = positions[i];
            Vector2cf sum;
            for (int k = 0; k <= degree; ++k) {
                std::complex<float> tmp = std::pow((1.0f+0if) + s*(omegas[k] - (1.0f+0if)), degree);
                sum.x() += Q(k, 0) * tmp;
                sum.y() += Q(k, 1) * tmp;
            }
            B(i, 0) = sum.x().real();
            B(i, 1) = sum.y().real();
        }

        std::vector<std::vector<Vector2f>> tmp(1);
        tmp[0] = ctrl_pts;
        return bezier_spline(tmp, B);
    }

    bezier_spline bezier_spline::bezier_curve(std::vector<Vector2f>& ctrl_pts, const float precision) {
        dbg_assert(precision < 1 && precision > 0, "spline percision must be in (0, 1)");

        const int n = static_cast<int>(std::round(1.0f / precision));

        dbg_assert(std::abs(n - (1.0f / precision)) < 0.0001,
                    "1/percision must be (close) to an integer, for arbitrary position values use the other bezier_curve method");

        std::vector<float> positions(n+1);
        for (std::size_t i = 0; i <= n; ++i) {
            positions[i] = std::min(i * precision, 1.0f);
        }

        return bezier_curve(ctrl_pts, positions);
    }

    std::vector<Vector2f> bezier_spline::transform_ctrl_pts(std::vector<Vector2f>& ctrl_pts) {
        std::vector<Vector2f> pts_out(ctrl_pts.size());
        for (std::size_t i = 0; i < ctrl_pts.size(); ++i) {
            pts_out[i] = Vector2f(ctrl_pts[i].x() + ctrl_pts[i].y(), ctrl_pts[i].x() - ctrl_pts[i].y());
        }
        return pts_out;
    }

    std::vector<std::complex<float>> bezier_spline::omega_table(const int degree) {
        using std::complex;

        std::vector<std::complex<float>> omegas(degree+1);
        omegas[0] = 1.0f+0if;
        for (std::size_t i = 1; i <= degree; ++i) {
            omegas[i] = pow(exp((std::complex<float>(std::numbers::pi) * -2if) / std::complex<float>(degree+1)), i);
        }

        return omegas;
    }
}
