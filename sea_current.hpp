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

            inline int num_pts() const;

            float arclength() const;
            std::vector<float> curvature() const;
            bezier_spline resample(std::vector<float>& arclens) const;


            static bezier_spline join_splines(std::vector<bezier_spline>& splines);

            static bezier_spline bezier_curve(const std::vector<Vector2f>& og_ctrl_pts, std::vector<float>& positions);
            static bezier_spline bezier_curve(std::vector<Vector2f>& ctrl_pts, float precision);

            static std::vector<std::complex<float>> omega_table(int degree);
    };

    bezier_spline::bezier_spline(std::vector<std::vector<Vector2f>>& ctrl_pts, Matrix<float, Dynamic, 2>& pts) : ctrl_pts(ctrl_pts), pts(pts) {}

    inline int bezier_spline::num_pts() const {
        return pts.rows();
    }

    // it's ok that this is static as a new matrix would have to be allocated either way
    bezier_spline bezier_spline::join_splines(std::vector<bezier_spline>& splines) {
        std::vector<std::vector<Vector2f>> ctrl_pts;

        int num_pts = 0;
        for (bezier_spline spline : splines) {
            num_pts += spline.num_pts();
            for (std::vector<Vector2f> cpv : spline.ctrl_pts) {
                ctrl_pts.push_back(cpv);
            }
        }

        Matrix<float, Dynamic, 2> joined_pts;
        joined_pts.setZero(num_pts, 2);

        int n = 0;
        for (bezier_spline spline : splines) {
            joined_pts.block(n, 0, spline.num_pts(), 2);
            n += spline.num_pts();
        }

        return bezier_spline(ctrl_pts, joined_pts);
    }

    bezier_spline bezier_spline::bezier_curve(const std::vector<Vector2f>& og_ctrl_pts, std::vector<float>& positions) {
        const int degree = og_ctrl_pts.size() - 1;

        //const std::vector<Vector2f> ctrl_pts = bezier_spline::transform_ctrl_pts(og_ctrl_pts);
        const std::vector<Vector2f> ctrl_pts = og_ctrl_pts;

        const std::vector<std::complex<float>> omegas = bezier_spline::omega_table(degree);

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
            for (int k = 0; k <= degree; ++k) {
                std::complex<float> tmp = std::pow((1.0f+0if) + s*(omegas[k] - (1.0f+0if)), degree);
                B(i, 0) += (Q(k, 0) * tmp).real();
                B(i, 1) += (Q(k, 1) * tmp).real();
            }
        }

        std::vector<std::vector<Vector2f>> tmp(1);
        tmp[0] = og_ctrl_pts;
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

    std::vector<std::complex<float>> bezier_spline::omega_table(const int degree) {
        using std::complex;

        std::vector<std::complex<float>> omegas(degree+1);
        omegas[0] = 1.0f+0if;
        for (std::size_t i = 1; i <= degree; ++i) {
            omegas[i] = pow(exp((complex<float>(std::numbers::pi) * -2if) / complex<float>(degree+1)), i);
        }

        return omegas;
    }
}
