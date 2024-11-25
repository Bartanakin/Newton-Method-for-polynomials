#pragma once
#include "Polynomial.h"
#include <cmath>
#include <random>

class NewtonMethod {
public:
    static constexpr double EPSILON = 1e-8;
    NewtonMethod(
        const Polynomial& p
    ) noexcept:
        p(p),
        dz({}) {
        this->init();
    }

    void init() {
        this->max = 0.;
        for (int i = 0; i < this->p.deg(); i++) {
            this->max += std::abs(this->p[i]);
        }

        this->max /= std::abs(this->p[this->p.deg()]);

        if (this->max < 1.) {
            this->max = 1.;
        }

        auto engine = std::default_random_engine(std::random_device{}());
        auto dis = std::uniform_real_distribution<double>(-this->max, this->max);
        this->z = ComplexNumber(dis(engine), dis(engine));
    }

    void iterate() {
        this->dz = this->p(this->z);
        this->det = this->dz.dx.re * this->dz.dy.im - this->dz.dy.re * this->dz.dx.im;
        if (this->det == 0.) {
            throw std::invalid_argument("Polynomial is not differentiable at " + std::string(this->z));
        }

        double a11, a12, a21, a22;
        a11 = this->dz.dy.im;
        a12 = -this->dz.dy.re;
        a21 = -this->dz.dx.im;
        a22 = this->dz.dx.re;

        this->z = this->z - ComplexNumber(a11 * this->dz.f.re + a12 * this->dz.f.im, a21 * this->dz.f.re + a22 * this->dz.f.im) * (1. / this->det);
    }

    void solve() {
        double a11, a12, a21, a22;
        double Mx, My;
        double M;
        double prev_norm = std::numeric_limits<double>::max();
        double norm = 0.;
        ComplexNumber prev;
        do {

            // randomize one more time
            if (this->z.abs() > 3. * this->max || norm > 1. && std::abs(norm - prev_norm) <= NewtonMethod::EPSILON) {
                this->init();
            }

            prev = this->z;
            this->iterate();

            M = this->det;
            Mx = this->dz.d2x2.re * this->dz.dy.im + this->dz.dx.re + this->dz.d2xy.im - this->dz.d2xy.re * this->dz.dx.im
                 - this->dz.dy.re * this->dz.d2x2.im;
            My = this->dz.d2xy.re * this->dz.dy.im + this->dz.dx.re + this->dz.d2y2.im - this->dz.d2y2.re * this->dz.dx.im
                 - this->dz.dy.re * this->dz.d2xy.im;

            a11 = (this->dz.d2xy.im * M - this->dz.dy.im * Mx) * this->dz.f.re + (this->dz.d2y2.im * M - this->dz.dy.im * My) * this->dz.f.im;
            // a12 = (-this->dz.d2xy.re * M + this->dz.dy.re * Mx) * this->dz.f.re - (this->dz.d2y2.re * M + this->dz.dy.re * My) * this->dz.f.im;
            // a21 = (-this->dz.d2x2.im * M + this->dz.dx.im * Mx) * this->dz.f.re - (this->dz.d2xy.im * M + this->dz.dx.im * My) * this->dz.f.im;
            a22 = (this->dz.d2x2.re * M - this->dz.dx.re * Mx) * this->dz.f.re + (this->dz.d2xy.re * M - this->dz.dx.re * My) * this->dz.f.im;

            prev_norm = norm;
            norm = (std::pow(a11, 2) + std::pow(a22, 2)) * std::pow(1. / M, 2);
        } while (std::abs(norm) >= 1. || (this->z - prev).abs() > 1e-8);
    }

    Polynomial getNextPolynomial() const {
        auto n = this->p.deg();
        if (std::abs(this->z.im) < 1e-6) {
            std::vector<double> W(this->p.deg());
            W[n - 1] = this->p[n];
            for (size_t i = n - 1; i >= 1; i--) {
                W[i - 1] = this->p[i] + this->z.re * W[i];
            }

            return Polynomial{std::move(W)};
        }

        if (n == 2) {
            return Polynomial{{this->p[2]}};
        }

        double alpha = -2. * this->z.re;
        double beta = std::pow(this->z.re, 2) + std::pow(this->z.im, 2);
        std::vector<double> W(this->p.deg() - 1);
        W[n - 2] = this->p[n];
        W[n - 3] = this->p[n - 1] - alpha * W[n - 2];
        for (size_t i = n - 2; i >= 2; i--) {
            W[i - 2] = this->p[i] - alpha * W[i - 1] - beta * W[i];
        }

        return Polynomial{std::move(W)};
    }

    const Polynomial& p;
    ComplexNumber z;
    DiffAtPoint dz;
    double det = 0.;
    double max;
};
