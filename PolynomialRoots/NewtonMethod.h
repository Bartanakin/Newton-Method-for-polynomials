#pragma once
#include "Polynomial.h"
#include <cmath>
#include <random>

class NewtonMethod {
public:
    static constexpr double EPSILON = 1e-8;

    NewtonMethod(
        const Polynomial& p,
        unsigned int seed
    ) noexcept:
        p(p),
        dz({}),
        seed(seed) {
        this->init();
    }

    NewtonMethod(
        const Polynomial& p,
        unsigned int seed,
        ComplexNumber z0
    ) noexcept:
        p(p),
        dz({}),
        z(z0),
        seed(seed) {}

    void init() {
        this->max = 0.;
        for (int i = 0; i < this->p.deg(); i++) {
            this->max += std::abs(this->p[i]);
        }

        this->max /= std::abs(this->p[this->p.deg()]);

        // An estimation of maximal value of a root
        if (this->max < 200.) {
            this->max = 200.;
        }

        auto engine = std::default_random_engine(this->seed);
        auto dis = std::uniform_real_distribution<double>(-this->max, this->max);
        this->z = ComplexNumber(dis(engine), dis(engine));
        // this->z = ComplexNumber(dis(engine), 0.);

        this->seed--;
    }

    void iterate() {
        // calculate value at point z along with the first and second derivatives
        this->dz = std::move(this->p(this->z, this->p.deg()));
        this->det = this->dz[1].absSquare();
        if (this->det == 0.) {
            throw std::invalid_argument("Polynomial is not differentiable at " + std::string(this->z));
        }

        // Coefficients of D^(-1)(x_0)
        double a11, a12, a21, a22;
        a11 = this->dz[1].re;
        a12 = this->dz[1].im;
        a21 = -this->dz[1].im;
        a22 = this->dz[1].re;

        // x_(n+1) = x_n - D^(-1)(x_x) * f(x_n)
        this->z =
            this->z - ComplexNumber(a11 * this->dz[0].re + a12 * this->dz[0].im, a21 * this->dz[0].re + a22 * this->dz[0].im) * (1. / this->det);
    }

    std::vector<ComplexNumber> solve(
        bool multiple
    ) {
        double prev_L_local = std::numeric_limits<double>::max();
        double L_local = 0.;
        auto L = std::numeric_limits<double>::max();
        double prev_error = std::numeric_limits<double>::max();
        bool calPrev = false;
        ComplexNumber prev = this->z;
        auto R = std::numeric_limits<double>::max();
        do {
            // randomize one more time
            if ((this->z - prev).abs() > 10. * prev_error + NewtonMethod::EPSILON || this->z.abs() > 3. * this->max
                || L_local > 1. && std::abs(L_local - prev_L_local) <= NewtonMethod::EPSILON) {
                this->init();
            }

            if (calPrev) {
                prev_error = (this->z - prev).abs();
                prev = this->z;
                prev_L_local = L_local;
            }

            calPrev = true;
            this->iterate();

            // DN(x_0) - if this value is greater than 1. then definitely we are not close enough to the fixed point yet
            L_local = this->dz[0].abs() * this->dz[2].abs() / this->dz[1].absSquare();
            if (L_local < 1.) {
                // set L that we are looking for:  DN(x_0) <= L = (1 + DN(x_0)) / 2 <= 1
                R = 2. * (prev - this->z).abs() / (1. - L_local);

                if (R == 0.) {
                    break; // We are so close to the solution that there is no use of looking for a better approximation
                }

                auto F = 0.;
                auto F_prime = this->dz[1].abs();
                auto F_bis = 0.;
                auto Ri = 1.;
                for (int i = 0; i <= this->p.deg(); i++) {
                    F += this->dz[i].abs() * Ri * (1. / DiffAtPoint::factorials[i]);
                    if (i + 1 <= this->p.deg() && i > 0) {
                        F_prime -= this->dz[i + 1].abs() * Ri * (1. / DiffAtPoint::factorials[i]);
                    }

                    if (i + 2 <= this->p.deg()) {
                        F_bis += this->dz[i + 2].abs() * Ri * (1. / DiffAtPoint::factorials[i]);
                    }

                    Ri *= R;
                }

                if (F_prime < 0.) {
                    continue;
                }

                // upper approximation of L in the B(x_0, R)n set
                L = F * F_bis / (F_prime * F_prime);
            }

            // check if L <= (1 + DN(x_0)) / 2 and require low error
        } while (L > 0.5 * (1. + L_local) || (this->z - prev).abs() > 1e-8);

        // Check if the root is a real one. If not then return two roots where one is the conjugate of the second (*)
        if (multiple && std::abs(this->z.im) > R) {
            return {this->z, this->z.conjugate()};
        }

        return {this->z};
    }

    // This method divides a polynomial
    Polynomial getNextPolynomial(
        std::vector<ComplexNumber> zs
    ) const {
        auto n = this->p.deg();

        // There is one root only hence I am sure it is a real one
        if (zs.size() == 1) {
            std::vector<double> W(this->p.deg());
            W[n - 1] = this->p[n];
            for (size_t i = n - 1; i >= 1; i--) {
                W[i - 1] = this->p[i] + zs[0].re * W[i];
            }

            return Polynomial{std::move(W)};
        }

        if (n == 2) {
            return Polynomial{{this->p[2]}};
        }

        // There are two roots, so I know that one is the conjugate of the second (look at (*))
        double alpha = -2. * zs[0].re;
        double beta = std::pow(zs[0].re, 2) + std::pow(zs[0].im, 2);
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
    unsigned int seed;
};
