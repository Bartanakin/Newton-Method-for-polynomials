#pragma once
#include <cmath>
#include "Matrix.h"

class JacobiRotation {
public:
    int p;
    int q;

    JacobiRotation() noexcept = default;

    JacobiRotation(
        int p,
        int q
    ):
        p(p),
        q(q) {
        if (p == q) {
            throw std::invalid_argument("JacobiRotation: p == q");
        }
    }

    // pair of sin(theta), cos(theta)
    std::pair<double, double> getAngle(const Matrix& M) const noexcept {
        auto p = this->p;
        auto q = this->q;
        auto A = M(p, p);
        auto D = M(q, q);
        if (A == D) {
            return {std::sin(M_PI / 4.), std::cos(M_PI / 4.)};
        }

        auto B = M(p, q);
        auto tan2theta = 2. * B / (A - D);
        auto tanTheta = (-1. + std::sqrt(1. + tan2theta * tan2theta)) / tan2theta;
        auto tanTheta2 = (-1. - std::sqrt(1. + tan2theta * tan2theta)) / tan2theta;
        auto sinTheta = tanTheta / std::sqrt(1 + tanTheta * tanTheta);
        auto cosTheta = std::sqrt(1. - sinTheta * sinTheta);
        auto sinTheta2 = tanTheta2 / std::sqrt(1 + tanTheta2 * tanTheta2);
        auto cosTheta2 = std::sqrt(1. - sinTheta2 * sinTheta2);

        return {sinTheta, cosTheta};
    }

    void print() const noexcept {
        std::cout << "Pivot: (" << this->p << ", " << this->q << ")" << std::endl;
    }
};