#pragma once
#include "DiffAtPoint.h"
#include <vector>

class Polynomial {
public:
    Polynomial(
        std::vector<double> a
    ) noexcept:
        a(std::move(a)) {}

    Polynomial(Polynomial&&) noexcept = default;
    Polynomial(const Polynomial& second) noexcept = delete;

    Polynomial& operator=(const Polynomial& rhs) = delete;
    Polynomial& operator=(Polynomial&& rhs) = default;

    DiffAtPoint operator()(
        ComplexNumber z,
        unsigned int order = 0
    ) const {
        DiffAtPoint x{std::move(std::vector<ComplexNumber>(order + 1, 0.))};
        auto init = std::vector<ComplexNumber>(order + 1, 0.);
        init[0] = z;
        if (init.size() > 1) {
            init[1] = 1.;
        }

        DiffAtPoint dz{std::move(init)};
        for (int i = this->a.size() - 1; i > 0; i--) {
            x += this->a[i];
            x *= dz;
        }

        x += this->a[0];

        return x;
    }

    operator std::string() const {
        std::string s = "";
        for (int i = this->a.size() - 1; i > 0; i--) {
            s += std::to_string(this->a[i]) + "z^" + std::to_string(i) + " + ";
        }

        s += std::to_string(this->a[0]);

        return s;
    }

    int deg() const { return this->a.size() - 1; }

    const double operator[](
        size_t i
    ) const {
        return this->a[i];
    }

private:
    std::vector<double> a; // coefficients
};
