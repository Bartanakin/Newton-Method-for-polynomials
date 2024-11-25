//
// Created by bartanakin on 11/18/24.
//

#pragma once
#include "ComplexNumber.h"

struct DiffAtPoint {
    ComplexNumber f;
    ComplexNumber dx;
    ComplexNumber dy;
    ComplexNumber d2x2;
    ComplexNumber d2xy;
    ComplexNumber d2y2;

    DiffAtPoint(
        ComplexNumber z,
        ComplexNumber dx,
        ComplexNumber dy,
        ComplexNumber d2x2,
        ComplexNumber d2xy,
        ComplexNumber d2y2
    ):
        f(z),
        dx(dx),
        dy(dy),
        d2x2(d2x2),
        d2xy(d2xy),
        d2y2(d2y2) {}

    DiffAtPoint(
        ComplexNumber z
    ):
        f(z),
        dx({}),
        dy({}),
        d2x2({}),
        d2xy({}),
        d2y2({}) {}

    DiffAtPoint(
        double real
    ):
        DiffAtPoint(ComplexNumber(real)) {}

    DiffAtPoint operator+(
        const DiffAtPoint& second
    ) const {
        return {
            second.f + this->f,
            second.dx + this->dx,
            second.dy + this->dy,
            second.d2x2 + this->d2x2,
            second.d2xy + this->d2xy,
            second.d2y2 + this->d2y2
        };
    }

    DiffAtPoint& operator+=(
        const DiffAtPoint& other
    ) {
        *this = *this + other;

        return *this;
    }

    DiffAtPoint operator*(
        const DiffAtPoint& second
    ) const {
        return {
            second.f * this->f,
            this->f * second.dx + this->dx * second.f,
            this->f * second.dy + this->dy * second.f,
            this->d2x2 * second.f + this->dx * second.dx + this->dx * second.dx + this->f * second.d2x2,
            this->d2xy * second.f + this->dx * second.dy + this->dy * second.dx + this->f * second.d2xy,
            this->d2y2 * second.f + this->dy * second.dy + this->dy * second.dy + this->f * second.d2y2,
        };
    }

    DiffAtPoint& operator*=(
        const DiffAtPoint& other
    ) {
        *this = *this * other;

        return *this;
    }

    DiffAtPoint operator*(
        const ComplexNumber& z
    ) const {
        return {
            z * this->f,
            z * this->dx,
            z * this->dy,
            z * this->d2x2,
            z * this->d2xy,
            z * this->d2y2
        };
    }

    operator std::string() const {
        return "f = " + std::string(this->f) + "\nd = \n[" + std::string(this->dx) + "\n "
        + std::string(this->dy) + "]\nd2 = \n[" + std::string(this->d2x2) + " " + std::string(this->d2xy) + "\n "
        + std::string(this->d2xy) + " " + std::string(this->d2y2) + "]\n";
    }
};
