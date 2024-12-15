#pragma once

#include <iostream>
#include <math.h>

struct ComplexNumber {
    double re;
    double im;

    ComplexNumber(
        double real,
        double imaginary
    ):
        re(real),
        im(imaginary) {}

    ComplexNumber(
        double real
    ):
        ComplexNumber(real, 0.) {}

    ComplexNumber():
        ComplexNumber(0.f) {}

    ComplexNumber operator+(
        const ComplexNumber& other
    ) const {
        return {this->re + other.re, this->im + other.im};
    }

    ComplexNumber operator+(
        double real
    ) const {
        return {this->re + real, this->im};
    }

    ComplexNumber& operator+=(
        const ComplexNumber& other
    ) {
        *this = *this + other;

        return *this;
    }

    ComplexNumber operator*(
        const ComplexNumber& other
    ) const {
        return {this->re * other.re - this->im * other.im, this->im * other.re + this->re * other.im};
    }

    ComplexNumber& operator*=(
        const ComplexNumber& other
    ) {
        *this = *this * other;

        return *this;
    }

    operator std::string() const { return std::to_string(this->re) + " + " + std::to_string(this->im) + "i"; }

    ComplexNumber operator-() const { return {-this->re, -this->im}; }

    ComplexNumber operator-(
        const ComplexNumber& second
    ) const {
        return *this + -second;
    }

    double abs() const { return std::sqrt(this->absSquare()); }

    double absSquare() const { return std::pow(this->re, 2) + std::pow(this->im, 2); }

    ComplexNumber conjugate() const { return {this->re, -this->im}; }
};

ComplexNumber operator""_i(
    long double imaginary
) {
    return {0., static_cast<double>(imaginary)};
}

ComplexNumber operator+(
    double real,
    ComplexNumber z
) {
    return z + real;
}
