//
// Created by bartanakin on 11/18/24.
//

#pragma once
#include "ComplexNumber.h"

struct DiffAtPoint {
    static std::vector<unsigned int> factorials;
    std::vector<ComplexNumber> z;

    DiffAtPoint(
        std::vector<ComplexNumber> z
    ):
        z(std::move(z)) {

        auto oldSize = factorials.size();
        if (factorials.size() < this->z.size()) {
            factorials.resize(this->z.size());
        }

        for (unsigned int i = oldSize; i < factorials.size(); i++) {
            factorials[i] = factorials[i - 1] * i;
        }
    }

    DiffAtPoint operator+(
        const DiffAtPoint& second
    ) const {
        if (second.z.size() != z.size()) {
            throw std::invalid_argument("size mismatch");
        }

        std::vector<ComplexNumber> temp(z.size());
        for (int i = 0; i < z.size(); i++) {
            temp[i] = z[i] + second.z[i];
        }

        return {std::move(temp)};
    }

    DiffAtPoint& operator+=(
        const DiffAtPoint& other
    ) {
        *this = *this + other;

        return *this;
    }

    DiffAtPoint& operator+=(
        const ComplexNumber& z
    ) {
        this->z[0] += z;

        return *this;
    }

    DiffAtPoint operator*(
        const DiffAtPoint& second
    ) const {
        if (second.z.size() != z.size()) {
            throw std::invalid_argument("size mismatch");
        }

        std::vector<ComplexNumber> temp(z.size());
        for (int k = 0; k < z.size(); k++) {
            ComplexNumber sum{};
            for (int i = 0; i <= k; ++i) {
                sum += z[i] * second.z[k - i];
            }

            temp[k] = sum;
        }

        return {std::move(temp)};
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

        std::vector<ComplexNumber> temp(this->z.size());
        for (int i = 0; i < this->z.size(); i++) {
            temp[i] = z * this->z[i];
        }

        return {std::move(temp)};
    }

    ComplexNumber operator[](size_t i) const {
        if (i >= z.size()) {
            throw std::invalid_argument("index out of range");
        }

        return this->z[i] * factorials[i];
    }

    operator std::string() const {
        std::string s = "[";
        for (int i = 0; i < this->z.size(); i++) {
            s += this->operator[](i);
            s += " ";
        }

        s += "]";

        return s;
    }
};

std::vector<unsigned int> DiffAtPoint::factorials = {1, 1};
