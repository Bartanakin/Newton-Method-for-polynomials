#pragma once
#include "JacobiRotation.h"
#include <random>

class PivotSelectionStrategyInterface {
public:
    virtual ~PivotSelectionStrategyInterface() = default;

    virtual JacobiRotation selectNextPivot(const Matrix& Mi) = 0;

    virtual std::string getStrategyName() const noexcept = 0;
};

class MaximalPivotSelectionStrategy: public PivotSelectionStrategyInterface {
public:
    JacobiRotation selectNextPivot(
        const Matrix& Mi
    ) override {
        auto max = std::numeric_limits<double>::min();
        JacobiRotation pivot(0, 1);
        for (int i = 0; i < Mi.getSize(); i++) {
            for (int j = i + 1; j < Mi.getSize(); j++) {
                if (max < std::abs(Mi(i, j))) {
                    pivot.p = i;
                    pivot.q = j;
                    max = std::abs(Mi(i, j));
                }
            }
        }

        return pivot;
    }

    std::string getStrategyName() const noexcept override {
        return "Jacobi method with maximal pivot selection";
    }
};

class CyclicPivotSelectionStrategy: public PivotSelectionStrategyInterface {
public:
    CyclicPivotSelectionStrategy() noexcept:
        curr_p(1),
        curr_q(0),
        M_ptr(nullptr) {}

    JacobiRotation selectNextPivot(
        const Matrix& Mi
    ) override {
        if (this->M_ptr == nullptr) {
            this->M_ptr = &Mi;
        } else if (this->M_ptr != &Mi) {
            throw std::invalid_argument("PivotSelectionStrategy: M_ptr is not the same as previously");
        }

        const int size = Mi.getSize();
        JacobiRotation pivot = {curr_p, curr_q};
        if (pivot.p > pivot.q) {
            pivot = {pivot.q, pivot.p};
        }

        ++curr_q;

        if (curr_q >= curr_p) {
            ++curr_p;
            curr_q = 0;
        }

        if (curr_p >= size) {
            curr_p = 1;
            curr_q = 0;
        }

        return pivot;
    }

    std::string getStrategyName() const noexcept override {
        return "Cyclic Jacobi method";
    }

private:
    int curr_p;
    int curr_q;
    const Matrix* M_ptr;
};

class RandomPivotSelectionStrategy: public PivotSelectionStrategyInterface {
public:
    RandomPivotSelectionStrategy(
        int seed
    ) noexcept:
        generator(seed) {}

    JacobiRotation selectNextPivot(
        const Matrix& Mi
    ) override {
        int p;
        int q;
        do {
            std::uniform_int_distribution<int> dist_p(0, Mi.getSize() - 1);
            p = dist_p(generator);

            std::uniform_int_distribution<int> dist_q(0, Mi.getSize() - 1);
            q = dist_q(generator);
        } while (p == q);

        if (p > q) {
            std::swap(p, q);
        }

        return {p, q};
    }

    std::string getStrategyName() const noexcept override {
        return "Jacobi method with random pivot selection";
    }

private:
    std::default_random_engine generator;
};
