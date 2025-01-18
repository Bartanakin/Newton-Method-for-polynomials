#pragma once

#include "JacobiRotation.h"
#include "PivotSelectionStrategy.h"
#include <memory>

class JacobiDiagonalizationMethod {
public:
    JacobiDiagonalizationMethod(
        const double EPS,
        const unsigned int maxIterations,
        std::unique_ptr<PivotSelectionStrategyInterface> pivotSelectionStrategy_ptr
    ) noexcept:
        EPS(EPS),
        MAX_ITERATIONS(maxIterations),
        pivotSelectionStrategy_ptr(std::move(pivotSelectionStrategy_ptr)) {}

    std::tuple<Matrix, Matrix, Matrix, unsigned int> run(
        const Matrix& M0
    ) const {
        const auto n = M0.getSize();
        Matrix R(n);

        // sum of all outdiagonal element squares
        double squareSum = 0.;
        for (int i = 0; i < n; i++) {
            R(i, i) = 1.;
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    continue;
                }

                squareSum += std::pow(M0(i, j), 2);
            }
        }

        auto Mi = M0;
        std::vector<std::tuple<JacobiRotation, std::tuple<double, double>>> rotations;
        std::vector<double> pElements(n);
        std::vector<double> qElements(n);

        unsigned int iter = 0;
        while (squareSum >= this->EPS && iter <= this->MAX_ITERATIONS) {
            auto pivot = this->pivotSelectionStrategy_ptr->selectNextPivot(Mi);
            const auto& p = pivot.p;
            const auto& q = pivot.q;
            if (std::abs(Mi(p,q)) < this->EPS / (n * n)) {
                continue;
            }

            auto [s, c] = pivot.getAngle(Mi);
            rotations.push_back({pivot, {s,c }});

            for (int i = 0; i < n; i++) {
                // M update calculation and buffering
                double temp_p;
                double temp_q;
                if (i == q) {
                    temp_p = (c*c - s*s) * Mi(p, q) + s * c * (Mi(q, q) - Mi(p, p));
                    temp_q = c * c * Mi(q, q) - 2. * s * c * Mi(p, q) + s * s * Mi(p, p);
                    squareSum += std::pow(temp_p, 2) - std::pow(Mi(p, i), 2);
                } else if (i == p) {
                    temp_p = c * c * Mi(p, p) + 2. * s * c * Mi(p, q) + s * s * Mi(q, q);
                    temp_q = (c*c - s*s) * Mi(p, q) + s * c * (Mi(q, q) - Mi(p, p));
                    squareSum += std::pow(temp_q, 2) - std::pow(Mi(i, q), 2);
                } else {
                    temp_p = c * Mi(i, p) + s * Mi(i, q);
                    temp_q = c * Mi(i, q) - s * Mi(i, p);
                }

                pElements[i] = temp_p;
                qElements[i] = temp_q;

                // R update
                double prev_rp = R(i, p);
                double prev_rq = R(i, q);
                R(i, p) = c * prev_rp + s * prev_rq;
                R(i, q) = c * prev_rq - s * prev_rp;
            }

            // updating entries
            for (int i = 0; i < n; i++) {
                Mi(i, p) = pElements[i];
                Mi(p, i) = pElements[i];

                Mi(i, q) = qElements[i];
                Mi(q, i) = qElements[i];
            }

            ++iter;
        }

        return {R.transpose(), Mi, R, iter};
    }

private:
    const double EPS;
    const unsigned int MAX_ITERATIONS;
    std::unique_ptr<PivotSelectionStrategyInterface> pivotSelectionStrategy_ptr;
};
