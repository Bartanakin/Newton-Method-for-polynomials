#include <array>
#include <chrono>

#include "JacobiDiagonalizationMethod.h"

int main() {
    auto matrix = Matrix::createFromStdcin();

    std::array<std::unique_ptr<PivotSelectionStrategyInterface>, 3> pivotSelectionStrategies;
    pivotSelectionStrategies[0] = std::make_unique<CyclicPivotSelectionStrategy>();
    pivotSelectionStrategies[1] = std::make_unique<MaximalPivotSelectionStrategy>();
    // pivotSelectionStrategies[2] = std::make_unique<RandomPivotSelectionStrategy>(std::random_device()());
    pivotSelectionStrategies[2] = std::make_unique<RandomPivotSelectionStrategy>(2137);

    for (int i = 0; i < pivotSelectionStrategies.size(); i++) {
        #define ITERATIONS 1000
        auto start = std::chrono::steady_clock::now();
        std::cout << pivotSelectionStrategies[i]->getStrategyName() << " ---------------------------->" << std::endl;
        auto jdm = JacobiDiagonalizationMethod(1e-6, ITERATIONS, std::move(pivotSelectionStrategies[i]));
        auto [Q_t, D, Q, last_iteration] = jdm.run(matrix);
        auto end = std::chrono::steady_clock::now();

        std::cout << "Time passed: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << " microseconds." << std::endl;
        if (last_iteration > ITERATIONS) {
            std::cout << "Diagonalization unsuccessful after " << last_iteration << " iterations." << std::endl;
            std::cout << std::endl;

            continue;
        }

        std::cout << "Diagonalization successful after " << last_iteration << " iterations." << std::endl;
        std::cout << "M_0: " << std::endl;
        matrix.display();

        std::cout << "M_0 = Q^{T}DQ, where: " << std::endl;
        std::cout << "Q^T " << std::endl;
        Q_t.display();

        std::cout << "D = " << std::endl;
        D.display();

        std::cout << "Q = " << std::endl;
        Q.display();

        std::cout << std::endl;
    }

    return 0;
}
