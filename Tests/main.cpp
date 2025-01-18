#include <iostream>
#include <cassert>
#include "../JacobiRotations/JacobiRotation.h" // Ensure this defines JacobiRotation properly
#include "../JacobiRotations/PivotSelectionStrategy.h" // Includes the strategies


// Function to print the selected pivot
int main() {
    // Create a test matrix
    Matrix testMatrix(4);
    testMatrix(0, 1) = 1.0;
    testMatrix(0, 2) = 2.0;
    testMatrix(0, 3) = -3.0;
    testMatrix(1, 2) = 4.0;
    testMatrix(1, 3) = 5.0;
    testMatrix(2, 3) = -6.0;

    std::cout << "Testing MaximalPivotSelectionStrategy" << std::endl;
    MaximalPivotSelectionStrategy maxStrategy;
    auto maxPivot = maxStrategy.selectNextPivot(testMatrix);
    maxPivot.print();
    assert(maxPivot.p == 2 && maxPivot.q == 3); // Verify largest element (6.0)

    std::cout << "\nTesting CyclicPivotSelectionStrategy" << std::endl;
    CyclicPivotSelectionStrategy cyclicStrategy;
    for (int i = 0; i < 6; ++i) { // Test cyclic iteration
        auto cyclicPivot = cyclicStrategy.selectNextPivot(testMatrix);
        cyclicPivot.print();
    }

    std::cout << "\nTesting RandomPivotSelectionStrategy" << std::endl;
    RandomPivotSelectionStrategy randomStrategy(42); // Use a fixed seed
    for (int i = 0; i < 3; ++i) { // Generate a few random pivots
        auto randomPivot = randomStrategy.selectNextPivot(testMatrix);
        randomPivot.print();
    }

    std::cout << "All tests completed successfully." << std::endl;
    return 0;
}
