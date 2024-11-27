#include <iostream>
#include <vector>
#include "NewtonMethod.h"

int main() {
    int d;
    std::cin >> d;
    std::string comment;
    std::getline(std::cin, comment);
    std::vector<double> coefficients(d + 1);
    for (int i = 0; i < coefficients.size(); i++) {
        std::cin >> coefficients[i];
    }

    Polynomial p_original{coefficients};
    Polynomial p{coefficients};

    std::cout << std::string(p) << std::endl;

    std::vector<ComplexNumber> roots(d);
    while (p.deg() > 1) {
        // Find an estimation of the root on a reduced polynomial
        NewtonMethod N(p);
        auto z = N.solve();
        // Use it as a starting point for the original polynomial
        NewtonMethod N_original{p_original, z};
        z = N_original.solve();
        if (std::abs(z.im) > 1e-6) {
            roots[p.deg() - 1] = z;
            roots[p.deg() - 2] = z.conjugate();
        } else {
            roots[p.deg() - 1] = z;
        }

        p = N.getNextPolynomial();
    }

    if (p.deg() == 1) {
        NewtonMethod N_original{p_original, {- p[0] / p[1]}};
        roots[p.deg() - 1] = N_original.solve();
    }

    for (int i = 0; i < d; i++) {
        std::cout << std::string(roots[i]) << ", ";
    }

    std::cout << std::endl;
}
