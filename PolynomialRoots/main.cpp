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

    // auto seed = std::random_device{}();
    unsigned int seed = 3065439501;
    std::cout << "Seed: " << seed << std::endl;
    std::vector<ComplexNumber> roots(d);
    int rootsCount = 0;
    while (p.deg() > 1) {
        // Find an estimation of the root on a reduced polynomial
        NewtonMethod N(p, seed);
        auto z = N.solve(false);
        // Use it as a starting point for the original polynomial
        NewtonMethod N_original{p_original, seed, z[0]};
        auto zs = N_original.solve(true);
        for (int i = rootsCount; i < rootsCount + zs.size(); i++) {
            roots[i] = zs[i - rootsCount];
        }

        rootsCount += zs.size();
        p = N.getNextPolynomial(zs);
    }

    if (p.deg() == 1) {
        NewtonMethod N_original{p_original, seed,{- p[0] / p[1]}};
        roots[roots.size() - 1] = N_original.solve(false)[0];
    }

    std::cout << "znaleziono " << roots.size() << " pierwiastkow." << std::endl;
    for (int i = 0; i < d; i++) {
        std::cout << std::string(roots[i]) << ", ";
    }

    std::cout << std::endl;
}
