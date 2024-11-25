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

    Polynomial p{coefficients};

    std::vector<ComplexNumber> roots(d);
    while (p.deg() > 1) {
        NewtonMethod N(p);
        N.solve();
        if (std::abs(N.z.im) > 1e-6) {
            roots[p.deg() - 1] = N.z;
            roots[p.deg() - 2] = N.z.conjugate();
        } else {
            roots[p.deg() - 1] = N.z;
        }

        p = N.getNextPolynomial();
    }

    if (p.deg() == 1) {
        roots[p.deg() - 1] = - p[0] / p[1];
    }

    for (int i = 0; i < d; i++) {
        std::cout << std::string(roots[i]) << ", ";
    }

    std::cout << std::endl;
}
