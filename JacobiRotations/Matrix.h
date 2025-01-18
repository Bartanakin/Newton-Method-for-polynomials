#pragma once

#include <iomanip>
#include <iostream>
#include <vector>
#include <stdexcept>

class Matrix {
public:
    Matrix(
        int n
    ):
        size(n) {
        matrix.resize(size, std::vector<double>(size, 0));
    }

    void display() const {
        int maxWidth = 0;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                int width = std::to_string(matrix[i][j]).length();
                if (width > maxWidth) {
                    maxWidth = width;
                }
            }
        }

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                std::cout << std::setw(maxWidth + 1) << matrix[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }

    Matrix transpose() const noexcept {
        Matrix m(this->getSize());
        for (int i = 0; i < m.getSize(); i++) {
            for (int j = 0; j < m.getSize(); j++) {
                m.matrix[i][j] = this->matrix[j][i];
            }
        }

        return m;
    }

    static Matrix createFromStdcin() {
        int n;
        std::cin >> n;
        std::string comment;
        std::getline(std::cin, comment);
        Matrix m(n);
        for (int i = 0; i < m.size; i++) {
            for (int j = 0; j < m.size; j++) {
                std::cin >> m.matrix[i][j];
            }
        }

        return m;
    }

    double& operator()(int row, int col) {
        validateIndices(row, col);

        return matrix[col][row];
    }

    double operator()(int row, int col) const {
        validateIndices(row, col);

        return matrix[col][row];
    }

    size_t getSize() const {
        return this->size;
    }

private:
    int size;
    std::vector<std::vector<double>> matrix;

    void validateIndices(int row, int col) const {
        if (row < 0 || row >= size || col < 0 || col >= size) {
            throw std::out_of_range("Index out of bounds");
        }
    }
};
