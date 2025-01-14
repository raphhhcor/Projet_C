#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <stdexcept>
#include <iostream>

class Matrix {
private:
    std::vector<std::vector<double>> data;
    size_t rows;
    size_t cols;

public:
    // --- Constructors ---
    Matrix(size_t rows, size_t cols);
    Matrix(const std::vector<std::vector<double>>& values);

    // --- Getters ---
    size_t getRows() const { return rows; }
    size_t getCols() const { return cols; }

    // --- To acces to the data ---
    double& operator()(size_t row, size_t col) {
        return data[row][col];
    }
    double operator()(size_t row, size_t col) const {
        return data[row][col];
    }

    // --- Display ---
    void print() const;

    // --- Matrices product ---
    Matrix matrix_product(const Matrix& other) const;

    // --- Tri-diagonal solver (Thomas algorithm) ---
    Matrix solveTridiagonal(const Matrix& b) const;

    // --- Remove the first and the last row of a matrix ---
    Matrix removeFirstAndLastRow() const;

    // ----------------------------------------------------------------------------
    // !! The following methods are here just to compute the inverse of a matrix but are not used in the PDE part !!
    // ----------------------------------------------------------------------------
    bool isSquare() const { return rows == cols; }
    double determinant() const;
    bool isInvertible() const;
    Matrix inverse() const;

private:
    double determinantRecursive(const std::vector<std::vector<double>>& mat) const;
    std::vector<std::vector<double>> getMinor(const std::vector<std::vector<double>>& mat, 
                                              size_t row, size_t col) const;
};

#endif // MATRIX_H
