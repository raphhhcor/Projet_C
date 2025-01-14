#include "matrix.h"
#include <cmath>
#include <stdexcept>

// ------------------- Constructors -----------------------
Matrix::Matrix(size_t rows, size_t cols)
    : rows(rows), cols(cols), data(rows, std::vector<double>(cols, 0.0))
{
}

Matrix::Matrix(const std::vector<std::vector<double>>& values)
    : data(values), rows(values.size()), cols(values[0].size())
{
}

// ------------------- Display ----------------
void Matrix::print() const {
    for (const auto& row : data) {
        for (double val : row) {
            std::cout << val << " ";
        }
        std::cout << '\n';
    }
}

// ----------------- Matrices product ---------------------
Matrix Matrix::matrix_product(const Matrix& other) const {
    if (this->cols != other.rows) {
        throw std::invalid_argument("Matrix dimensions do not match for multiplication.");
    }
    Matrix result(this->rows, other.cols);

    for (size_t i = 0; i < this->rows; ++i) {
        for (size_t j = 0; j < other.cols; ++j) {
            double sum = 0.0;
            for (size_t k = 0; k < this->cols; ++k) {
                sum += (*this)(i, k) * other(k, j);
            }
            result(i, j) = sum;
        }
    }
    return result;
}

// ----------------- Tri-diagonal solver (Thomas algorithm) -------------
Matrix Matrix::solveTridiagonal(const Matrix& b) const {
    if (rows != cols) {
        throw std::runtime_error("Matrix must be square for tridiagonal solve.");
    }
    if (b.getRows() != rows || b.getCols() != 1) {
        throw std::runtime_error("Right-hand side vector must be n x 1.");
    }

    size_t n = rows;
    Matrix x(n, 1);  // solution

    // We extract the diags
    std::vector<double> a(n-1, 0.0); 
    std::vector<double> c(n-1, 0.0);
    std::vector<double> d(n,   0.0);
    std::vector<double> rhs(n, 0.0);

    for (size_t i = 0; i < n; ++i) {
        d[i] = data[i][i];
        rhs[i] = b(i, 0);
    }
    for (size_t i = 0; i < n-1; ++i) {
        a[i] = data[i+1][i];
        c[i] = data[i][i+1];
    }

    // 1) Forward elimination
    for (size_t i = 1; i < n; ++i) {
        double w = a[i-1] / d[i-1];
        d[i]   -= w * c[i-1];    
        rhs[i] -= w * rhs[i-1];
    }

    // 2) Back substitution
    x(n-1, 0) = rhs[n-1] / d[n-1];
    for (int i = (int)n - 2; i >= 0; --i) {
        x(i, 0) = (rhs[i] - c[i] * x(i+1, 0)) / d[i];
    }

    return x;
}

// ----------------- removeFirstAndLastRow -----------------
Matrix Matrix::removeFirstAndLastRow() const {
    if (rows <= 2) {
        throw std::logic_error("Matrix must have more than two rows to remove the first and last rows.");
    }

    Matrix result(rows - 2, cols);
    for (size_t i = 1; i < rows - 1; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            result(i - 1, j) = data[i][j];
        }
    }
    return result;
}

// -------------------------------------------------------------------
//  Methods only for the first question but are not used
// -------------------------------------------------------------------
bool Matrix::isInvertible() const {
    if (!isSquare()) return false;
    double det = determinant();
    return (std::fabs(det) > 1e-14);
}

double Matrix::determinant() const {
    if (!isSquare()) {
        throw std::logic_error("Determinant is defined only for square matrices.");
    }
    return determinantRecursive(data);
}

double Matrix::determinantRecursive(const std::vector<std::vector<double>>& mat) const {
    size_t n = mat.size();
    if (n == 1) return mat[0][0];
    if (n == 2) {
        return mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0];
    }
    double det = 0.0;
    for (size_t col = 0; col < n; ++col) {
        double sign = ((col % 2) == 0) ? 1.0 : -1.0;
        det += sign * mat[0][col] * determinantRecursive(getMinor(mat, 0, col));
    }
    return det;
}

std::vector<std::vector<double>> Matrix::getMinor(const std::vector<std::vector<double>>& mat, 
                                                  size_t row, size_t col) const {
    std::vector<std::vector<double>> minor;
    minor.reserve(mat.size() - 1);
    for (size_t i = 0; i < mat.size(); ++i) {
        if (i == row) continue;
        std::vector<double> newRow;
        newRow.reserve(mat[i].size() - 1);
        for (size_t j = 0; j < mat[i].size(); ++j) {
            if (j == col) continue;
            newRow.push_back(mat[i][j]);
        }
        minor.push_back(std::move(newRow));
    }
    return minor;
}

Matrix Matrix::inverse() const {
    if (!isInvertible()) {
        throw std::logic_error("Matrix is not invertible (determinant ~ 0).");
    }

    size_t n = rows;
    Matrix adj(n, n);
    double det = determinant();

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            double sign = ((i + j) % 2 == 0) ? 1.0 : -1.0;
            adj(j, i) = sign * determinantRecursive(getMinor(data, i, j));
        }
    }
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            adj(i, j) /= det;
        }
    }
    return adj;
}
