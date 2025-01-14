#include <iostream>
#include "pde.h"
#include "matrix.h"

// Constructeur
PDE::PDE(const std::vector<double>& time_grid, const std::vector<double>& space_grid,
         const std::vector<double>& min_boundary_condition, const std::vector<double>& max_boundary_condition,
         const std::vector<double>& terminal_condition, const Matrix& a, const Matrix& b, const Matrix& c, const Matrix& d)
    : time_grid(time_grid), space_grid(space_grid), 
      min_boundary_condition(min_boundary_condition), max_boundary_condition(max_boundary_condition),
      terminal_condition(terminal_condition), a(a), b(b), c(c), d(d) {}

// Function for P
Matrix PDE::buildMatrixP(double delta_t, double delta_x, double theta, double coef_i) {
    int m = space_grid.size();
    Matrix P(m - 2, m - 2);

    for (int i = 0; i < m - 2; ++i) {
        P(i, i) = a(coef_i, i + 1) - (1.0 / delta_t + (2.0 * theta * c(coef_i, i + 1)) / (delta_x * delta_x));
        if (i < m - 3) {
            P(i, i + 1) = (b(coef_i, i + 1) / (2.0 * delta_x)) + (theta * c(coef_i, i + 1) / (delta_x * delta_x));
        }
        if (i > 0) {
            P(i, i - 1) = -(b(coef_i, i + 1) / (2.0 * delta_x)) + (theta * c(coef_i, i + 1) / (delta_x * delta_x));
        }
    }

    return P;
}

// Function for Q
Matrix PDE::buildMatrixQ(double delta_t, double delta_x, double theta, double coef_i) const {
    int m = space_grid.size();
    Matrix Q(m - 2, m - 2);

    for (size_t i = 0; i < m - 2; ++i) {
        Q(i, i) = 1 / delta_t - (2 * (1 - theta) * c(coef_i, i + 1)) / (delta_x * delta_x);
        if (i < m - 3) {
            Q(i, i + 1) = ((1 - theta) * c(coef_i, i + 1)) / (delta_x * delta_x);
        }
        if (i > 0) {
            Q(i, i - 1) = ((1 - theta) * c(coef_i, i + 1)) / (delta_x * delta_x);
        }
    }

    return Q;
}

// Function for V
Matrix PDE::buildMatrixV(double delta_t, double delta_x, double theta, double coef_i) const {
    size_t m = d.getCols();
    Matrix V(m - 2, 1);

    for (size_t i = 0; i < m - 2; ++i) {
        double term1 = d(coef_i, i + 1);
        double term2 = 0.0;
        if (i == 0) {
            term2 = (-b(coef_i, i + 1) / (2.0 * delta_x) + theta * c(coef_i, i + 1) / (delta_x * delta_x)) * min_boundary_condition[coef_i] + (1 - theta) * c(coef_i, i + 1) / (delta_x * delta_x) * min_boundary_condition[coef_i + 1];
        }
        if (i == m - 3) {
            term2 = (b(coef_i, i + 1) / (2.0 * delta_x) + theta * c(coef_i, i + 1) / (delta_x * delta_x)) * max_boundary_condition[coef_i] + (1 - theta) * c(coef_i, i + 1) / (delta_x * delta_x) * max_boundary_condition[coef_i + 1];
        }
        V(i, 0) = term1 + term2;
    }
    return V;
}


Matrix PDE::initializeMatrixU() const {
    // Create a matrix with `space_steps` rows and `time_steps` columns
    double time_steps = time_grid.size();
    double space_steps = space_grid.size();
    Matrix U(space_steps, time_steps);

    // Set the last column to the terminal_condition vector
    for (size_t i = 0; i < space_steps; ++i) {
        U(i, time_steps - 1) = terminal_condition[i];
    }
    // Set the first row to the min_boundary_condition vector
    for (size_t i = 0; i < time_steps; ++i) {
        U(0, i) = min_boundary_condition[i];
    }
    // Set the last row to the max_boundary_condition vector
    for (size_t i = 0; i < time_steps; ++i) {
        U(space_steps - 1, i) = max_boundary_condition[i];
    }
    // Return the initialized matrix
    return U;
}

// Compute a specific column of U
void PDE::compute_prec_colu_U(size_t coef_i, Matrix& U, double delta_t, double delta_x, double theta) {
    if ((coef_i >= time_grid.size() - 1) || (coef_i < 0)) {
        throw std::out_of_range("Invalid column index for computation.");
    }
    Matrix P = buildMatrixP(delta_t, delta_x, theta, coef_i);
    Matrix Q = buildMatrixQ(delta_t, delta_x, theta, coef_i);
    Matrix V = buildMatrixV(delta_t, delta_x, theta, coef_i);

    // Retrieve U_next (column coef_i + 1)
    size_t m = U.getRows();
    Matrix U_next(m - 2, 1);
    for (size_t i = 0; i < m - 2; ++i) {
        U_next(i, 0) = U(i + 1, coef_i + 1);
    }

    // Compute Q * U_next
    Matrix Q_U_next = Q.matrix_product(U_next);

    // Add Q * U_next and V
    if (Q_U_next.getRows() != V.getRows() || Q_U_next.getCols() != V.getCols()) {
        throw std::logic_error("Dimensions of Q * U_next and V do not match.");
    }
    Matrix Q_U_next_plus_V(Q_U_next.getRows(), Q_U_next.getCols());
    for (size_t i = 0; i < Q_U_next.getRows(); ++i) {
        for (size_t j = 0; j < Q_U_next.getCols(); ++j) {
            Q_U_next_plus_V(i, j) = Q_U_next(i, j) + V(i, j);
        }
    }
    Matrix U_prec_col = P.solveTridiagonal(Q_U_next_plus_V);

    for (size_t i = 0; i < U_prec_col.getRows(); ++i) {
        U_prec_col(i, 0) = - U_prec_col(i, 0);
    }

    for (size_t i = 0; i < U_prec_col.getRows(); ++i) {
        U(i+1, coef_i) = U_prec_col(i, 0);
    }
}

void PDE::compute_all_columns_U(Matrix& U, double delta_t, double delta_x, double theta) {
    size_t n = time_grid.size();

    if (U.getCols() != n || U.getRows() != space_grid.size()) {
        throw std::logic_error("Matrix U is not correctly initialized.");
    }

    for (int i = n - 2; i >= 0; --i) {
        compute_prec_colu_U(i, U, delta_t, delta_x, theta);
    }
}