#include "functions.h"
#include <cmath>
#include <iostream>

// Standard normal cumulative distribution function
double norm_cdf(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2));
}

// Black-Scholes price calculation with continuous dividend yield
double black_scholes_price(char option_type, 
                           double S,     // Spot price
                           double K,     // Strike price
                           double T,     // Time to maturity
                           double r,     // Risk-free interest rate
                           double sigma, // Volatility
                           double q      // Continuous dividend yield
) {
    // Calculate d1 and d2
    double d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) 
                / (sigma * std::sqrt(T));
    double d2 = d1 - sigma * std::sqrt(T);

    // Call option
    if (option_type == 'C' || option_type == 'c') {
        return S * std::exp(-q * T) * norm_cdf(d1) 
             - K * std::exp(-r * T) * norm_cdf(d2);
    } 
    // Put option
    else if (option_type == 'P' || option_type == 'p') {
        return K * std::exp(-r * T) * norm_cdf(-d2) 
             - S * std::exp(-q * T) * norm_cdf(-d1);
    } 
    // Error handling for invalid option types
    else {
        std::cerr << "Error: Invalid option type. Use 'C' for Call or 'P' for Put." 
                  << std::endl;
        return -1.0; // Return -1 to indicate an error
    }
}

// Function for the time grid
std::vector<double> create_time_grid(double total_time, int num_steps) {
    std::vector<double> time_grid(num_steps + 1);
    double delta_t = total_time / num_steps;

    for (int i = 0; i <= num_steps; ++i) {
        time_grid[i] = i * delta_t;
    }

    return time_grid;
}

// Function of the space grid
std::vector<double> create_space_grid(double spot, double lambda, double sigma, int num_points, double maturity) {
    std::vector<double> space_grid(num_points + 1);
    double delta_x = (2 * lambda * sigma * std::sqrt(maturity)) / num_points;

    for (int i = 0; i <= num_points; ++i) {
        space_grid[i] = spot - lambda * sigma * std::sqrt(maturity) + i * delta_x;
    }

    return space_grid;
}

Matrix computeMatrixA(
    const std::vector<double>& time_grid,
    const std::vector<double>& spacegrid,
    const std::vector<double>& rate_grid
) {
    size_t m = time_grid.size();
    size_t n = spacegrid.size();
    Matrix A(m, n);

    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            A(i, j) = -rate_grid[i];
        }
    }

    return A;
}

Matrix computeMatrixB(
    const std::vector<double>& time_grid,
    const std::vector<double>& spacegrid,
    const std::vector<double>& rate_grid,
    const std::vector<double>& dividend_grid,
    const Matrix& dupireVol
) {
    size_t m = time_grid.size();
    size_t n = spacegrid.size();
    Matrix B(m, n);

    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            B(i, j) = (rate_grid[i] - dividend_grid[i]) * spacegrid[j];
        }
    }

    return B;
}

Matrix computeMatrixC(
    const std::vector<double>& time_grid,
    const std::vector<double>& spacegrid,
    const Matrix& dupireVol
) {
    size_t m = time_grid.size();
    size_t n = spacegrid.size();
    Matrix C(m, n);

    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            double sigma = dupireVol(j, i);
            C(i, j) = 0.5 * sigma * sigma * spacegrid[j] * spacegrid[j];
        }
    }

    return C;
}

Matrix computeMatrixD(
    const std::vector<double>& time_grid,
    const std::vector<double>& spacegrid
) {
    size_t m = time_grid.size();
    size_t n = spacegrid.size();
    Matrix D(m, n);

    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            D(i, j) = 0.0;
        }
    }

    return D;
}

std::vector<double> create_call_terminal_condition(
    const std::vector<double>& space_grid, 
    const std::vector<double>& rate_grid,
    const std::vector<double>& dividend_grid,
    double spot, double strike, double maturity) {
    size_t n = rate_grid.size();
    double delta_t = maturity / (n - 1);

    std::vector<double> terminal_condition(space_grid.size());
    for (size_t i = 0; i < space_grid.size(); ++i) {
        terminal_condition[i] = std::max(0.0, space_grid[i] - strike);
    }

    return terminal_condition;
}

std::vector<double> create_call_min_boundary_condition(
    const std::vector<double>& time_grid,
    const std::vector<double>& space_grid, 
    const std::vector<double>& rate_grid,
    const std::vector<double>& dividend_grid,
    double strike, double maturity) {
    size_t n = rate_grid.size();
    double delta_t = maturity / (n - 1);
    std::vector<double> min_boundary_condition(n);
    for (int i = 0; i < n; i++) {
        min_boundary_condition[i] = std::max(space_grid[0] - strike * exp(-(rate_grid[i] - dividend_grid[i]) * time_grid[i]), 0.0);
    }
    return min_boundary_condition;
}

std::vector<double> create_call_max_boundary_condition(
    const std::vector<double>& time_grid,
    const std::vector<double>& space_grid, 
    const std::vector<double>& rate_grid,
    const std::vector<double>& dividend_grid,
    double strike, double maturity) {
    size_t n = rate_grid.size();
    double delta_t = maturity / (n - 1);
    std::vector<double> max_boundary_condition(n);
    for (int i = 0; i < n; i++) {
        max_boundary_condition[i] = (space_grid[space_grid.size()-1] - strike * exp(-(rate_grid[i] - dividend_grid[i]) * (maturity - time_grid[i])));
    }
    return max_boundary_condition;
}