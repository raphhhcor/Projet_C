#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include "matrix.h"
#include <vector>

// Standard normal cumulative distribution function
double norm_cdf(double x);

// Black-Scholes price calculation
double black_scholes_price(char option_type, double S, double K, double T, double r, double sigma, double q);

// Function for the time grid
std::vector<double> create_time_grid(double total_time, int num_steps);

// Function of the space grid
std::vector<double> create_space_grid(double spot, double lambda, double sigma, int num_points, double maturity);

Matrix computeMatrixA(
    const std::vector<double>& time_grid,
    const std::vector<double>& spacegrid,
    const std::vector<double>& rate_grid
);

Matrix computeMatrixB(
    const std::vector<double>& time_grid,
    const std::vector<double>& spacegrid,
    const std::vector<double>& rate_grid,
    const std::vector<double>& dividend_grid,
    const Matrix& dupireVol
);

Matrix computeMatrixC(
    const std::vector<double>& time_grid,
    const std::vector<double>& spacegrid,
    const Matrix& dupireVol
);

Matrix computeMatrixD(
    const std::vector<double>& time_grid,
    const std::vector<double>& spacegrid
);

std::vector<double> create_call_terminal_condition(
    const std::vector<double>& space_grid, 
    const std::vector<double>& rate_grid,
    const std::vector<double>& dividend_grid,
    double spot, double strike, double maturity
);

std::vector<double> create_call_min_boundary_condition(
    const std::vector<double>& time_grid,
    const std::vector<double>& space_grid, 
    const std::vector<double>& rate_grid,
    const std::vector<double>& dividend_grid,
    double strike, double maturity
);

std::vector<double> create_call_max_boundary_condition(
    const std::vector<double>& time_grid,
    const std::vector<double>& space_grid, 
    const std::vector<double>& rate_grid,
    const std::vector<double>& dividend_grid,
    double strike, double maturity
);

#endif // FUNCTIONS_H
