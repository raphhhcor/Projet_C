#include <iostream>
#include "pde.h"
#include "matrix.h"
#include "functions.h"


Matrix compute_pde_price_main(
    // Input
    int n_timesteps,
    double spot,
    double strike ,
    double maturity,
    double lambda,
    double sigma,
    int m_steps,
    double theta,
    double rate,
    double div
) {
        std::vector<double> rate_grid(n_timesteps + 1, rate);      // Constant rate
    std::vector<double> dividend_grid(n_timesteps + 1, div);  // Div constant

    // Grid init
    std::vector<double> time_grid = create_time_grid(maturity, n_timesteps);
    std::vector<double> space_grid = create_space_grid(spot, lambda, sigma, m_steps, maturity);

    double delta_t = maturity / n_timesteps;
    double delta_x = (2 * lambda * sigma * std::sqrt(maturity)) / m_steps;

    Matrix DupireVol(space_grid.size(), time_grid.size()); // Correctly initialize with dimensions

    // Fill DupireVol with appropriate values
    for (size_t i = 0; i < space_grid.size(); ++i) {
        for (size_t j = 0; j < time_grid.size(); ++j) {
            // Compute the volatility based on space_grid and time_grid
            DupireVol(i, j) = sigma;
        }
    }
    // Calcul des matrices a, b, c, d
    Matrix A = computeMatrixA(time_grid, space_grid, rate_grid);
    Matrix B = computeMatrixB(time_grid, space_grid, rate_grid, dividend_grid, DupireVol);
    Matrix C = computeMatrixC(time_grid, space_grid, DupireVol);
    Matrix D = computeMatrixD(time_grid, space_grid);

    // Condition terminale
    std::vector<double> terminal_condition = create_call_terminal_condition(space_grid, rate_grid, dividend_grid, spot, strike, maturity);

    std::vector<double> min_boundary = create_call_min_boundary_condition(time_grid, space_grid, rate_grid, dividend_grid, strike, maturity);

    std::vector<double> max_boundary = create_call_max_boundary_condition(time_grid, space_grid, rate_grid, dividend_grid, strike, maturity);

    // Création d'une instance de PDE
    PDE pde(time_grid, space_grid, min_boundary, max_boundary, terminal_condition, A, B, C, D);


    // Construire la matrice U
    Matrix U = pde.initializeMatrixU();

    pde.compute_all_columns_U(U, delta_t, delta_x, theta);

    return U;
}

int main() {
    // Input
    int n_timesteps = 400;
    double spot = 100.0;
    double strike = 100.0;
    double maturity = 1.0;
    double lambda = 130.0;
    double sigma = 0.2;
    int m_steps = 400;
    double theta = 0.5;
    double rate = 0.05;
    double div = 0.02;

    Matrix U = compute_pde_price_main(
        n_timesteps,
        spot,
        strike ,
        maturity,
        lambda,
        sigma,
        m_steps,
        theta,
        rate,
        div
    );

    std::cout << "PDE - ATM Call Price = " << U(ceil(m_steps / 2), 0) << std::endl;

    double call_price = black_scholes_price('C', spot, strike, maturity, rate, sigma, div);
    std::cout << "BS - ATM Call Price = " << call_price << std::endl;

    // Parameters for grid refinement
    std::vector<int> timeSteps = {50, 100, 200, 400}; // Different N_t values
    std::vector<int> spaceSteps = {50, 100, 200, 400}; // Different M values

    // Print header
    std::cout << "TimeSteps\tSpaceSteps\tComputedPrice\tAnalyticalPrice\tError\n";

    for (int M : spaceSteps) {
        for (int N_t : timeSteps) {
            Matrix U = compute_pde_price_main(
                N_t,
                spot,
                strike ,
                maturity,
                lambda,
                sigma,
                M,
                theta,
                rate,
                div
            );
            double computedPrice = U(ceil(M / 2), 0);

            // Compute error
            double error = std::abs(computedPrice - call_price);

            // Print results to terminal
            std::cout << N_t << "\t\t" << M << "\t\t" << computedPrice
                      << "\t\t" << call_price << "\t\t" << error << std::endl;
        }
    }

    std::cout << "Convergence analysis completed.\n";

    return 0;
}