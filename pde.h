#include "matrix.h"

class PDE {
public:
    // Constructors
    PDE(const std::vector<double>& time_grid, const std::vector<double>& space_grid,
        const std::vector<double>& min_boundary_condition, const std::vector<double>& max_boundary_condition,
        const std::vector<double>& terminal_condition, const Matrix& a, const Matrix& b, const Matrix& c, const Matrix& d);

    // Function for P
    Matrix buildMatrixP(double delta_t, double delta_x, double theta, double coef_i);
    // Function for Q
    Matrix buildMatrixQ(double delta_t, double delta_x, double theta, double coef_i) const;
    // Function for V
    Matrix buildMatrixV(double delta_t, double delta_x, double theta, double coef_i) const;
    
    Matrix initializeMatrixU() const;

    // Function to compute a specific column of U
    void compute_prec_colu_U(size_t coef_i, Matrix& U, double delta_t, double delta_x, double theta);

    void compute_all_columns_U(Matrix& U, double delta_t, double delta_x, double theta);
private:
    std::vector<double> time_grid;
    std::vector<double> space_grid;
    std::vector<double> min_boundary_condition;
    std::vector<double> max_boundary_condition;
    std::vector<double> terminal_condition;
    Matrix a, b, c, d;
};
