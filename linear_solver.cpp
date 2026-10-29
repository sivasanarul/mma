// SPDX-License-Identifier: GPL-3.0-or-later

#include "linear_solver.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

void solve_linear_system(
    const std::vector<double>& matrix,
    const std::vector<double>& right_hand_side,
    std::vector<double>& solution
) {
    const std::size_t dimension = right_hand_side.size();
    if (dimension == 0 || matrix.size() != dimension * dimension) {
        throw std::invalid_argument("linear system dimensions are inconsistent");
    }

    std::vector<double> coefficients = matrix;
    std::vector<double> values = right_hand_side;
    solution.assign(dimension, 0.0);

    double largest_entry = 0.0;
    for (const double value : coefficients) {
        largest_entry = std::max(largest_entry, std::abs(value));
    }
    const double pivot_tolerance =
        std::numeric_limits<double>::epsilon() *
        std::max(1.0, largest_entry) * static_cast<double>(dimension);

    for (std::size_t column = 0; column < dimension; ++column) {
        std::size_t pivot_row = column;
        double pivot_size = std::abs(coefficients[column * dimension + column]);
        for (std::size_t row = column + 1; row < dimension; ++row) {
            const double candidate = std::abs(coefficients[row * dimension + column]);
            if (candidate > pivot_size) {
                pivot_size = candidate;
                pivot_row = row;
            }
        }
        if (pivot_size <= pivot_tolerance) {
            throw std::runtime_error("linear system is singular to working precision");
        }

        if (pivot_row != column) {
            for (std::size_t entry = column; entry < dimension; ++entry) {
                std::swap(
                    coefficients[column * dimension + entry],
                    coefficients[pivot_row * dimension + entry]
                );
            }
            std::swap(values[column], values[pivot_row]);
        }

        const double pivot = coefficients[column * dimension + column];
        for (std::size_t row = column + 1; row < dimension; ++row) {
            const double factor = coefficients[row * dimension + column] / pivot;
            coefficients[row * dimension + column] = 0.0;
            for (std::size_t entry = column + 1; entry < dimension; ++entry) {
                coefficients[row * dimension + entry] -=
                    factor * coefficients[column * dimension + entry];
            }
            values[row] -= factor * values[column];
        }
    }

    for (std::size_t offset = 0; offset < dimension; ++offset) {
        const std::size_t row = dimension - 1 - offset;
        double value = values[row];
        for (std::size_t column = row + 1; column < dimension; ++column) {
            value -= coefficients[row * dimension + column] * solution[column];
        }
        solution[row] = value / coefficients[row * dimension + row];
    }
}
