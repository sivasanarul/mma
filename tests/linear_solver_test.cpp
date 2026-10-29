// SPDX-License-Identifier: GPL-3.0-or-later

#include "linear_solver.h"

#include <cmath>
#include <iostream>
#include <vector>

int main() {
    // The leading zero requires a row swap, exercising partial pivoting.
    const std::vector<double> matrix = {
        0.0, 2.0, 1.0,
        1.0, -2.0, -3.0,
        2.0, 3.0, 1.0,
    };
    const std::vector<double> right_hand_side = {7.0, -12.0, 11.0};
    std::vector<double> solution;
    solve_linear_system(matrix, right_hand_side, solution);

    const std::vector<double> expected = {1.0, 2.0, 3.0};
    for (std::size_t index = 0; index < expected.size(); ++index) {
        if (std::abs(solution[index] - expected[index]) > 1.0e-12) {
            std::cerr << "unexpected solution at index " << index << '\n';
            return 1;
        }
    }
    return 0;
}
