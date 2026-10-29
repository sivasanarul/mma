// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef MMA_LINEAR_SOLVER_H
#define MMA_LINEAR_SOLVER_H

#include <vector>

void solve_linear_system(
    const std::vector<double>& matrix,
    const std::vector<double>& right_hand_side,
    std::vector<double>& solution
);

#endif
