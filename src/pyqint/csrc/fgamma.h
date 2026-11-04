/**************************************************************************
 *   This file is part of PyQInt.                                         *
 *                                                                        *
 *   Author: Ivo Filot <ivo@ivofilot.nl>                                  *
 *                                                                        *
 *   PyQInt is free software:                                             *
 *   you can redistribute it and/or modify it under the terms of the      *
 *   GNU General Public License as published by the Free Software         *
 *   Foundation, either version 3 of the License, or (at your option)     *
 *   any later version.                                                   *
 *                                                                        *
 *   PyQInt is distributed in the hope that it will be useful,            *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty          *
 *   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.              *
 *   See the GNU General Public License for more details.                 *
 *                                                                        *
 *   You should have received a copy of the GNU General Public License    *
 *   along with this program.  If not, see http://www.gnu.org/licenses/.  *
 *                                                                        *
 **************************************************************************/

#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include <vector>
#include <limits>
#include <algorithm>

/**
 * @brief Boys function evaluator based on a tabulated Taylor expansion.
 *
 * Provides:
 *  - F_ν(T) to near machine precision for all orders and arguments
 *  - a Taylor expansion around tabulated points for T < T_TABLE
 *  - erf-based upward recurrence for large T, where it is stable
 *  - efficient block evaluation using stable recurrence
 */
class BoysFunction {
public:
    // Tunable parameters
    static constexpr double T_TABLE  = 40.0;    // tabulated range [0, T_TABLE)
    static constexpr double DT_TABLE = 0.05;    // spacing of the tabulated points
    static constexpr int    NTAYLOR  = 8;       // terms of the Taylor expansion

    /**
     * @brief Construct a Boys-function evaluator.
     *
     * @param nu_max  Maximum Boys order to precompute.
     */
    explicit BoysFunction(int nu_max);

    /**
     * @brief Evaluate the Boys function F_ν(T).
     *
     * @param nu  Boys order (nu ≥ 0)
     * @param T   Argument (T ≥ 0)
     */
    double operator()(int nu, double T) const;

    /**
     * @brief Compute F_0(T) … F_νmax(T).
     *
     * @param nu_max  Maximum order
     * @param T       Argument
     * @param F       Output array (size ≥ nu_max+1)
     */
    void compute_block(int nu_max, double T, double* F) const;

    /**
     * @brief Return the maximum order stored in the table.
     *
     * @return maximum Boys order
     */
    int max_order() const noexcept;

private:
    // ---- Table state ----
    int nu_tab_max_;            // highest order evaluated from the table
    int nrows_;                 // orders stored: 0 … nu_tab_max_ + NTAYLOR - 1
    int npoints_;               // tabulated points T_i = i * DT_TABLE
    std::vector<double> Ftab_;  // F_ν(T_i), stored per point

    // ---- Internal helpers ----
    /**
     * @brief Initialize the table.
     *
     * @param nu_max  Maximum order to precompute
     */
    void init_table(int nu_max);

    /**
     * @brief Evaluate F_ν(T) by a Taylor expansion around the nearest
     *        tabulated point.
     *
     * @param nu  Boys order (≤ nu_tab_max_)
     * @param T   Argument (< T_TABLE)
     *
     * @return Boys function value
     */
    double taylor(int nu, double T) const noexcept;

    // ---- Exact reference implementation ----
    /**
     * @brief Evaluate the Boys function without the table.
     *
     * @param nu  Boys order
     * @param T   Argument
     *
     * @return exact value
     */
    static double Fgamma_exact(int nu, double T);
};
