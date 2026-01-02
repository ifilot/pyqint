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
 * @brief Boys function evaluator with interpolation and exact fallbacks.
 *
 * Provides:
 *  - exact evaluation for small / large T
 *  - interpolation-based acceleration for moderate T
 *  - efficient block evaluation using stable recurrence
 */
class BoysFunction {
public:
    // Tunable parameters
    static constexpr double T_SMALL = 1e-6;
    static constexpr double T_LARGE = 30.0;
    static constexpr int    NTABLE  = 2048;

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

    /// Maximum order stored in the interpolation table
    int max_order() const noexcept;

private:
    // ---- Interpolation table state ----
    int nu_tab_max_;
    double logT_min_;
    double logT_max_;
    double inv_dlogT_;
    std::vector<double> Htab_;

    // ---- Internal helpers ----
    void init_table(int nu_max);
    double interpolate(int nu, double T) const noexcept;
    int idx(int nu, int i) const noexcept;

    // ---- Exact reference implementations ----
    static double Fgamma_exact(int nu, double T);
    static void   Fgamma_block_exact(int nu_max, double T, double* F);
};