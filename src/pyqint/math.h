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

/*
 * @fn ipow
 * @brief Fast integer power for small non-negative exponents
 *
 * Computes x^n for integer n using explicit multiplication instead of std::pow.
 * This is intended for small angular momentum exponents (l, m, n) in GTO
 * evaluation, where std::pow(double, int) is prohibitively expensive and
 * inhibits inlining and vectorization.
 *
 * For n = 0..4, the result is fully unrolled. For larger n, a simple
 * multiplication loop is used.
 *
 * @param x   Base value
 * @param n   Non-negative integer exponent
 *
 * @return x raised to the power n
 *
 * @note This function assumes n >= 0. Behavior is undefined for negative n.
 * @note Marked noexcept and inline to enable aggressive optimization.
 */
inline double ipow(double x, int n) noexcept {
    switch (n) {
        case 0: return 1.0;
        case 1: return x;
        case 2: return x * x;
        case 3: return x * x * x;
        case 4: { double x2 = x * x; return x2 * x2; }
        default:
            double r = 1.0;
            for (int i = 0; i < n; ++i) r *= x;
            return r;
    }
}