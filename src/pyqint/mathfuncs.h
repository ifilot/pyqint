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

/**
 * @brief Integer-based exponentiation using binary exponentiation
 *
 * Computes an integer power of a floating-point base using an efficient
 * binary exponentiation algorithm.
 *
 * Negative exponents are supported via reciprocal evaluation.
 *
 * @param base Base value
 * @param exp Integer exponent
 * @return base raised to the power exp
 */
inline double ipow(double base, int exp) {
    if (exp < 0) {
        return 1.0 / ipow(base, -exp);
    }
    double result = 1.0;
    while (exp > 0) {
        if (exp & 1)
            result *= base;
        base *= base;
        exp >>= 1;
    }
    return result;
}

/**
 * @brief Returns (-1)^n as a double
 *
 * Utility function commonly used for phase factors and parity checks.
 *
 * @param n Integer exponent
 * @return -1.0 if n is odd, +1.0 if n is even
 */
inline double sign_pow(int n) noexcept {
    return (n & 1) ? -1.0 : 1.0;
}