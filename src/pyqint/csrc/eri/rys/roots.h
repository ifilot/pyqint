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

#include "fgamma.h"

namespace integrals::eri::rys {

/**
 * @brief Compute the Gaussian quadrature rule for exp(-T t²) on [0, 1].
 *
 * With x = t², the rule integrates polynomials in x through degree 2n - 1.
 *
 * @param n              Number of roots, at least one
 * @param T              Finite, nonnegative Boys argument
 * @param boys_function  Boys function evaluator
 * @param roots          Output array of n nodes in the x = t² convention
 * @param weights        Output array of n weights summing to F_0(T)
 */
void compute_roots(int n, double T, const BoysFunction& boys_function,
                   double* roots, double* weights);

}  // namespace integrals::eri::rys
