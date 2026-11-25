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

#include <array>
#include <cstddef>

#include "cgf.h"
#include "fgamma.h"

namespace integrals::eri::rys {

/**
 * @brief Evaluate a primitive ERI using Rys quadrature.
 *
 * @param gto_a          First Cartesian Gaussian primitive
 * @param gto_b          Second Cartesian Gaussian primitive
 * @param gto_c          Third Cartesian Gaussian primitive
 * @param gto_d          Fourth Cartesian Gaussian primitive
 * @param boys_function  Boys function evaluator
 *
 * @return Unnormalized electron repulsion integral (ab | cd)
 */
double evaluate(const GTO& gto_a, const GTO& gto_b,
                const GTO& gto_c, const GTO& gto_d,
                const BoysFunction& boys_function);

struct ShellComponent {
    std::array<const CGF*, 4> contractions;
    std::size_t result_index;
};

/**
 * @brief Evaluate Cartesian components of one verified shell quartet together.
 *
 * @param components     Quartets sharing centers, radial contractions, and angular orders
 * @param count          Number of Cartesian components
 * @param boys_function  Boys function evaluator
 * @param results        Packed tensor to receive the contracted integrals
 */
void evaluate_shell_quartet(const ShellComponent* components, std::size_t count,
                            const BoysFunction& boys_function, double* results);

}  // namespace integrals::eri::rys
