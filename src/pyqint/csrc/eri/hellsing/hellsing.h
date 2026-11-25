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

#include "cgf.h"
#include "fgamma.h"
#include "hellsing_cache.h"
#include "eri/primitive_pair.h"

namespace integrals::eri::hellsing {

/**
 * @brief Evaluate a primitive ERI by building the Hellsing terms on demand.
 *
 * @param gto_a          First Cartesian Gaussian primitive
 * @param gto_b          Second Cartesian Gaussian primitive
 * @param gto_c          Third Cartesian Gaussian primitive
 * @param gto_d          Fourth Cartesian Gaussian primitive
 * @param boys_function  Boys function evaluator
 *
 * @return Unnormalized electron repulsion integral (ab | cd)
 */
double evaluate_direct(
    const GTO& gto_a, const GTO& gto_b,
    const GTO& gto_c, const GTO& gto_d,
    const BoysFunction& boys_function);

/**
 * @brief Evaluate a primitive ERI using cached one-dimensional Hellsing kernels.
 *
 * @param gto_a          First Cartesian Gaussian primitive
 * @param gto_b          Second Cartesian Gaussian primitive
 * @param gto_c          Third Cartesian Gaussian primitive
 * @param gto_d          Fourth Cartesian Gaussian primitive
 * @param boys_function  Boys function evaluator
 * @param cache          Kernel cache covering all Cartesian powers in the quartet
 *
 * @return Unnormalized electron repulsion integral (ab | cd)
 */
double evaluate_cached(
    const GTO& gto_a, const GTO& gto_b,
    const GTO& gto_c, const GTO& gto_d,
    const BoysFunction& boys_function, const HellsingCacheTable1D& cache);

/**
 * @brief Evaluate a primitive ERI using prepared Gaussian products for both pairs.
 */
double evaluate_prepared(
    const GTO& gto_a, const GTO& gto_b, const GTO& gto_c, const GTO& gto_d,
    const PrimitivePair& bra, const PrimitivePair& ket, const Vec3& ab, const Vec3& cd,
    const BoysFunction& boys_function, const HellsingCacheTable1D& cache);

}  // namespace integrals::eri::hellsing
