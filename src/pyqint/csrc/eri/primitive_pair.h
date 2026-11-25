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
#include "gaussian_integrals.h"

namespace integrals::eri {

struct PrimitivePair {
    Vec3 product_center;
    double exponent_sum;
    double attenuation;
};

/**
 * @brief Form a Gaussian product for a primitive pair at a known separation.
 */
inline PrimitivePair prepare_pair(const GTO& a, const GTO& b, double distance_squared) {
    const double alpha_a = a.get_alpha();
    const double alpha_b = b.get_alpha();
    const double exponent_sum = alpha_a + alpha_b;
    return {gaussian::gaussian_product_center(alpha_a, a.get_position(), alpha_b, b.get_position()),
            exponent_sum, std::exp(-alpha_a * alpha_b * distance_squared / exponent_sum)};
}

}  // namespace integrals::eri
