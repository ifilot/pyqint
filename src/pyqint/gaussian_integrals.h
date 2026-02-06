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

#include "vec3.h"

namespace integrals::gaussian {

/**
 * @brief Computes the Gaussian product center of two primitive Gaussians
 *
 * @param double alpha1  Exponent of first Gaussian
 * @param const Vec3& a  Center of first Gaussian
 * @param double alpha2  Exponent of second Gaussian
 * @param const Vec3& b  Center of second Gaussian
 *
 * Calculates the position of the Gaussian product center P given by:
 * P = (alpha1*A + alpha2*B) / (alpha1 + alpha2)
 *
 * @return Vec3 Position vector of the Gaussian product center
 */
Vec3 gaussian_product_center(double alpha1, const Vec3& a,
                              double alpha2, const Vec3& b);

/**
 * @brief Computes the binomial coefficient
 *
 * @param int a  Upper index
 * @param int b  Lower index
 *
 * Calculates the binomial coefficient:
 * (a choose b) = a! / (b! (a-b)!)
 *
 * @return double Value of the binomial coefficient
 */
double binomial(int a, int b);

/**
 * @brief Computes binomial prefactor used in Gaussian integral recursion
 *
 * @param int s      Summation index
 * @param int ia     Angular momentum component of first Gaussian
 * @param int ib     Angular momentum component of second Gaussian
 * @param double xpa Distance P_x - A_x
 * @param double xpb Distance P_x - B_x
 *
 * Evaluates the binomial prefactor appearing in Hermite Gaussian
 * expansion formulas for overlap integrals.
 *
 * @return double Value of the binomial prefactor
 */
double binomial_prefactor(int s, int ia, int ib, double xpa, double xpb);

/**
 * @brief Computes one-dimensional overlap integral component
 *
 * @param int l1     Angular momentum of first Gaussian along axis
 * @param int l2     Angular momentum of second Gaussian along axis
 * @param double x1  Distance P_x - A_x
 * @param double x2  Distance P_x - B_x
 * @param double gamma Sum of Gaussian exponents (alpha1 + alpha2)
 *
 * Computes the 1D overlap term used in separable Gaussian integrals.
 *
 * @return double Value of the 1D overlap contribution
 */
double overlap_1d(int l1, int l2, double x1, double x2, double gamma);

/**
 * @brief Computes overlap integral between two primitive Gaussian orbitals
 *
 * @param double alpha1  Exponent of first Gaussian
 * @param unsigned int l1 Angular momentum in x for first Gaussian
 * @param unsigned int m1 Angular momentum in y for first Gaussian
 * @param unsigned int n1 Angular momentum in z for first Gaussian
 * @param const Vec3& a  Center of first Gaussian
 * @param double alpha2  Exponent of second Gaussian
 * @param unsigned int l2 Angular momentum in x for second Gaussian
 * @param unsigned int m2 Angular momentum in y for second Gaussian
 * @param unsigned int n2 Angular momentum in z for second Gaussian
 * @param const Vec3& b  Center of second Gaussian
 *
 * Calculates the value of < gto1 | gto2 > for primitive Gaussian type orbitals.
 *
 * @return double Value of the 3D overlap integral
 */
double overlap_gto(double alpha1, unsigned int l1, unsigned int m1, unsigned int n1, const Vec3& a,
                   double alpha2, unsigned int l2, unsigned int m2, unsigned int n2, const Vec3& b);

}  // namespace normalization
