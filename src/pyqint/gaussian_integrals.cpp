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

#define _USE_MATH_DEFINES
#include <cmath>

#include "factorials.h"
#include "gaussian_integrals.h"

namespace integrals::gaussian {

/**
 * @brief Computes the Gaussian product center of two primitive Gaussians
 *
 * @param alpha1  Exponent of first Gaussian
 * @param a       Center of first Gaussian
 * @param alpha2  Exponent of second Gaussian
 * @param b       Center of second Gaussian
 *
 * Calculates the position of the Gaussian product center P.
 *
 * @return Vec3 Position vector of the Gaussian product center
 */
Vec3 gaussian_product_center(double alpha1, const Vec3& a,
                             double alpha2, const Vec3& b) {
    return (alpha1 * a + alpha2 * b) / (alpha1 + alpha2);
}

/**
 * @brief Computes the binomial coefficient
 *
 * @param a  Upper index
 * @param b  Lower index
 *
 * Returns 1 if indices are outside valid range.
 *
 * @return double Value of the binomial coefficient
 */
double binomial(int a, int b) {
    if ((a < 0) || (b < 0) || (a - b < 0)) {
        return 1.0;
    }
    return factorial(a) / (factorial(b) * factorial(a - b));
}

/**
 * @brief Computes binomial prefactor used in Gaussian integral recursion
 *
 * @param s    Summation index
 * @param ia   Angular momentum component of first Gaussian
 * @param ib   Angular momentum component of second Gaussian
 * @param xpa  Distance P_x - A_x
 * @param xpb  Distance P_x - B_x
 *
 * @return double Value of the binomial prefactor
 */
double binomial_prefactor(int s, int ia, int ib,
                          double xpa, double xpb) {
    double sum = 0.0;

    for (int t = 0; t < s + 1; t++) {
        if ((s - ia <= t) && (t <= ib)) {
            sum += binomial(ia, s - t) *
                   binomial(ib, t) *
                   std::pow(xpa, ia - s + t) *
                   std::pow(xpb, ib - t);
        }
    }

    return sum;
}

/**
 * @brief Computes one-dimensional overlap integral component
 *
 * @param l1     Angular momentum of first Gaussian along axis
 * @param l2     Angular momentum of second Gaussian along axis
 * @param x1     Distance P_x - A_x
 * @param x2     Distance P_x - B_x
 * @param gamma  Sum of Gaussian exponents
 *
 * @return double Value of the 1D overlap contribution
 */
double overlap_1d(int l1, int l2, double x1, double x2, double gamma) {
    double sum = 0.0;

    for (int i = 0; i < (1 + std::floor(0.5 * (l1 + l2))); i++) {
        sum += binomial_prefactor(2 * i, l1, l2, x1, x2) *
               (i == 0 ? 1 : double_factorial(2 * i - 1)) /
               std::pow(2 * gamma, i);
    }

    return sum;
}

/**
 * @brief Computes overlap integral between two primitive Gaussian orbitals
 *
 * @param alpha1  Exponent of first Gaussian
 * @param l1      Angular momentum in x for first Gaussian
 * @param m1      Angular momentum in y for first Gaussian
 * @param n1      Angular momentum in z for first Gaussian
 * @param a       Center of first Gaussian
 * @param alpha2  Exponent of second Gaussian
 * @param l2      Angular momentum in x for second Gaussian
 * @param m2      Angular momentum in y for second Gaussian
 * @param n2      Angular momentum in z for second Gaussian
 * @param b       Center of second Gaussian
 *
 * @return double Value of the 3D overlap integral
 */
double overlap_gto(double alpha1, unsigned int l1, unsigned int m1, unsigned int n1, const Vec3& a,
                   double alpha2, unsigned int l2, unsigned int m2, unsigned int n2, const Vec3& b) {

    double rab2 = (a - b).norm2();
    double gamma = alpha1 + alpha2;
    Vec3 p = gaussian_product_center(alpha1, a, alpha2, b);

    double pre = std::pow(M_PI / gamma, 1.5) * std::exp(-alpha1 * alpha2 * rab2 / gamma);
    double wx = overlap_1d(l1, l2, p[0] - a[0], p[0] - b[0], gamma);
    double wy = overlap_1d(m1, m2, p[1] - a[1], p[1] - b[1], gamma);
    double wz = overlap_1d(n1, n2, p[2] - a[2], p[2] - b[2], gamma);

    return pre * wx * wy * wz;
}

}  // namespace normalization
