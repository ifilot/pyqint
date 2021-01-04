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
 * Incomplete Gamma Function
 *
 * Used in the evaluation of the  two-electron integrals.
 *
 */

/*
 * The functions below were extracted from:
 *
 * Numerical Recipes
 * William H. Press, Saul A. Teukolsky, William T.,
 * Vetterling and Brian P. Flannery
 * 3rd edition page 261
 * ISBN-13: 978-0521880688
 */

#include <iostream>
#include <limits>
#include <cmath>

class GammaInc {
public:
    double Fgamma(const double m, double x) const;

    /*
     * @fn gamm_inc
     * @brief Calculates the incomplete gamma function P(a,x)
     *
     *         gamma(a,x)
     *         ----------
     *            G(a)
     *
     * @param a     "squared width" of the IGF
     * @param x     Upper bound of the integral in the IGF
     *
     * returns double value of the incomplete Gamma Function
     */
    double gamm_inc(const double a, const double x) const;

    /*
     * @fn gamm_inc
     * @brief Calculates the incomplete gamma function P(a,x)
     *
     * This routine selects the best function to use in the
     * evaluation of the Incomplete Gamma Function (IGF).
     *
     * @param a     "squared width" of the IGF
     * @param x     Upper bound of the integral in the IGF
     *
     * returns double value of the incomplete Gamma Function
     */
    double gammp(const double m, const double x) const;

private:
    /*
     * @fn gser
     * @brief Gamma function P(a,x) evaluated by its series representation
     *
     * @param a     "squared width" of the IGF
     * @param x     Upper bound of the integral in the IGF
     *
     * returns double value of the incomplete Gamma Function
     */
    double gser(const double a, const double x) const;


    double gammln(const double xx) const;

    /*
     * @fn gcf
     * @brief Gamma function P(a,x) evaluated by its continued fraction representation
     *
     * @param a     "squared width" of the IGF
     * @param x     Upper bound of the integral in the IGF
     *
     * returns double value of the incomplete Gamma Function
     */
    double gcf(const double a, const double x) const;

    /*
     * @fn gammpapprox
     * @brief Incomplete Gamma function P(a,x) or Q(a,x) evaluated by quadrature
     *
     * Returns P(a,x) or Q(a,x), when psig is 1 or 0, respectively.
     *
     * @param a     "squared width" of the IGF
     * @param x     Upper bound of the integral in the IGF
     * @param psig  Whether to evaluate P(a,x) or Q(a,x)
     *
     * returns double value of the incomplete Gamma Function
     */
    double gammpapprox(double a, double x, int psig) const;
};
