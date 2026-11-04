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

class Plotter {
private:

public:
    /**
     * @brief Construct a Plotter instance.
     */
    Plotter();

    /**
     * @brief Evaluate the wavefunction on a grid.
     *
     * @param grid   Flattened grid coordinates (x, y, z, x, y, z, ...)
     * @param coeff  Basis coefficients
     * @param cgfs   Contracted Gaussian functions
     *
     * @return Wavefunction values per grid point
     */
    std::vector<double> plot_wavefunction(const std::vector<double>& grid, 
                                          const std::vector<double>& coeff, 
                                          const std::vector<CGF>& cgfs) const;

    /**
     * @brief Evaluate the wavefunction gradient on a grid.
     *
     * @param grid   Flattened grid coordinates (x, y, z, x, y, z, ...)
     * @param coeff  Basis coefficients
     * @param cgfs   Contracted Gaussian functions
     *
     * @return Gradient values per grid point (x, y, z, ...)
     */
    std::vector<double> plot_gradient(const std::vector<double>& grid, 
                                      const std::vector<double>& coeff, 
                                      const std::vector<CGF>& cgfs) const;

    /**
     * @brief Evaluate a single basis function on a grid.
     *
     * @param grid  Flattened grid coordinates (x, y, z, x, y, z, ...)
     * @param cgf   Contracted Gaussian function
     *
     * @return Basis function values per grid point
     */
    std::vector<double> plot_basis_function(const std::vector<double>& grid, 
                                            const CGF& cgf) const;
private:
};
