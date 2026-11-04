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

#include "plotter.h"

/**
 * @brief Construct a Plotter instance.
 */
Plotter::Plotter() {}

/**
 * @brief Evaluate the wavefunction on a grid.
 *
 * @param grid   Flattened grid coordinates (x, y, z, x, y, z, ...)
 * @param coeff  Basis coefficients
 * @param cgfs   Contracted Gaussian functions
 *
 * @return Wavefunction values per grid point
 */
std::vector<double> Plotter::plot_wavefunction(const std::vector<double>& grid, 
                                               const std::vector<double>& coeff, 
                                               const std::vector<CGF>& cgfs) const {
    std::vector<double> results(grid.size() / 3, 0.0);

    #pragma omp parallel for schedule(dynamic, 64)
    for(int i=0; i<(int)grid.size(); i+=3) { // have to use signed int for MSVC OpenMP here
        for(unsigned int j=0; j<coeff.size(); j++) {
            results[i/3] += coeff[j] * cgfs[j].get_amp(grid[i], grid[i+1], grid[i+2]);
        }
    }

    return results;
}

/**
 * @brief Evaluate a single basis function on a grid.
 *
 * @param grid  Flattened grid coordinates (x, y, z, x, y, z, ...)
 * @param cgf   Contracted Gaussian function
 *
 * @return Basis function values per grid point
 */
std::vector<double> Plotter::plot_basis_function(const std::vector<double>& grid, 
                                                 const CGF& cgf) const {
    std::vector<double> results(grid.size() / 3, 0.0);

    #pragma omp parallel for schedule(dynamic, 64)
    for(int i=0; i<(int)grid.size(); i+=3) { // have to use signed int for MSVC OpenMP here
        results[i/3] = cgf.get_amp(grid[i], grid[i+1], grid[i+2]);
    }

    return results;
}

/**
 * @brief Evaluate a set of basis functions on a grid.
 *
 * All basis functions are evaluated in a single parallel loop over the grid
 * points. Each thread evaluates every basis function for its share of the
 * grid points, such that the cost of starting the parallel region is paid
 * only once instead of once per basis function.
 *
 * @param grid  Flattened grid coordinates (x, y, z, x, y, z, ...)
 * @param cgfs  Contracted Gaussian functions
 *
 * @return Basis function values, ordered as [basis function][grid point]
 */
std::vector<double> Plotter::plot_basis_functions(const std::vector<double>& grid,
                                                  const std::vector<CGF>& cgfs) const {
    const size_t npts = grid.size() / 3;
    std::vector<double> results(cgfs.size() * npts, 0.0);

    #pragma omp parallel for schedule(static)
    for(int i=0; i<(int)npts; i++) { // have to use signed int for MSVC OpenMP here
        for(size_t j=0; j<cgfs.size(); j++) {
            results[j * npts + i] = cgfs[j].get_amp(grid[3*i], grid[3*i+1], grid[3*i+2]);
        }
    }

    return results;
}

/**
 * @brief Evaluate the gradients of a set of basis functions on a grid.
 *
 * @param grid  Flattened grid coordinates (x, y, z, x, y, z, ...)
 * @param cgfs  Contracted Gaussian functions
 *
 * @return Gradient values, ordered as [basis function][grid point][x, y, z]
 */
std::vector<double> Plotter::plot_basis_function_gradients(const std::vector<double>& grid,
                                                           const std::vector<CGF>& cgfs) const {
    const size_t npts = grid.size() / 3;
    std::vector<double> results(cgfs.size() * npts * 3, 0.0);

    #pragma omp parallel for schedule(static)
    for(int i=0; i<(int)npts; i++) { // have to use signed int for MSVC OpenMP here
        for(size_t j=0; j<cgfs.size(); j++) {
            const auto grad = cgfs[j].get_grad(grid[3*i], grid[3*i+1], grid[3*i+2]);
            for(unsigned int k=0; k<3; k++) {
                results[(j * npts + i) * 3 + k] = grad[k];
            }
        }
    }

    return results;
}

/**
 * @brief Evaluate the wavefunction gradient on a grid.
 *
 * @param grid   Flattened grid coordinates (x, y, z, x, y, z, ...)
 * @param coeff  Basis coefficients
 * @param cgfs   Contracted Gaussian functions
 *
 * @return Gradient values per grid point (x, y, z, ...)
 */
std::vector<double> Plotter::plot_gradient(const std::vector<double>& grid, 
                                           const std::vector<double>& coeff, 
                                           const std::vector<CGF>& cgfs) const {
    std::vector<double> results(grid.size(), 0.0);

    #pragma omp parallel for schedule(dynamic, 64)
    for(int i=0; i<(int)grid.size(); i+=3) { // have to use signed int for MSVC OpenMP here
        for(unsigned int j=0; j<coeff.size(); j++) {
            const auto grad = cgfs[j].get_grad(grid[i], grid[i+1], grid[i+2]);
            for(unsigned int k=0; k<3; k++) {
                results[i+k] += coeff[j] * grad[k];
            }
        }
    }

    return results;
}
