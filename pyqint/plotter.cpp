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

Plotter::Plotter() {}

std::vector<double> Plotter::plot_wavefunction(const std::vector<double>& grid, 
                                               const std::vector<double>& coeff, 
                                               const std::vector<CGF>& cgfs) const {
    std::vector<double> results(grid.size() / 3, 0.0);

    #pragma omp parallel for
    for(int i=0; i<(int)grid.size(); i+=3) { // have to use signed int for MSVC OpenMP here
        for(unsigned int j=0; j<coeff.size(); j++) {
            results[i/3] += coeff[j] * cgfs[j].get_amp(grid[i], grid[i+1], grid[i+2]);
        }
    }

    return results;
}

std::vector<double> Plotter::plot_gradient(const std::vector<double>& grid, 
                                           const std::vector<double>& coeff, 
                                           const std::vector<CGF>& cgfs) const {
    std::vector<double> results(grid.size(), 0.0);

    #pragma omp parallel for
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