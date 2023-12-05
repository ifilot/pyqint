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
    Plotter();

    std::vector<double> plot_wavefunction(const std::vector<double>& grid, 
                                          const std::vector<double>& coeff, 
                                          const std::vector<CGF>& cgfs) const;

    std::vector<double> plot_gradient(const std::vector<double>& grid, 
                                      const std::vector<double>& coeff, 
                                      const std::vector<CGF>& cgfs) const;

private:
};
