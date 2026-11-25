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

#include "roots.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace integrals::eri::rys {

namespace {

constexpr int grid_size = 48;

struct GridPoint {
    double x;
    double weight;
};

/**
 * @brief Build the auxiliary Gauss–Legendre grid for t on [0, 1].
 *
 * Store squared nodes so the Stieltjes recurrence operates in x = t².
 */
std::vector<GridPoint> make_grid(int size) {
    std::vector<GridPoint> grid(size);
    for (int i = 0; i < (size + 1) / 2; ++i) {
        double z = std::cos(M_PI * (i + 0.75) / (size + 0.5));
        double derivative;
        for (;;) {
            double polynomial = 1.0;
            double previous = 0.0;
            for (int k = 1; k <= size; ++k) {
                const double next = ((2.0 * k - 1.0) * z * polynomial -
                                     (k - 1.0) * previous) / k;
                previous = polynomial;
                polynomial = next;
            }
            derivative = size * (previous - z * polynomial) / (1.0 - z * z);
            const double next = z - polynomial / derivative;
            if (std::abs(next - z) <= std::numeric_limits<double>::epsilon()) {
                break;
            }
            z = next;
        }
        const double low = (1.0 - z) * 0.5;
        const double high = (1.0 + z) * 0.5;
        const double weight = 1.0 / ((1.0 - z * z) * derivative * derivative);
        grid[i] = {low * low, weight};
        grid[size - 1 - i] = {high * high, weight};
    }
    return grid;
}

/**
 * @brief Build the Jacobi matrix of the weighted grid by Stieltjes recurrence.
 *
 * Carry the square root of the grid weights in the polynomial vectors, so
 * their Euclidean inner products are the weighted polynomial inner products.
 *
 * @return Zeroth moment of the grid
 */
double stieltjes(int n, double T, const std::vector<GridPoint>& grid,
                 double* diagonal, double* off_diagonal) {
    const size_t size = grid.size();
    std::array<double, 2 * grid_size> local;
    std::vector<double> extra;
    double* previous = local.data();
    if (size > grid_size) {
        extra.resize(2 * size);
        previous = extra.data();
    }
    double* current = previous + size;
    std::fill_n(previous, size, 0.0);

    double mass = 0.0;
    for (size_t i = 0; i < size; ++i) {
        const double w = grid[i].weight * std::exp(-T * grid[i].x);
        current[i] = std::sqrt(w);
        mass += w;
    }
    const double normalizer = 1.0 / std::sqrt(mass);
    for (size_t i = 0; i < size; ++i) {
        current[i] *= normalizer;
    }

    double beta = 0.0;
    for (int k = 0; k < n; ++k) {
        double alpha = 0.0;
        for (size_t i = 0; i < size; ++i) {
            alpha += grid[i].x * current[i] * current[i];
        }
        diagonal[k] = alpha;
        if (k + 1 == n) {
            break;
        }

        double norm_squared = 0.0;
        for (size_t i = 0; i < size; ++i) {
            const double next = (grid[i].x - alpha) * current[i] - beta * previous[i];
            previous[i] = next;
            norm_squared += next * next;
        }
        beta = std::sqrt(norm_squared);
        off_diagonal[k] = beta;
        const double inv_beta = 1.0 / beta;
        for (size_t i = 0; i < size; ++i) {
            previous[i] *= inv_beta;
        }
        std::swap(previous, current);
    }
    off_diagonal[n - 1] = 0.0;
    return mass;
}

/**
 * @brief Diagonalize the Jacobi matrix by implicit-shift QL.
 *
 * The Golub–Welsch rule needs only the first component of each eigenvector.
 * Apply the QL rotations to that row and discard the other components.
 */
void diagonalize(int n, double* diagonal, double* off_diagonal, double* first_row) {
    std::fill_n(first_row, n, 0.0);
    first_row[0] = 1.0;
    constexpr double eps = std::numeric_limits<double>::epsilon();

    for (int left = 0; left < n; ++left) {
        for (int iteration = 0;; ++iteration) {
            int right = left;
            while (right + 1 < n &&
                   std::abs(off_diagonal[right]) >
                       eps * (std::abs(diagonal[right]) + std::abs(diagonal[right + 1]))) {
                ++right;
            }
            if (right == left) {
                break;
            }
            if (iteration == 64) {
                throw std::runtime_error("Rys quadrature eigenvalues did not converge");
            }

            const double shift = (diagonal[left + 1] - diagonal[left]) /
                                 (2.0 * off_diagonal[left]);
            const double signed_norm = std::copysign(std::hypot(shift, 1.0), shift);
            double g = diagonal[right] - diagonal[left] +
                       off_diagonal[left] / (shift + signed_norm);
            double sine = 1.0;
            double cosine = 1.0;
            double correction = 0.0;
            int i = right - 1;
            for (; i >= left; --i) {
                const double f = sine * off_diagonal[i];
                const double b = cosine * off_diagonal[i];
                const double norm = std::hypot(f, g);
                off_diagonal[i + 1] = norm;
                if (norm == 0.0) {
                    diagonal[i + 1] -= correction;
                    break;
                }

                sine = f / norm;
                cosine = g / norm;
                g = diagonal[i + 1] - correction;
                const double r = (diagonal[i] - g) * sine + 2.0 * cosine * b;
                correction = sine * r;
                diagonal[i + 1] = g + correction;
                g = cosine * r - b;

                const double next = first_row[i + 1];
                first_row[i + 1] = sine * first_row[i] + cosine * next;
                first_row[i] = cosine * first_row[i] - sine * next;
            }
            if (i < left) {
                diagonal[left] -= correction;
                off_diagonal[left] = g;
            }
            off_diagonal[right] = 0.0;
        }
    }
}

/**
 * @brief Solve the two-point rule using Boys moments through order three.
 */
void two_roots(double T, const BoysFunction& boys_function, double* roots, double* weights) {
    double moments[4];
    boys_function.compute_block(3, T, moments);

    const double mean = moments[1] / moments[0];
    const double second = moments[2] / moments[0];
    const double variance = second - mean * mean;
    const double other = (moments[3] / moments[0] - 2.0 * mean * second +
                          mean * mean * mean) / variance;
    const double gap = std::sqrt((other - mean) * (other - mean) + 4.0 * variance);
    roots[0] = 0.5 * (mean + other - gap);
    roots[1] = 0.5 * (mean + other + gap);
    weights[0] = moments[0] * (roots[1] - mean) / gap;
    weights[1] = moments[0] - weights[0];
}

}  // namespace

void compute_roots(int n, double T, const BoysFunction& boys_function,
                   double* roots, double* weights) {
    if (n == 1) {
        double moments[2];
        boys_function.compute_block(1, T, moments);
        roots[0] = moments[1] / moments[0];
        weights[0] = moments[0];
        return;
    }

    // With y = T x, the measure tends to the generalized Laguerre measure
    // y^(-1/2) exp(-y). The cutoff puts the truncated tail beyond the moments
    // needed for a rule of order n, to double precision for typical ERIs.
    const bool use_laguerre = T >= 4.0 * n + 64.0;
    if (n == 2 && !use_laguerre) {
        two_roots(T, boys_function, roots, weights);
        return;
    }

    std::array<double, 16> local;
    std::vector<double> extra;
    double* off_diagonal = local.data();
    if (n > static_cast<int>(local.size())) {
        extra.resize(n);
        off_diagonal = extra.data();
    }

    double mass;
    if (use_laguerre) {
        mass = 0.5 * std::sqrt(M_PI) / std::sqrt(T);
        for (int i = 0; i < n; ++i) {
            roots[i] = 2.0 * i + 0.5;
            off_diagonal[i] = std::sqrt((i + 1.0) * (i + 0.5));
        }
        off_diagonal[n - 1] = 0.0;
    } else {
        static const auto grid = make_grid(grid_size);
        if (2 * n + 24 <= grid_size) {
            mass = stieltjes(n, T, grid, roots, off_diagonal);
        } else {
            mass = stieltjes(n, T, make_grid(2 * n + 24), roots, off_diagonal);
        }
    }

    diagonalize(n, roots, off_diagonal, weights);
    const double scale = use_laguerre ? 1.0 / T : 1.0;
    for (int i = 0; i < n; ++i) {
        roots[i] *= scale;
        weights[i] = mass * weights[i] * weights[i];
    }
}

}  // namespace integrals::eri::rys
