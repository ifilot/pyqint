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

#include "integrals.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

namespace {

struct DensityPair {
    std::array<size_t, 2> basis;
    size_t orientations;
    double coulomb_weight;
};

}  // namespace

void Integrator::contract_packed_tei(size_t n, const double* values, const double* density_a,
                                     const double* density_b, double* coulomb, double* exchange_a,
                                     double* exchange_b) const {
    if (n == 0) {
        return;
    }

    const size_t matrix_size = n * n;
    const size_t pair_count = n * (n + 1) / 2;
    const size_t quartet_count = pair_count * (pair_count + 1) / 2;
    std::vector<DensityPair> pairs;
    pairs.reserve(pair_count);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            const size_t ij = i * n + j;
            const size_t ji = j * n + i;
            double weight = density_a[ij];
            if (density_b) {
                weight += density_b[ij];
            }
            if (i != j) {
                weight += density_a[ji];
                if (density_b) {
                    weight += density_b[ji];
                }
            }
            pairs.push_back({{i, j}, i == j ? 1u : 2u, weight});
        }
    }

    int thread_count = 1;
    #ifdef _OPENMP
    const size_t useful_threads = std::max<size_t>(1, quartet_count / 50000);
    thread_count = static_cast<int>(std::min<size_t>(omp_get_max_threads(), useful_threads));
    #endif
    const size_t workspace_size = pair_count + (density_b ? 2 : 1) * matrix_size;
    std::vector<double> partials(thread_count * workspace_size, 0.0);

    #pragma omp parallel num_threads(thread_count) if(thread_count > 1)
    {
        int thread = 0;
        #ifdef _OPENMP
        thread = omp_get_thread_num();
        #endif
        double* j_part = partials.data() + thread * workspace_size;
        double* ka_part = j_part + pair_count;
        double* kb_part = ka_part + matrix_size;

        const auto row_count = static_cast<std::ptrdiff_t>(pair_count);
        #pragma omp for schedule(static, 16)
        for (std::ptrdiff_t row = 0; row < row_count; ++row) {
            const auto p = static_cast<size_t>(row);
            const auto& bra = pairs[p];
            const double* row_values = values + p * (p + 1) / 2;
            for (size_t q = 0; q <= p; ++q) {
                const double value = row_values[q];
                if (value == 0.0) {
                    continue;
                }
                const auto& ket = pairs[q];
                const bool swap_pairs = p != q;
                j_part[p] += value * ket.coulomb_weight;
                if (swap_pairs) {
                    j_part[q] += value * bra.coulomb_weight;
                }

                // Include each distinct within-pair orientation and its pair swap.
                for (size_t a = 0; a < bra.orientations; ++a) {
                    const size_t i = bra.basis[a];
                    const size_t j = bra.basis[1 - a];
                    for (size_t b = 0; b < ket.orientations; ++b) {
                        const size_t k = ket.basis[b];
                        const size_t l = ket.basis[1 - b];
                        ka_part[i * n + k] += value * density_a[j * n + l];
                        if (density_b) {
                            kb_part[i * n + k] += value * density_b[j * n + l];
                        }
                        if (swap_pairs) {
                            ka_part[k * n + i] += value * density_a[l * n + j];
                            if (density_b) {
                                kb_part[k * n + i] += value * density_b[l * n + j];
                            }
                        }
                    }
                }
            }
        }
    }

    std::fill_n(coulomb, matrix_size, 0.0);
    std::fill_n(exchange_a, matrix_size, 0.0);
    if (density_b) {
        std::fill_n(exchange_b, matrix_size, 0.0);
    }
    for (int thread = 0; thread < thread_count; ++thread) {
        const double* j_part = partials.data() + thread * workspace_size;
        const double* ka_part = j_part + pair_count;
        const double* kb_part = ka_part + matrix_size;
        for (size_t p = 0; p < pair_count; ++p) {
            const auto& [i, j] = pairs[p].basis;
            coulomb[i * n + j] += j_part[p];
        }
        for (size_t i = 0; i < matrix_size; ++i) {
            exchange_a[i] += ka_part[i];
        }
        if (density_b) {
            for (size_t i = 0; i < matrix_size; ++i) {
                exchange_b[i] += kb_part[i];
            }
        }
    }
    for (const auto& pair : pairs) {
        const auto& [i, j] = pair.basis;
        coulomb[j * n + i] = coulomb[i * n + j];
    }
}
