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

#include "rys.h"

#include "roots.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace integrals::eri::rys {

namespace {

struct Axis {
    int a, b, c, d;
    double pa, qc, ab, cd, pq;
};

/**
 * @brief Evaluate the Cartesian polynomial along one axis at a quadrature node.
 *
 * Raise powers about centers A and C using the moments of the coupled
 * Gaussian, then transfer powers to B and D with the horizontal recurrence.
 */
double axis_factor(const Axis& axis, double root, double half_inv_p, double half_inv_q,
                   double q_fraction, double p_fraction, double covariance,
                   double* moments) {
    const int bra = axis.a + axis.b;
    const int ket = axis.c + axis.d;
    if (bra + ket == 0) {
        return 1.0;
    }

    const double bra_shift = q_fraction * root;
    const double ket_shift = p_fraction * root;
    const double mean_a = axis.pa - bra_shift * axis.pq;
    const double mean_c = axis.qc + ket_shift * axis.pq;
    const double variance_a = (1.0 - bra_shift) * half_inv_p;
    const double variance_c = (1.0 - ket_shift) * half_inv_q;
    const int stride = ket + 1;

    moments[0] = 1.0;
    if (ket > 0) {
        moments[1] = mean_c;
        for (int j = 1; j < ket; ++j) {
            moments[j + 1] = mean_c * moments[j] + j * variance_c * moments[j - 1];
        }
    }
    for (int i = 0; i < bra; ++i) {
        double* row = moments + i * stride;
        double* next = row + stride;
        const double adjustment = i * variance_a;
        next[0] = mean_a * row[0];
        if (i > 0) {
            next[0] += adjustment * (row - stride)[0];
        }
        for (int j = 1; j <= ket; ++j) {
            next[j] = mean_a * row[j] + j * covariance * row[j - 1];
            if (i > 0) {
                next[j] += adjustment * (row - stride)[j];
            }
        }
    }

    for (int b = 0; b < axis.b; ++b) {
        for (int i = axis.a; i < bra - b; ++i) {
            double* row = moments + i * stride;
            for (int j = axis.c; j <= ket; ++j) {
                row[j] = row[stride + j] + axis.ab * row[j];
            }
        }
    }
    double* row = moments + axis.a * stride;
    for (int d = 0; d < axis.d; ++d) {
        for (int j = axis.c; j < ket - d; ++j) {
            row[j] = row[j + 1] + axis.cd * row[j];
        }
    }
    return row[axis.c];
}

struct QuartetGeometry {
    Vec3 ab, cd, pa, qc, pq;
    double p, q, inv_pq, p_fraction, q_fraction, T, prefactor;
};

/**
 * @brief Compute the geometry and prefactor shared by a primitive shell quartet.
 */
QuartetGeometry prepare_geometry(const GTO& gto_a, const GTO& gto_b,
                                 const GTO& gto_c, const GTO& gto_d) {
    const Vec3& a = gto_a.get_position();
    const Vec3& b = gto_b.get_position();
    const Vec3& c = gto_c.get_position();
    const Vec3& d = gto_d.get_position();
    const double alpha_a = gto_a.get_alpha();
    const double alpha_b = gto_b.get_alpha();
    const double alpha_c = gto_c.get_alpha();
    const double alpha_d = gto_d.get_alpha();
    const double p = alpha_a + alpha_b;
    const double q = alpha_c + alpha_d;
    const double inv_pq = 1.0 / (p + q);
    const double p_fraction = p * inv_pq;
    const double q_fraction = q * inv_pq;
    const Vec3 ab = a - b;
    const Vec3 cd = c - d;
    const Vec3 pa = (-alpha_b / p) * ab;
    const Vec3 qc = (-alpha_d / q) * cd;
    const Vec3 pq = (a - c) + pa - qc;
    const double T = (p * q_fraction) * pq.norm2();

    constexpr double two_pi_to_five_halves = 34.986836655249725692525643359742;
    const double prefactor = two_pi_to_five_halves / (p * q * std::sqrt(p + q)) *
        std::exp(-alpha_a * (alpha_b / p) * ab.norm2() -
                 alpha_c * (alpha_d / q) * cd.norm2());
    return {ab, cd, pa, qc, pq, p, q, inv_pq, p_fraction, q_fraction, T, prefactor};
}

using AxisPowers = std::array<int, 4>;

Axis make_axis(const QuartetGeometry& geometry, const AxisPowers& powers, int dim) {
    return {powers[0], powers[1], powers[2], powers[3], geometry.pa[dim], geometry.qc[dim],
            geometry.ab[dim], geometry.cd[dim], geometry.pq[dim]};
}

}  // namespace

double evaluate(const GTO& gto_a, const GTO& gto_b,
                const GTO& gto_c, const GTO& gto_d,
                const BoysFunction& boys_function) {
    const auto [ab, cd, pa, qc, pq, p, q, inv_pq, p_fraction, q_fraction, T, prefactor] =
        prepare_geometry(gto_a, gto_b, gto_c, gto_d);
    if (prefactor == 0.0 || std::isinf(T)) {
        return 0.0;
    }

    const std::array<Axis, 3> axes{{
        {static_cast<int>(gto_a.get_l()), static_cast<int>(gto_b.get_l()),
         static_cast<int>(gto_c.get_l()), static_cast<int>(gto_d.get_l()),
         pa[0], qc[0], ab[0], cd[0], pq[0]},
        {static_cast<int>(gto_a.get_m()), static_cast<int>(gto_b.get_m()),
         static_cast<int>(gto_c.get_m()), static_cast<int>(gto_d.get_m()),
         pa[1], qc[1], ab[1], cd[1], pq[1]},
        {static_cast<int>(gto_a.get_n()), static_cast<int>(gto_b.get_n()),
         static_cast<int>(gto_c.get_n()), static_cast<int>(gto_d.get_n()),
         pa[2], qc[2], ab[2], cd[2], pq[2]},
    }};
    int total_order = 0;
    size_t moment_count = 0;
    for (const auto& axis : axes) {
        total_order += axis.a + axis.b + axis.c + axis.d;
        moment_count = std::max(moment_count,
                                static_cast<size_t>(axis.a + axis.b + 1) *
                                static_cast<size_t>(axis.c + axis.d + 1));
    }
    if (total_order == 0) {
        return prefactor * boys_function(0, T);
    }

    const int root_count = total_order / 2 + 1;
    std::array<double, 32> local_rule;
    std::vector<double> large_rule;
    double* roots = local_rule.data();
    if (2 * static_cast<size_t>(root_count) > local_rule.size()) {
        large_rule.resize(2 * static_cast<size_t>(root_count));
        roots = large_rule.data();
    }
    double* weights = roots + root_count;
    compute_roots(root_count, T, boys_function, roots, weights);

    std::array<double, 64> local_moments;
    std::vector<double> large_moments;
    double* moments = local_moments.data();
    if (moment_count > local_moments.size()) {
        large_moments.resize(moment_count);
        moments = large_moments.data();
    }

    const double half_inv_p = 0.5 / p;
    const double half_inv_q = 0.5 / q;
    double sum = 0.0;
    for (int i = 0; i < root_count; ++i) {
        const double root = roots[i];
        const double covariance = 0.5 * inv_pq * root;
        double term = weights[i];
        for (const auto& axis : axes) {
            term *= axis_factor(axis, root, half_inv_p, half_inv_q,
                                q_fraction, p_fraction, covariance, moments);
        }
        sum += term;
    }
    return prefactor * sum;
}

void evaluate_shell_quartet(const ShellComponent* components, std::size_t count,
                            const BoysFunction& boys_function, double* results) {
    // Catalog the distinct angular recurrences needed on each axis.
    std::array<std::vector<AxisPowers>, 3> axis_powers;
    std::vector<std::array<size_t, 3>> component_axes(count);
    size_t moment_count = 0;
    for (size_t i = 0; i < count; ++i) {
        const auto& cgfs = components[i].contractions;
        std::array<AxisPowers, 3> cartesian{};
        for (size_t j = 0; j < cgfs.size(); ++j) {
            const GTO& gto = cgfs[j]->get_gto(0);
            cartesian[0][j] = gto.get_l();
            cartesian[1][j] = gto.get_m();
            cartesian[2][j] = gto.get_n();
        }
        for (size_t dim = 0; dim < axis_powers.size(); ++dim) {
            auto& catalog = axis_powers[dim];
            const auto& powers = cartesian[dim];
            const auto existing = std::find(catalog.begin(), catalog.end(), powers);
            component_axes[i][dim] = existing - catalog.begin();
            if (existing == catalog.end()) {
                catalog.push_back(powers);
                moment_count = std::max(moment_count, static_cast<size_t>(powers[0] + powers[1] + 1) *
                                                     static_cast<size_t>(powers[2] + powers[3] + 1));
            }
        }
    }

    const auto& [a, b, c, d] = components[0].contractions;
    int total_order = 0;
    for (const CGF* cgf : components[0].contractions) {
        const GTO& gto = cgf->get_gto(0);
        total_order += gto.get_l() + gto.get_m() + gto.get_n();
    }
    const int root_count = total_order / 2 + 1;
    std::vector<double> roots(root_count);
    std::vector<double> weights(root_count);
    std::vector<double> moments(moment_count);
    std::array<std::vector<Axis>, 3> axes;
    std::array<std::vector<double>, 3> axis_values;
    for (size_t dim = 0; dim < axes.size(); ++dim) {
        axes[dim].resize(axis_powers[dim].size());
        axis_values[dim].resize(axis_powers[dim].size());
    }
    std::vector<double> primitive_sums(count);
    std::vector<double> contracted_sums(count, 0.0);

    for (size_t i = 0; i < a->size(); ++i) {
        for (size_t j = 0; j < b->size(); ++j) {
            for (size_t k = 0; k < c->size(); ++k) {
                for (size_t l = 0; l < d->size(); ++l) {
                    const QuartetGeometry geometry = prepare_geometry(
                        a->get_gto(i), b->get_gto(j), c->get_gto(k), d->get_gto(l));
                    if (geometry.prefactor == 0.0 || std::isinf(geometry.T)) {
                        continue;
                    }
                    compute_roots(root_count, geometry.T, boys_function, roots.data(), weights.data());
                    for (size_t dim = 0; dim < axes.size(); ++dim) {
                        for (size_t term = 0; term < axes[dim].size(); ++term) {
                            axes[dim][term] = make_axis(geometry, axis_powers[dim][term], dim);
                        }
                    }

                    std::fill(primitive_sums.begin(), primitive_sums.end(), 0.0);
                    const double half_inv_p = 0.5 / geometry.p;
                    const double half_inv_q = 0.5 / geometry.q;
                    for (int r = 0; r < root_count; ++r) {
                        const double root = roots[r];
                        const double covariance = 0.5 * geometry.inv_pq * root;
                        for (size_t dim = 0; dim < axes.size(); ++dim) {
                            for (size_t term = 0; term < axes[dim].size(); ++term) {
                                axis_values[dim][term] = axis_factor(axes[dim][term], root,
                                    half_inv_p, half_inv_q, geometry.q_fraction, geometry.p_fraction,
                                    covariance, moments.data());
                            }
                        }
                        for (size_t m = 0; m < count; ++m) {
                            const auto& [x, y, z] = component_axes[m];
                            primitive_sums[m] += weights[r] * axis_values[0][x] *
                                                 axis_values[1][y] * axis_values[2][z];
                        }
                    }

                    for (size_t m = 0; m < count; ++m) {
                        const auto& [c1, c2, c3, c4] = components[m].contractions;
                        const double coefficients = c1->get_coefficient_gto(i) * c2->get_coefficient_gto(j) *
                                                    c3->get_coefficient_gto(k) * c4->get_coefficient_gto(l);
                        const double norms = c1->get_norm_gto(i) * c2->get_norm_gto(j) *
                                             c3->get_norm_gto(k) * c4->get_norm_gto(l);
                        contracted_sums[m] += coefficients * norms * (geometry.prefactor * primitive_sums[m]);
                    }
                }
            }
        }
    }

    for (size_t m = 0; m < count; ++m) {
        const auto& [c1, c2, c3, c4] = components[m].contractions;
        results[components[m].result_index] = contracted_sums[m] *
            c1->get_contraction_norm() * c2->get_contraction_norm() *
            c3->get_contraction_norm() * c4->get_contraction_norm();
    }
}

}  // namespace integrals::eri::rys
