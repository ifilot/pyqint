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

#include "orchestrator.h"

namespace integrals::eri {

namespace {

bool coincident(const Vec3& a, const Vec3& b) {
    return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
}

std::array<unsigned int, 3> powers(const GTO& gto) {
    return {gto.get_l(), gto.get_m(), gto.get_n()};
}

unsigned int angular_momentum(const std::array<unsigned int, 3>& cartesian) {
    return cartesian[0] + cartesian[1] + cartesian[2];
}

bool same_radial_contraction(const CGF& a, const CGF& b) {
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a.get_gto(i).get_alpha() != b.get_gto(i).get_alpha() ||
            a.get_coefficient_gto(i) != b.get_coefficient_gto(i)) {
            return false;
        }
    }
    return true;
}

}  // namespace

EriOrchestrator::EriOrchestrator(const std::vector<CGF>& cgfs) {
    contractions.reserve(cgfs.size());
    for (const auto& cgf : cgfs) {
        Contraction info{cgf.size(), std::nullopt, std::nullopt};
        if (info.primitive_count != 0) {
            const GTO& first = cgf.get_gto(0);
            info.shared_center = first.get_position();
            info.cartesian_powers = powers(first);
            for (std::size_t i = 1; i < info.primitive_count; ++i) {
                const GTO& primitive = cgf.get_gto(i);
                // Fast paths require coincident primitive centers.
                if (!coincident(primitive.get_position(), first.get_position())) {
                    info.shared_center.reset();
                    info.cartesian_powers.reset();
                    break;
                }
                if (powers(primitive) != powers(first)) {
                    info.cartesian_powers.reset();
                }
            }
        }
        contractions.push_back(info);
    }
}

void EriOrchestrator::prepare_shells(const std::vector<CGF>& cgfs) {
    std::vector<std::size_t> representatives;
    for (std::size_t i = 0; i < cgfs.size(); ++i) {
        auto& info = contractions[i];
        info.shell_id = i;
        if (!info.shared_center || !info.cartesian_powers) {
            continue;
        }

        for (const auto j : representatives) {
            const auto& candidate = contractions[j];
            if (info.primitive_count == candidate.primitive_count &&
                angular_momentum(*info.cartesian_powers) == angular_momentum(*candidate.cartesian_powers) &&
                coincident(*info.shared_center, *candidate.shared_center) &&
                same_radial_contraction(cgfs[i], cgfs[j])) {
                info.shell_id = j;
                break;
            }
        }
        if (info.shell_id == i) {
            representatives.push_back(i);
        }
    }
}

void EriOrchestrator::prepare_pairs(const std::vector<CGF>& cgfs) {
    const std::size_t n = contractions.size();
    pairs.resize(n * (n + 1) / 2);
    pair_routing.resize(pairs.size());
    // Bound scratch storage when contractions are unusually large.
    std::size_t remaining = 64 * 1024 * 1024 / sizeof(PrimitivePair);

    for (std::size_t j = 0; j < n; ++j) {
        const auto& b = contractions[j];
        if (!b.shared_center) {
            continue;
        }
        for (std::size_t i = 0; i <= j; ++i) {
            const auto& a = contractions[i];
            const std::size_t index = j * (j + 1) / 2 + i;
            pair_routing[index] = describe_pair(a, b);
            if (!a.shared_center || a.primitive_count > remaining / b.primitive_count) {
                continue;
            }
            auto& pair = pairs[index].emplace();
            pair.displacement = *a.shared_center - *b.shared_center;
            const double distance_squared = pair.displacement.norm2();
            pair.primitives.reserve(a.primitive_count * b.primitive_count);
            for (std::size_t p = 0; p < a.primitive_count; ++p) {
                for (std::size_t q = 0; q < b.primitive_count; ++q) {
                    pair.primitives.push_back(prepare_pair(cgfs[i].get_gto(p), cgfs[j].get_gto(q), distance_squared));
                }
            }
            remaining -= pair.primitives.size();
        }
    }
}

const EriOrchestrator::Pair* EriOrchestrator::pair(std::size_t index) const {
    const auto& pair = pairs[index];
    return pair ? &*pair : nullptr;
}

EriOrchestrator::PairRouting EriOrchestrator::describe_pair(const Contraction& a, const Contraction& b) {
    PairRouting route;
    if (!a.cartesian_powers || !b.cartesian_powers) {
        return route;
    }

    route.center = a.shared_center;
    const unsigned int order_a = angular_momentum(*a.cartesian_powers);
    const unsigned int order_b = angular_momentum(*b.cartesian_powers);
    route.angular_order = order_a + order_b;
    route.all_s_or_p = order_a <= 1 && order_b <= 1;
    route.all_d_or_higher = order_a >= 2 && order_b >= 2;
    for (unsigned int axis = 0; axis < 3; ++axis) {
        const unsigned int bit = 1u << axis;
        const unsigned int power = (*a.cartesian_powers)[axis] + (*b.cartesian_powers)[axis];
        if (power != 0) {
            route.occupied_axes |= bit;
        }
        if (power % 2 == 1) {
            route.odd_axes |= bit;
        }
        if ((*a.shared_center)[axis] == (*b.shared_center)[axis]) {
            route.common_planes |= bit;
        }
    }
    return route;
}

EriOrchestrator::QuartetRoute EriOrchestrator::plan_quartet(std::size_t ij, std::size_t kl) const {
    const auto& bra = pair_routing[ij];
    const auto& ket = pair_routing[kl];
    QuartetRoute route{baseline_engine(EriEngine::Auto)};
    if (!bra.center || !ket.center) {
        return route;
    }

    // Odd pair parities combine with XOR. A zero requires all four centers in the same plane.
    const unsigned int odd_common_planes = (bra.odd_axes ^ ket.odd_axes) & bra.common_planes & ket.common_planes;
    for (unsigned int axis = 0; axis < 3; ++axis) {
        if ((odd_common_planes & (1u << axis)) && (*bra.center)[axis] == (*ket.center)[axis]) {
            route.zero_by_parity = true;
            return route;
        }
    }

    const unsigned int total = bra.angular_order + ket.angular_order;
    const unsigned int axes = bra.occupied_axes | ket.occupied_axes;
    if (total == 0) {
        route.engine = EriEngine::Rys;
    } else if (bra.all_d_or_higher && ket.all_d_or_higher && (axes & (axes - 1)) == 0) {
        route.engine = EriEngine::Rys;
        route.min_shell_components = 2;
    } else if (bra.all_s_or_p && ket.all_s_or_p) {
        // Minimum group sizes validated on planar shell quartets with 1–4 p shells.
        constexpr std::array<std::size_t, 5> measured_components{0, 2, 5, 14, 41};
        route.min_shell_components = measured_components[total];
    }
    return route;
}

}  // namespace integrals::eri
