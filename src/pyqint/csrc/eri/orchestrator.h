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

#include <array>
#include <cstddef>
#include <optional>
#include <vector>

#include "primitive_pair.h"

enum class EriEngine {
    Hellsing,
    HellsingDirect,
    Huzinaga,
    HuzinagaScalar,
    Rys,
    Auto
};

namespace integrals::eri {

/** @brief Timings and route counts for an ERI tensor build. */
struct EriBuildProfile {
    double pair_setup = 0.0;
    double quartet_planning = 0.0;
    double integral_evaluation = 0.0;
    double tensor_assembly = 0.0;
    std::size_t parity_zero_quartets = 0;
    std::size_t prepared_hellsing_quartets = 0;
    std::size_t rys_quartets = 0;
    std::size_t rys_shell_batches = 0;
};

/**
 * @brief Prepare the contraction metadata used by automatic ERI dispatch.
 */
class EriOrchestrator {
public:
    struct Pair {
        Vec3 displacement;
        std::vector<PrimitivePair> primitives;
    };

    struct Contraction {
        std::size_t primitive_count;
        std::optional<Vec3> shared_center;
        std::optional<std::array<unsigned int, 3>> cartesian_powers;
        std::size_t shell_id = 0;
    };

    struct QuartetRoute {
        EriEngine engine;
        std::size_t min_shell_components = 0;
        bool zero_by_parity = false;
    };

    /**
     * @brief Inspect the primitive centers of every contraction once.
     *
     * @param cgfs  Contracted functions in basis order
     */
    explicit EriOrchestrator(const std::vector<CGF>& cgfs);

    /**
     * @brief Prepare primitive pairs for single-center contractions in basis order.
     */
    void prepare_pairs(const std::vector<CGF>& cgfs);

    /**
     * @brief Identify Cartesian shell components by their verified radial contractions.
     */
    void prepare_shells(const std::vector<CGF>& cgfs);

    /**
     * @brief Look up a prepared CGF pair by its packed index.
     */
    const Pair* pair(std::size_t index) const;

    /**
     * @brief Return the metadata for one contraction in the prepared basis.
     */
    const Contraction& contraction(std::size_t i) const { return contractions[i]; }

    /**
     * @brief Choose the kernel using the angular and reflection properties of two pairs.
     *
     * @param ij  Packed bra pair index
     * @param kl  Packed ket pair index
     *
     * @return Route and shell batching requirements for the quartet
     */
    QuartetRoute plan_quartet(std::size_t ij, std::size_t kl) const;

    /**
     * @brief Resolve the baseline engine for the requested ERI mode.
     */
    static EriEngine baseline_engine(EriEngine requested) {
        return requested == EriEngine::Auto ? EriEngine::Hellsing : requested;
    }

private:
    struct PairRouting {
        std::optional<Vec3> center;
        unsigned int angular_order = 0;
        unsigned int occupied_axes = 0;
        unsigned int odd_axes = 0;
        unsigned int common_planes = 0;
        bool all_d_or_higher = false;
        bool all_s_or_p = false;
    };

    /**
     * @brief Summarize the angular and reflection properties of a verified CGF pair.
     */
    static PairRouting describe_pair(const Contraction& a, const Contraction& b);

    std::vector<Contraction> contractions;
    std::vector<std::optional<Pair>> pairs;
    std::vector<PairRouting> pair_routing;
};

}  // namespace integrals::eri
