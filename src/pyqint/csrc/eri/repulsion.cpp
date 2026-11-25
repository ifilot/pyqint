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
#include "hellsing/hellsing.h"
#include "huzinaga/huzinaga.h"
#include "rys/rys.h"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <exception>
#include <memory>
#include <new>
#include <stdexcept>

using integrals::eri::EriBuildProfile;
using integrals::eri::EriOrchestrator;

namespace {

using Clock = std::chrono::steady_clock;
using integrals::eri::rys::ShellComponent;

using BasisPair = std::array<size_t, 2>;

struct QuartetJob {
    size_t result_index;
    size_t bra;
    size_t ket;
    EriEngine engine;
};

struct ShellPairComponent {
    size_t packed_index;
    BasisPair basis_indices;
};

struct ShellCandidate {
    QuartetJob fallback;
    ShellComponent component;
};

struct ShellBatch {
    size_t start;
    size_t count;
};

class QuartetPlanningMemoryError : public std::bad_alloc {
public:
    QuartetPlanningMemoryError(size_t requested, size_t available)
        : message("ERI quartet planning needs at least " + std::to_string(requested) +
                  " bytes of workspace, with " + std::to_string(available) + " bytes available.") {}

    const char* what() const noexcept override { return message.c_str(); }

private:
    std::string message;
};

struct QuartetPlanningMemory {
    size_t budget;
    size_t used = 0;
};

/** @brief Charge live planning buffers, including overlapping vector reallocations. */
template<class T>
struct QuartetPlanningAllocator {
    using value_type = T;
    QuartetPlanningMemory* memory;

    explicit QuartetPlanningAllocator(QuartetPlanningMemory* memory) noexcept : memory(memory) {}

    template<class U>
    QuartetPlanningAllocator(const QuartetPlanningAllocator<U>& other) noexcept : memory(other.memory) {}

    T* allocate(size_t count) {
        const size_t bytes = count * sizeof(T);
        if (bytes > memory->budget - memory->used) {
            throw QuartetPlanningMemoryError(memory->used + bytes, memory->budget);
        }
        T* buffer = std::allocator<T>{}.allocate(count);
        memory->used += bytes;
        return buffer;
    }

    void deallocate(T* buffer, size_t count) noexcept {
        std::allocator<T>{}.deallocate(buffer, count);
        memory->used -= count * sizeof(T);
    }

    template<class U>
    bool operator==(const QuartetPlanningAllocator<U>& other) const noexcept { return memory == other.memory; }

    template<class U>
    bool operator!=(const QuartetPlanningAllocator<U>& other) const noexcept { return !(*this == other); }
};

template<class T>
using PlanningVector = std::vector<T, QuartetPlanningAllocator<T>>;

constexpr size_t triangular(size_t n) {
    return n * (n + 1) / 2;
}

/**
 * @brief Group CGF pairs by canonical shell pair and orient each pair to match.
 */
std::vector<std::vector<ShellPairComponent>> group_shell_pairs(
    const std::vector<BasisPair>& basis_pairs, const EriOrchestrator& orchestrator) {
    std::vector<std::vector<ShellPairComponent>> groups(basis_pairs.size());
    for (size_t index = 0; index < basis_pairs.size(); ++index) {
        auto [i, j] = basis_pairs[index];
        size_t shell_i = orchestrator.contraction(i).shell_id;
        size_t shell_j = orchestrator.contraction(j).shell_id;
        if (shell_i > shell_j) {
            std::swap(shell_i, shell_j);
            std::swap(i, j);
        }
        groups[triangular(shell_j) + shell_i].push_back({index, {i, j}});
    }

    groups.erase(std::remove_if(groups.begin(), groups.end(),
                                [](const auto& group) { return group.empty(); }), groups.end());
    return groups;
}

/**
 * @brief Contract a primitive ERI kernel over four Gaussian functions.
 */
template<class Kernel>
double contract_quartet(const CGF& a, const CGF& b, const CGF& c, const CGF& d,
                         const std::array<std::size_t, 4>& counts, Kernel evaluate) {
    const auto& [count1, count2, count3, count4] = counts;
    double sum = 0.0;
    for (size_t i = 0; i < count1; ++i) {
        const GTO& g1 = a.get_gto(i);
        for (size_t j = 0; j < count2; ++j) {
            const GTO& g2 = b.get_gto(j);
            for (size_t k = 0; k < count3; ++k) {
                const GTO& g3 = c.get_gto(k);
                for (size_t l = 0; l < count4; ++l) {
                    const GTO& g4 = d.get_gto(l);
                    const double coeff = g1.get_coefficient() * g2.get_coefficient() *
                                         g3.get_coefficient() * g4.get_coefficient();
                    sum += g1.get_norm() * g2.get_norm() * g3.get_norm() * g4.get_norm() * coeff *
                           evaluate(g1, g2, g3, g4, i, j, k, l);
                }
            }
        }
    }
    return sum * a.get_contraction_norm() * b.get_contraction_norm() *
                 c.get_contraction_norm() * d.get_contraction_norm();
}

}  // namespace

/**
 * @brief Evaluate all two-electron integrals within the basis set.
 *
 * @param cgfs  Contracted Gaussian functions
 * @param planner_budget  Bytes available for ERI quartet planning
 *
 * @return dense two-electron tensor in row-major order
 */
std::vector<double> Integrator::evaluate_tei(const std::vector<CGF>& cgfs, size_t planner_budget) const {
    std::vector<double> tetensor(std::pow(cgfs.size(), 4));
    this->calculate_two_electron_integrals(cgfs, tetensor, nullptr, planner_budget);
    return tetensor;
}

std::vector<double> Integrator::evaluate_tei(const std::vector<CGF>& cgfs, EriBuildProfile& profile,
                                            size_t planner_budget) const {
    std::vector<double> tetensor(std::pow(cgfs.size(), 4));
    this->calculate_two_electron_integrals(cgfs, tetensor, &profile, planner_budget);
    return tetensor;
}

/**
 * @brief Compute unique ERIs in teindex order.
 *
 * @param cgfs      Contracted Gaussian functions
 * @param tedouble  Zero-initialized storage for the unique ERIs
 * @param profile   Optional stage timings and route counts
 * @param planner_budget  Bytes available for ERI quartet planning
 */
void Integrator::evaluate_packed_tei(const std::vector<CGF>& cgfs, double* tedouble,
                                     EriBuildProfile* profile, size_t planner_budget) const {
    const size_t sz = cgfs.size();
    const size_t pair_count = triangular(sz);
    if (profile) {
        *profile = {};
    }

    if (EriOrchestrator::baseline_engine(eri_engine) == EriEngine::Hellsing) {
        unsigned int max_l = 0;
        for (const auto& cgf : cgfs) {
            max_l = std::max(max_l, cgf.max_primitive_l());
        }
        hellsing_cache.ensure_lmax(max_l);
    }

    std::optional<EriOrchestrator> orchestrator;
    if (eri_engine == EriEngine::Auto) {
        const auto start = profile ? Clock::now() : Clock::time_point{};
        orchestrator.emplace(cgfs);
        orchestrator->prepare_pairs(cgfs);
        orchestrator->prepare_shells(cgfs);
        if (profile) {
            profile->pair_setup = std::chrono::duration<double>(Clock::now() - start).count();
        }
    }

    const auto planning_start = profile ? Clock::now() : Clock::time_point{};
    // Generate each canonical CGF pair in packed order.
    std::vector<BasisPair> basis_pairs;
    std::vector<size_t> pair_offsets;
    basis_pairs.reserve(pair_count);
    pair_offsets.reserve(pair_count);
    for (size_t j = 0; j < sz; ++j) {
        for (size_t i = 0; i <= j; ++i) {
            pair_offsets.push_back(triangular(basis_pairs.size()));
            basis_pairs.push_back({i, j});
        }
    }

    QuartetPlanningMemory memory{planner_budget};
    PlanningVector<QuartetJob> jobs{QuartetPlanningAllocator<QuartetJob>(&memory)};
    PlanningVector<ShellComponent> shell_components{QuartetPlanningAllocator<ShellComponent>(&memory)};
    PlanningVector<ShellBatch> shell_batches{QuartetPlanningAllocator<ShellBatch>(&memory)};
    size_t parity_zero_quartets = 0;
    size_t rys_quartets = 0;
    if (!orchestrator) {
        jobs.reserve(triangular(pair_count));
        for (size_t ket = 0; ket < pair_count; ++ket) {
            const size_t offset = pair_offsets[ket];
            for (size_t bra = 0; bra <= ket; ++bra) {
                jobs.push_back({offset + bra, bra, ket, eri_engine});
            }
        }
    } else {
        // Walk the shell-pair triangle. Within a diagonal block, visit its lower triangle.
        const auto groups = group_shell_pairs(basis_pairs, *orchestrator);
        PlanningVector<ShellCandidate> candidates{QuartetPlanningAllocator<ShellCandidate>(&memory)};
        for (size_t outer = 0; outer < groups.size(); ++outer) {
            const auto& ket_group = groups[outer];
            for (size_t inner = 0; inner <= outer; ++inner) {
                const auto& bra_group = groups[inner];
                candidates.clear();
                size_t minimum_components = 0;
                for (size_t k = 0; k < ket_group.size(); ++k) {
                    const auto& ket = ket_group[k];
                    const size_t end = inner == outer ? k + 1 : bra_group.size();
                    for (size_t i = 0; i < end; ++i) {
                        const auto& bra = bra_group[i];
                        const auto route = orchestrator->plan_quartet(bra.packed_index, ket.packed_index);
                        if (route.zero_by_parity) {
                            ++parity_zero_quartets;
                            continue;
                        }

                        const size_t first = std::min(bra.packed_index, ket.packed_index);
                        const size_t second = std::max(bra.packed_index, ket.packed_index);
                        const size_t index = pair_offsets[second] + first;
                        const QuartetJob job{index, first, second, route.engine};
                        if (route.min_shell_components == 0) {
                            rys_quartets += route.engine == EriEngine::Rys;
                            jobs.push_back(job);
                            continue;
                        }

                        minimum_components = route.min_shell_components;
                        const auto& [a, b] = bra.basis_indices;
                        const auto& [c, d] = ket.basis_indices;
                        candidates.push_back({job, {{&cgfs[a], &cgfs[b], &cgfs[c], &cgfs[d]}, index}});
                    }
                }

                if (candidates.empty()) {
                    continue;
                }
                if (candidates.size() >= minimum_components) {
                    rys_quartets += candidates.size();
                    shell_batches.push_back({shell_components.size(), candidates.size()});
                    for (const auto& candidate : candidates) {
                        shell_components.push_back(candidate.component);
                    }
                } else {
                    for (const auto& candidate : candidates) {
                        rys_quartets += candidate.fallback.engine == EriEngine::Rys;
                        jobs.push_back(candidate.fallback);
                    }
                }
            }
        }
    }

    // Keep exceptions on their OpenMP thread until the parallel region ends.
    const auto evaluation_start = profile ? Clock::now() : Clock::time_point{};
    if (profile) {
        profile->quartet_planning = std::chrono::duration<double>(evaluation_start - planning_start).count();
    }
    std::exception_ptr eri_error;
    size_t prepared_quartets = 0;
    const auto job_count = static_cast<std::ptrdiff_t>(jobs.size() + shell_batches.size());
    #pragma omp parallel for schedule(dynamic) reduction(+:prepared_quartets)
    for (std::ptrdiff_t s = 0; s < job_count; ++s) {
        try {
            const auto index = static_cast<size_t>(s);
            if (index >= jobs.size()) {
                const auto& batch = shell_batches[index - jobs.size()];
                integrals::eri::rys::evaluate_shell_quartet(
                    shell_components.data() + batch.start, batch.count, boys_function, tedouble);
                continue;
            }

            const auto& job = jobs[index];
            const auto& [i, j] = basis_pairs[job.bra];
            const auto& [k, l] = basis_pairs[job.ket];
            const bool use_pairs = orchestrator && job.engine == EriEngine::Hellsing;
            const auto* bra = use_pairs ? orchestrator->pair(job.bra) : nullptr;
            const auto* ket = use_pairs ? orchestrator->pair(job.ket) : nullptr;
            tedouble[job.result_index] = this->repulsion(cgfs[i], cgfs[j], cgfs[k], cgfs[l], job.engine, bra, ket);
            prepared_quartets += profile && bra && ket;
        } catch (...) {
            #pragma omp critical
            {
                if (!eri_error) {
                    eri_error = std::current_exception();
                }
            }
        }
    }
    if (eri_error) {
        std::rethrow_exception(eri_error);
    }

    if (profile) {
        profile->integral_evaluation = std::chrono::duration<double>(Clock::now() - evaluation_start).count();
        profile->parity_zero_quartets = parity_zero_quartets;
        profile->prepared_hellsing_quartets = prepared_quartets;
        profile->rys_quartets = eri_engine == EriEngine::Rys ? jobs.size() : rys_quartets;
        profile->rys_shell_batches = shell_batches.size();
    }
}

void Integrator::calculate_two_electron_integrals(const std::vector<CGF>& cgfs,
                                                  std::vector<double>& tetensor,
                                                  EriBuildProfile* profile, size_t planner_budget) const {
    const size_t sz = cgfs.size();
    std::vector<double> tedouble(triangular(triangular(sz)), 0.0);
    this->evaluate_packed_tei(cgfs, tedouble.data(), profile, planner_budget);
    const auto assembly_start = profile ? Clock::now() : Clock::time_point{};

    // Expand packed ERIs by rows using the dense-to-packed map of basis pairs.
    const size_t dense_pair_count = sz * sz;
    std::vector<size_t> pair_ids(dense_pair_count);
    std::vector<size_t> pair_offsets(triangular(sz));
    for (size_t j = 0; j < sz; ++j) {
        for (size_t i = 0; i <= j; ++i) {
            const size_t id = triangular(j) + i;
            pair_ids[i * sz + j] = pair_ids[j * sz + i] = id;
            pair_offsets[id] = triangular(id);
        }
    }

    const auto row_count = static_cast<std::ptrdiff_t>(dense_pair_count);
    #pragma omp parallel for schedule(static) if(dense_pair_count >= 1024)
    for (std::ptrdiff_t row = 0; row < row_count; ++row) {
        double* output = tetensor.data() + static_cast<size_t>(row) * dense_pair_count;
        const size_t bra = pair_ids[row];
        const size_t offset = pair_offsets[bra];
        for (size_t col = 0; col < dense_pair_count; ++col) {
            const size_t ket = pair_ids[col];
            output[col] = tedouble[ket <= bra ? offset + ket : pair_offsets[ket] + bra];
        }
    }
    if (profile) {
        profile->tensor_assembly = std::chrono::duration<double>(Clock::now() - assembly_start).count();
    }
}

/**************************************************************************
 * TWO-ELECTRON INTEGRALS
 **************************************************************************/

/**
 * @brief Calculate two-electron repulsion integral of four CGFs.
 *
 * @param cgf1  Contracted Gaussian Function
 * @param cgf2  Contracted Gaussian Function
 * @param cgf3  Contracted Gaussian Function
 * @param cgf4  Contracted Gaussian Function
 *
 * Calculates the value of < cgf1 | cgf2 | cgf3 | cgf4 >
 *
 * @return value of the repulsion integral
 */
double Integrator::repulsion(const CGF& cgf1,const CGF& cgf2,const CGF& cgf3,const CGF& cgf4) const {
    return this->repulsion(cgf1, cgf2, cgf3, cgf4, EriOrchestrator::baseline_engine(eri_engine));
}

double Integrator::repulsion(const CGF& cgf1, const CGF& cgf2, const CGF& cgf3, const CGF& cgf4,
                             EriEngine engine,
                             const EriOrchestrator::Pair* bra, const EriOrchestrator::Pair* ket) const {
    const std::array<size_t, 4> counts{cgf1.size(), cgf2.size(), cgf3.size(), cgf4.size()};
    if (engine == EriEngine::Hellsing && bra && ket) {
        return contract_quartet(cgf1, cgf2, cgf3, cgf4, counts,
            [&](const GTO& a, const GTO& b, const GTO& c, const GTO& d,
                size_t i, size_t j, size_t k, size_t l) {
                return integrals::eri::hellsing::evaluate_prepared(
                    a, b, c, d, bra->primitives[i * counts[1] + j],
                    ket->primitives[k * counts[3] + l],
                    bra->displacement, ket->displacement, boys_function, hellsing_cache);
            });
    }
    return contract_quartet(cgf1, cgf2, cgf3, cgf4, counts,
        [&](const GTO& a, const GTO& b, const GTO& c, const GTO& d, size_t, size_t, size_t, size_t) {
            return this->repulsion(a, b, c, d, engine);
        });
}

/**
 * @brief Calculate derivative of the two-electron repulsion integral for four CGFs.
 *
 * @param cgf1     Contracted Gaussian Function
 * @param cgf2     Contracted Gaussian Function
 * @param cgf3     Contracted Gaussian Function
 * @param cgf4     Contracted Gaussian Function
 * @param nucleus  Nucleus coordinates
 * @param coord    Cartesian direction (0=x,1=y,2=z)
 *
 * Calculates the value of d/dcx < cgf1 | cgf2 | cgf3 | cgf4 >
 *
 * @return value of the repulsion integral derivative
 */
double Integrator::repulsion_deriv(const CGF& cgf1, const CGF& cgf2, const CGF& cgf3, const CGF& cgf4,
    const Vec3& nucleus, unsigned int coord) const {
    return this->repulsion_deriv(prepare_cgf_derivative(cgf1, nucleus), prepare_cgf_derivative(cgf2, nucleus),
                                 prepare_cgf_derivative(cgf3, nucleus), prepare_cgf_derivative(cgf4, nucleus), coord);
}

double Integrator::repulsion_deriv(const CgfDerivative& a, const CgfDerivative& b,
                                    const CgfDerivative& c, const CgfDerivative& d, unsigned int coord) const {
    if ((a.none_move() && b.none_move() && c.none_move() && d.none_move()) ||
        (a.all_move() && b.all_move() && c.all_move() && d.all_move())) {
        return 0.0;
    }

    double sum = 0.0;
    for (unsigned int i = 0; i < a.cgf.size(); ++i) {
        const GTO& g1 = a.cgf.get_gto(i);
        const bool a_moves = a.moves(i);
        const double weight1 = g1.get_norm() * g1.get_coefficient();
        for (unsigned int j = 0; j < b.cgf.size(); ++j) {
            const GTO& g2 = b.cgf.get_gto(j);
            const bool b_moves = b.moves(j);
            const double weight2 = weight1 * g2.get_norm() * g2.get_coefficient();
            for (unsigned int k = 0; k < c.cgf.size(); ++k) {
                const GTO& g3 = c.cgf.get_gto(k);
                const bool c_moves = c.moves(k);
                const double weight3 = weight2 * g3.get_norm() * g3.get_coefficient();
                for (unsigned int l = 0; l < d.cgf.size(); ++l) {
                    const bool d_moves = d.moves(l);
                    // A primitive quartet is invariant under a common translation.
                    if (a_moves == b_moves && b_moves == c_moves && c_moves == d_moves) {
                        continue;
                    }

                    const GTO& g4 = d.cgf.get_gto(l);
                    const double weight4 = weight3 * g4.get_norm() * g4.get_coefficient();
                    const double t1 = a_moves ? this->repulsion_deriv(g1, g2, g3, g4, coord) : 0.0;
                    const double t2 = b_moves ? this->repulsion_deriv(g2, g1, g3, g4, coord) : 0.0;
                    const double t3 = c_moves ? this->repulsion_deriv(g3, g4, g1, g2, coord) : 0.0;
                    const double t4 = d_moves ? this->repulsion_deriv(g4, g3, g1, g2, coord) : 0.0;
                    sum += weight4 * (t1 + t2 + t3 + t4);
                }
            }
        }
    }

    double result = sum * a.cgf.get_contraction_norm() * b.cgf.get_contraction_norm() *
                          c.cgf.get_contraction_norm() * d.cgf.get_contraction_norm();
    const double log_norm_deriv = a.log_norm_deriv[coord] + b.log_norm_deriv[coord] +
                                  c.log_norm_deriv[coord] + d.log_norm_deriv[coord];
    if (log_norm_deriv != 0.0) {
        result += log_norm_deriv * this->repulsion(a.cgf, b.cgf, c.cgf, d.cgf);
    }
    return result;
}

/**
 * @brief Calculate two-electron repulsion integral of four GTOs.
 *
 * @param gto1  Gaussian Type Orbital
 * @param gto2  Gaussian Type Orbital
 * @param gto3  Gaussian Type Orbital
 * @param gto4  Gaussian Type Orbital
 *
 * Calculates the value of < gto1 | gto2 | gto3 | gto4 >
 *
 * @return value of the repulsion integral
 */
double Integrator::repulsion(const GTO &gto1, const GTO &gto2, const GTO &gto3, const GTO &gto4) const {
    return this->repulsion(gto1, gto2, gto3, gto4, EriOrchestrator::baseline_engine(eri_engine));
}

double Integrator::repulsion(const GTO& gto1, const GTO& gto2, const GTO& gto3, const GTO& gto4,
                             EriEngine engine) const {
    using namespace integrals::eri;

    switch (engine) {
        case EriEngine::Hellsing:
            return hellsing::evaluate_cached(gto1, gto2, gto3, gto4, boys_function, hellsing_cache);
        case EriEngine::HellsingDirect:
            return hellsing::evaluate_direct(gto1, gto2, gto3, gto4, boys_function);
        case EriEngine::Huzinaga:
            return huzinaga::evaluate_block(gto1, gto2, gto3, gto4, boys_function);
        case EriEngine::HuzinagaScalar:
            return huzinaga::evaluate_scalar(gto1, gto2, gto3, gto4, boys_function);
        case EriEngine::Rys:
            return rys::evaluate(gto1, gto2, gto3, gto4, boys_function);
        case EriEngine::Auto:
            break;
    }
    throw std::logic_error("Invalid electron repulsion engine");
}

/**
 * @brief Calculate repulsion derivative integral of four GTOs.
 *
 * @param gto1   Gaussian Type Orbital
 * @param gto2   Gaussian Type Orbital
 * @param gto3   Gaussian Type Orbital
 * @param gto4   Gaussian Type Orbital
 * @param coord  Cartesian direction (0=x,1=y,2=z)
 *
 * Calculates the value of < d/dx gto1 | gto2 | gto3 | gto4 >
 *
 * @return value of the overlap derivative integral
 */
double Integrator::repulsion_deriv(const GTO& gto1, const GTO& gto2, const GTO &gto3, const GTO &gto4, unsigned int coord) const {
    std::array<unsigned int, 3> gto_ang = {gto1.get_l(), gto1.get_m(), gto1.get_n()};

    if(gto_ang[coord] != 0) {

        gto_ang[coord] += 1; // calculate l+1 term
        GTO gto1_new1(gto1.get_coefficient(), gto1.get_position()[0], gto1.get_position()[1], gto1.get_position()[2],
                      gto1.get_alpha(), gto_ang[0], gto_ang[1], gto_ang[2]);
        double term_plus = this->repulsion(gto1_new1, gto2, gto3, gto4);

        gto_ang[coord] -= 2; // calculate l-1 term
        GTO gto1_new2(gto1.get_coefficient(), gto1.get_position()[0], gto1.get_position()[1], gto1.get_position()[2],
                      gto1.get_alpha(), gto_ang[0], gto_ang[1], gto_ang[2]);
        double term_min = this->repulsion(gto1_new2, gto2, gto3, gto4);

        gto_ang[coord] += 1; // recover l

        return 2.0 * gto1.get_alpha() * term_plus - gto_ang[coord] * term_min;
    } else { // s-type GTO
        gto_ang[coord] += 1;
        GTO gto1_new(gto1.get_coefficient(), gto1.get_position()[0], gto1.get_position()[1], gto1.get_position()[2],
                     gto1.get_alpha(), gto_ang[0], gto_ang[1], gto_ang[2]);
        double term = this->repulsion(gto1_new, gto2, gto3, gto4);
        return 2.0 * gto1.get_alpha() * term;
    }
}

/**
 * @brief Compute the flattened two-electron integral index.
 *
 * @param i  basis function index
 * @param j  basis function index
 * @param k  basis function index
 * @param l  basis function index
 *
 * @return flattened index in the packed TEI array
 */
size_t Integrator::teindex(size_t i, size_t j, size_t k, size_t l) const {
    if(i < j) {
        std::swap(i,j);
    }
    if(k < l) {
        std::swap(k,l);
    }

    size_t ij = i * (i + 1) / 2 + j;
    size_t kl = k * (k + 1) / 2 + l;

    if(ij < kl) {
        std::swap(ij,kl);
    }

    return ij * (ij + 1) / 2 + kl;
}

/**
 * @brief Ensure the Hellsing cache supports the basis set.
 *
 * @param cgf1  Contracted Gaussian Function
 * @param cgf2  Contracted Gaussian Function
 * @param cgf3  Contracted Gaussian Function
 * @param cgf4  Contracted Gaussian Function
 */
void Integrator::ensure_hellsing_cache(const CGF& cgf1, const CGF& cgf2, const CGF& cgf3, const CGF& cgf4) {
    if (EriOrchestrator::baseline_engine(eri_engine) != EriEngine::Hellsing) {
        return;
    }
    unsigned int lmax = std::max({cgf1.max_primitive_l(),
                                  cgf2.max_primitive_l(),
                                  cgf3.max_primitive_l(),
                                  cgf4.max_primitive_l()});

    // one more to also include derivatives
    this->hellsing_cache.ensure_lmax(lmax+1);
}
