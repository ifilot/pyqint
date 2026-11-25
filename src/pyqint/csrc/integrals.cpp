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
#include "gaussian_integrals.h"

#include <exception>

using integrals::eri::EriOrchestrator;

/**
 * @brief Construct an Integrator instance and capture build-time information.
 *
 * The compile information (compiler, version, OpenMP, date, time) is
 * captured once at construction time using compile-time macros.
 *
 * @param lmax        Maximum angular momentum to cache initially
 * @param nu_max      Maximum Boys function order to tabulate
 * @param eri_engine  Primitive electron repulsion engine
 */
Integrator::Integrator(int lmax, int nu_max, EriEngine eri_engine)
    : boys_function(nu_max),
      eri_engine(eri_engine),
      hellsing_cache(EriOrchestrator::baseline_engine(eri_engine) == EriEngine::Hellsing ? lmax : 0)
{
    // ---------------------------------------------------------------------
    // Compile date & time
    // ---------------------------------------------------------------------
#ifdef __DATE__
    compile_date = __DATE__;
#else
    compile_date = "unknown";
#endif

#ifdef __TIME__
    compile_time = __TIME__;
#else
    compile_time = "unknown";
#endif

    // ---------------------------------------------------------------------
    // Compiler type & version
    // ---------------------------------------------------------------------
#if defined(_MSC_VER)
    compiler_type    = "MSVC";
    compiler_version = std::to_string(_MSC_VER);

#elif defined(__clang__)
    compiler_type    = "Clang";
    compiler_version = __clang_version__;

#elif defined(__GNUC__)
    compiler_type    = "GCC";
    compiler_version = __VERSION__;

#else
    compiler_type    = "unknown";
    compiler_version = "unknown";
#endif

    // ---------------------------------------------------------------------
    // OpenMP version
    // ---------------------------------------------------------------------
#ifdef _OPENMP
    openmp_version = std::to_string(_OPENMP);
#else
    openmp_version = "disabled";
#endif
}

/**
 * @brief Evaluate all one- and two-electron integrals for a basis set.
 *
 * @param cgfs     Contracted Gaussian functions
 * @param charges  Nuclear charges
 * @param px       Nuclear x coordinates
 * @param py       Nuclear y coordinates
 * @param pz       Nuclear z coordinates
 * @param planner_budget  Bytes available for ERI quartet planning
 *
 * @return packed vector of overlap, kinetic, nuclear, and TEI values
 */
std::vector<double> Integrator::evaluate_cgfs(const std::vector<CGF>& cgfs,
                                              const std::vector<int>& charges,
                                              const std::vector<double>& px,
                                              const std::vector<double>& py,
                                              const std::vector<double>& pz,
                                              size_t planner_budget) const {
    auto results = this->evaluate_one_electron_integrals(cgfs, charges, px, py, pz);
    const auto tetensor = this->evaluate_tei(cgfs, planner_budget);
    results.insert(results.end(), tetensor.begin(), tetensor.end());
    return results;
}

std::vector<double> Integrator::evaluate_one_electron_integrals(const std::vector<CGF>& cgfs,
                                                                const std::vector<int>& charges,
                                                                const std::vector<double>& px,
                                                                const std::vector<double>& py,
                                                                const std::vector<double>& pz) const {
    const size_t sz = cgfs.size();

    std::vector<double> S(sz*sz, 0.0);
    std::vector<double> T(sz*sz, 0.0);
    std::vector<double> V(sz*sz, 0.0);
    std::vector<double> Vn(charges.size() * sz * sz, 0.0);

    // calculate the integral values using the integrator class
    #pragma omp parallel for schedule(dynamic)
    for(int i=0; i<(int)sz; i++) {  // have to use signed int for MSVC OpenMP here
        for(unsigned int j=i; j<sz; j++) {
            S[i*sz+j] = S[j*sz+i] = this->overlap(cgfs[i], cgfs[j]);
            T[i*sz+j] = T[j*sz+i] = this->kinetic(cgfs[i], cgfs[j]);
            for(unsigned int k=0; k<charges.size(); k++) {
                // this is to avoid the race condition
                Vn[k*sz*sz + i*sz+j] = Vn[k*sz*sz + j*sz+i] = this->nuclear(cgfs[i], cgfs[j], Vec3(px[k], py[k], pz[k]), charges[k]);
            }
        }
    }

    // combine all nuclear attraction integrals
    for(unsigned int i=0; i<sz; i++) {
        for(unsigned int j=0; j<sz; j++) {
            for(unsigned int k=0; k<charges.size(); k++) {
                V[i*sz+j] += Vn[k*sz*sz + i*sz+j];
            }
        }
    }

    // Matrices occupy consecutive slices in the Python interface.
    std::vector<double> results;
    results.reserve(3 * sz * sz);
    results.insert(results.end(), S.begin(), S.end());
    results.insert(results.end(), T.begin(), T.end());
    results.insert(results.end(), V.begin(), V.end());
    return results;
}

Integrator::CgfDerivative Integrator::prepare_cgf_derivative(const CGF& cgf, const Vec3& nucleus,
                                                            const std::optional<Vec3>& shared_center) const {
    if (shared_center) {
        return {cgf, same_center(*shared_center, nucleus), {}, {}};
    }

    const bool first_moves = cgf.size() != 0 && same_center(cgf.get_gto(0).get_position(), nucleus);
    CgfDerivative result{cgf, first_moves, {}, {}};
    unsigned int i = 1;
    while (i < cgf.size() && same_center(cgf.get_gto(i).get_position(), nucleus) == first_moves) {
        ++i;
    }
    if (i >= cgf.size()) {
        return result;
    }

    result.primitive_moves.resize(cgf.size(), first_moves);
    for (; i < cgf.size(); ++i) {
        result.primitive_moves[i] = same_center(cgf.get_gto(i).get_position(), nucleus);
    }

    // N = <phi|phi>^(-1/2). Only moving/fixed pairs contribute to d(log N).
    // Counting each such pair once cancels the 1/2 from differentiation.
    const double norm = cgf.get_contraction_norm();
    for (unsigned int i = 0; i < cgf.size(); ++i) {
        if (!result.moves(i)) {
            continue;
        }
        const GTO& moving = cgf.get_gto(i);
        for (unsigned int j = 0; j < cgf.size(); ++j) {
            if (result.moves(j)) {
                continue;
            }
            const GTO& fixed = cgf.get_gto(j);
            const double weight = (norm * moving.get_norm() * moving.get_coefficient()) *
                                  (norm * fixed.get_norm() * fixed.get_coefficient());
            for (unsigned int coord = 0; coord < 3; ++coord) {
                result.log_norm_deriv[coord] -= weight * this->overlap_deriv_gto(moving, fixed, coord);
            }
        }
    }
    return result;
}

double Integrator::two_center_deriv(const CgfDerivative& bra, const CgfDerivative& ket, unsigned int coord,
                                    ContractedIntegral integral, PrimitiveDerivative derivative) const {
    if ((bra.none_move() && ket.none_move()) || (bra.all_move() && ket.all_move())) {
        return 0.0;
    }

    double sum = 0.0;
    for (unsigned int i = 0; i < bra.cgf.size(); ++i) {
        const bool bra_moves = bra.moves(i);
        const GTO& g1 = bra.cgf.get_gto(i);
        const double weight1 = g1.get_norm() * g1.get_coefficient();
        for (unsigned int j = 0; j < ket.cgf.size(); ++j) {
            if (bra_moves == ket.moves(j)) {
                continue;
            }
            const GTO& g2 = ket.cgf.get_gto(j);
            const double deriv = bra_moves ? (this->*derivative)(g1, g2, coord) :
                                             (this->*derivative)(g2, g1, coord);
            sum += weight1 * g2.get_norm() * g2.get_coefficient() * deriv;
        }
    }

    double result = sum * bra.cgf.get_contraction_norm() * ket.cgf.get_contraction_norm();
    const double log_norm_deriv = bra.log_norm_deriv[coord] + ket.log_norm_deriv[coord];
    if (log_norm_deriv != 0.0) {
        result += log_norm_deriv * (this->*integral)(bra.cgf, ket.cgf);
    }
    return result;
}

/**
 * @brief Evaluate geometric derivatives for all cgfs in buffer.
 *
 * @param cgfs     Contracted Gaussian functions
 * @param charges  Nuclear charges
 * @param px       Nuclear x coordinates
 * @param py       Nuclear y coordinates
 * @param pz       Nuclear z coordinates
 *
 * @return packed vector of overlap, kinetic, nuclear, and TEI derivatives
 */
std::vector<double> Integrator::evaluate_geometric_derivatives(const std::vector<CGF>& cgfs,
                                                               const std::vector<int>& charges,
                                                               const std::vector<double>& px,
                                                               const std::vector<double>& py,
                                                               const std::vector<double>& pz) const {
    if (cgfs.empty() || charges.empty()) {
        return {};
    }
    if (EriOrchestrator::baseline_engine(eri_engine) == EriEngine::Hellsing) {
        unsigned int max_l = 0;
        for (const auto& cgf : cgfs) {
            max_l = std::max(max_l, cgf.max_primitive_l());
        }
        // Derivatives can raise one Cartesian power by one.
        hellsing_cache.ensure_lmax(max_l + 1);
    }

    const size_t sz = cgfs.size();
    std::optional<EriOrchestrator> orchestrator;
    if (eri_engine == EriEngine::Auto) {
        orchestrator.emplace(cgfs);
    }
    const std::optional<Vec3> no_shared_center;

    std::vector<Vec3> nuclei;
    nuclei.reserve(charges.size());
    std::vector<CgfDerivative> derivatives;
    derivatives.reserve(charges.size() * sz);
    for (size_t n = 0; n < charges.size(); ++n) {
        nuclei.emplace_back(px[n], py[n], pz[n]);
        for (size_t i = 0; i < sz; ++i) {
            const auto& center = orchestrator ? orchestrator->contraction(i).shared_center : no_shared_center;
            derivatives.push_back(prepare_cgf_derivative(cgfs[i], nuclei.back(), center));
        }
    }

    // One matrix per integral, nucleus, and Cartesian direction.
    const size_t datasz = charges.size() * 3 * sz * sz;
    std::vector<double> S(datasz, 0.0);
    std::vector<double> T(datasz, 0.0);
    std::vector<double> V(datasz, 0.0);

    // calculate the integral values using the integrator class
    for(unsigned int n=0; n<charges.size(); n++) { // loop over nuclei
        for(unsigned int k=0; k<3; k++) { // loop over directions

            // build container for summation over nuclei
            std::vector<double> Vn(charges.size() * sz * sz, 0.0);

            #pragma omp parallel for schedule(dynamic)
            for(int i=0; i<(int)sz; i++) {  // have to use signed int for MSVC OpenMP here
                for(int j=0; j<(int)sz; j++) {
                    const size_t idx = (n*3+k) * sz * sz + i * sz + j;
                    const CgfDerivative& bra = derivatives[n * sz + i];
                    const CgfDerivative& ket = derivatives[n * sz + j];
                    S[idx] = this->two_center_deriv(bra, ket, k, &Integrator::overlap, &Integrator::overlap_deriv_gto);
                    T[idx] = this->two_center_deriv(bra, ket, k, &Integrator::kinetic, &Integrator::kinetic_deriv_gto);

                    for(unsigned int l=0; l<charges.size(); l++) {
                        Vn[l*sz*sz + i*sz + j] = this->nuclear_deriv(
                            bra, ket, nuclei[l], charges[l], same_center(nuclei[l], nuclei[n]), k);
                    }
                }
            }

            // combine all nuclear attraction integrals
            for(unsigned int i=0; i<sz; i++) {
                for(unsigned int j=0; j<sz; j++) {
                    for(unsigned int l=0; l<charges.size(); l++) {
                        V[(n*3+k) * sz * sz + i * sz + j] += Vn[l*sz*sz + i*sz + j];
                    }
                }
            }
        }
    }

    // // calculate all two-electron integrals
    const size_t max_teints = this->teindex(sz-1,sz-1,sz-1,sz-1) + 1;
    std::vector<double> tedouble(max_teints * charges.size() * 3, -1.0);

    // it is more efficient to first 'unroll' the fourfold nested loop
    // into a single vector of jobs to execute
    std::vector<std::array<size_t, 7>> jobs;
    jobs.reserve(tedouble.size());
    for(size_t n=0; n<charges.size(); n++) {        // nuclear derivatives
        for(size_t d=0; d<3; d++) {                 // directions
            for(size_t i=0; i<sz; i++) {
                for(size_t j=0; j<sz; j++) {
                    size_t ij = i*(i+1)/2 + j;
                    for(size_t k=0; k<sz; k++) {
                        for(size_t l=0; l<sz; l++) {
                            size_t kl = k * (k+1)/2 + l;
                            if(ij <= kl) {
                                size_t idx = this->teindex(i,j,k,l) + (n*3 + d) * max_teints;

                                if(idx >= tedouble.size()) {
                                    throw std::runtime_error("Process tried to access illegal array position");
                                }

                                if(tedouble[idx] < 0.0) {
                                    tedouble[idx] = 1.0;
                                    jobs.emplace_back(std::array<size_t, 7>({idx, i, j, k, l, n, d}));
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Keep exceptions on their OpenMP thread until the parallel region ends.
    std::exception_ptr eri_error;
    #pragma omp parallel for schedule(dynamic)
    for (std::ptrdiff_t s = 0; s < static_cast<std::ptrdiff_t>(jobs.size()); ++s) {
        const size_t idx = jobs[s][0];
        const size_t i = jobs[s][1];
        const size_t j = jobs[s][2];
        const size_t k = jobs[s][3];
        const size_t l = jobs[s][4];
        const size_t n = jobs[s][5];
        const size_t d = jobs[s][6];
        try {
            tedouble[idx] = this->repulsion_deriv(
                derivatives[n * sz + i], derivatives[n * sz + j],
                derivatives[n * sz + k], derivatives[n * sz + l], d);
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

    // package everything into results vector, will be unpacked in
    // connected Python class
    std::vector<double> results(charges.size() * 3 * sz * sz * 3 + tedouble.size());
    for(size_t n=0; n<charges.size(); n++) {        // loop over nuclei
        for(size_t d=0; d<3; d++) {                 // loop over directions
            size_t spos = (n * 3 + d) * sz * sz;
            std::memcpy(&results[spos], &S[(n*3+d) * sz * sz], sz * sz * sizeof(double));

            size_t tpos = spos + charges.size() * 3 * sz * sz;
            std::memcpy(&results[tpos], &T[(n*3+d) * sz * sz], sz * sz * sizeof(double));

            size_t vpos = tpos + charges.size() * 3 * sz * sz;
            std::memcpy(&results[vpos], &V[(n*3+d) * sz * sz], sz * sz * sizeof(double));
        }
    }
    std::memcpy(&results[charges.size() * 3 * sz * sz * 3], tedouble.data(), tedouble.size() * sizeof(double));

    return results;
}

/**************************************************************************
 * OVERLAP INTEGRALS
 **************************************************************************/

/**
 * @brief Calculate overlap integral of two CGFs.
 *
 * @param cgf1  Contracted Gaussian Function 1
 * @param cgf2  Contracted Gaussian Function 2
 *
 * Calculates the value of < cgf1 | cgf2 >
 *
 * @return value of the overlap integral
 */
double Integrator::overlap(const CGF& cgf1, const CGF& cgf2) const {
    double sum = 0.0;

    // loop over all GTOs inside the CGF, calculate the overlap integrals
    // and sum all the integral values
    for (unsigned int k = 0; k < cgf1.size(); k++) {
        for (unsigned int l = 0; l < cgf2.size(); l++) {
            sum += cgf1.get_norm_gto(k) *
                   cgf2.get_norm_gto(l) *
                   cgf1.get_coefficient_gto(k) *
                   cgf2.get_coefficient_gto(l) *
                   this->overlap_gto(cgf1.get_gto(k), cgf2.get_gto(l));
        }
    }
    const auto norm1 = cgf1.get_contraction_norm();
    const auto norm2 = cgf2.get_contraction_norm();

    return sum * norm1 * norm2;
}

/**
 * @brief Calculate overlap integral of two GTOs.
 *
 * @param gto1  Gaussian Type Orbital 1
 * @param gto2  Gaussian Type Orbital 2
 *
 * Calculates the value of < gto1 | gto2 >
 *
 * @return value of the overlap integral
 */
double Integrator::overlap_gto(const GTO& gto1, const GTO& gto2) const {
    return this->overlap(gto1.get_alpha(), gto1.get_l(), gto1.get_m(), gto1.get_n(), gto1.get_position(),
                         gto2.get_alpha(), gto2.get_l(), gto2.get_m(), gto2.get_n(), gto2.get_position());
}

/**
 * @brief Calculate the geometric derivative of the overlap integral.
 *
 * @param cgf1     Contracted Gaussian Function 1
 * @param cgf2     Contracted Gaussian Function 2
 * @param nucleus  nucleus position
 * @param coord    Cartesian direction (0=x,1=y,2=z)
 *
 * Calculates the value of d/dx < cgf1 | cgf2 >
 *
 * @return value of the overlap derivative
 */
double Integrator::overlap_deriv(const CGF& cgf1,
                                 const CGF& cgf2,
                                 const Vec3& nucleus,
                                 unsigned int coord) const {
    return this->two_center_deriv(prepare_cgf_derivative(cgf1, nucleus), prepare_cgf_derivative(cgf2, nucleus),
                                  coord, &Integrator::overlap, &Integrator::overlap_deriv_gto);
}

/**
 * @brief Calculate the geometric derivative of the overlap integral for two GTOs.
 *
 * @param gto1   Gaussian Type Orbital
 * @param gto2   Gaussian Type Orbital
 * @param coord  Cartesian direction (0=x,1=y,2=z)
 *
 * @return value of the overlap derivative
 */
double Integrator::overlap_deriv_gto(const GTO& gto1,
                                     const GTO& gto2,
                                     unsigned int coord) const {
    std::array<unsigned int, 3> gto_ang = {gto1.get_l(),
                                           gto1.get_m(),
                                           gto1.get_n()};

    if(gto_ang[coord] != 0) {
        gto_ang[coord] += 1; // calculate l+1 term
        double term_plus = this->overlap(gto1.get_alpha(),
                                         gto_ang[0], gto_ang[1], gto_ang[2],
                                         gto1.get_position(),
                                         gto2.get_alpha(),
                                         gto2.get_l(), gto2.get_m(), gto2.get_n(),
                                         gto2.get_position());

        gto_ang[coord] -= 2; // calculate l-1 term

        double term_min = this->overlap(gto1.get_alpha(),
                                        gto_ang[0], gto_ang[1], gto_ang[2],
                                        gto1.get_position(),
                                        gto2.get_alpha(),
                                        gto2.get_l(), gto2.get_m(), gto2.get_n(),
                                        gto2.get_position());
        gto_ang[coord] += 1; // recover l

        return 2.0 * gto1.get_alpha() * term_plus - gto_ang[coord] * term_min;
    } else { // s-type GTO
        gto_ang[coord] += 1;
        double term1 = this->overlap(gto1.get_alpha(),
                                     gto_ang[0], gto_ang[1], gto_ang[2],
                                     gto1.get_position(),
                                     gto2.get_alpha(),
                                     gto2.get_l(), gto2.get_m(), gto2.get_n(),
                                     gto2.get_position());
        return 2.0 * gto1.get_alpha() * term1;
    }
}

/**************************************************************************
 * DIPOLE INTEGRALS
 **************************************************************************/

std::vector<double> Integrator::evaluate_dipoles(const std::vector<CGF>& cgfs,
                                                 double x, double y, double z) const {
    const size_t n = cgfs.size();
    const std::array<double, 3> origin{x, y, z};
    std::vector<double> values(3 * n * n);
    #pragma omp parallel for schedule(dynamic) if(n >= 16)
    for (std::ptrdiff_t row = 0; row < static_cast<std::ptrdiff_t>(n); ++row) {
        const auto i = static_cast<size_t>(row);
        for (size_t j = 0; j <= i; ++j) {
            for (unsigned int axis = 0; axis < origin.size(); ++axis) {
                const size_t offset = axis * n * n;
                const double value = this->dipole(cgfs[i], cgfs[j], axis, origin[axis]);
                values[offset + i * n + j] = values[offset + j * n + i] = value;
            }
        }
    }
    return values;
}

/**
 * @brief Calculate dipole integral of two CGFs.
 *
 * @param cgf1  Contracted Gaussian Function
 * @param cgf2  Contracted Gaussian Function
 * @param cc    Cartesian direction (0=x,1=y,2=z)
 * @param ref   Reference position
 *
 * Calculates the value of < cgf1 | (cc - ref) | cgf2 >
 *
 * @return value of the dipole integral
 */
double Integrator::dipole(const CGF& cgf1,
                          const CGF& cgf2,
                          unsigned int cc,
                          double ref) const {
    double sum = 0.0;

    // loop over all GTOs inside the CGF, calculate the dipole integrals
    // and sum all the integral values
    for(unsigned int k = 0; k < cgf1.size(); k++) {
        for(unsigned int l = 0; l < cgf2.size(); l++) {
            sum += cgf1.get_norm_gto(k) *
                   cgf2.get_norm_gto(l) *
                   cgf1.get_coefficient_gto(k) *
                   cgf2.get_coefficient_gto(l) *
                   this->dipole_gto(cgf1.get_gto(k), cgf2.get_gto(l), cc, ref);
        }
    }

    const auto norm1 = cgf1.get_contraction_norm();
    const auto norm2 = cgf2.get_contraction_norm();

    return sum * norm1 * norm2;
}

/**
 * @brief Calculate dipole integral of two GTOs.
 *
 * @param gto1  Gaussian Type Orbital
 * @param gto2  Gaussian Type Orbital
 * @param cc    Cartesian direction (0=x,1=y,2=z)
 * @param ref   Reference position
 *
 * Calculates the value of < gto1 | (cc - ref) | gto2 >
 *
 * @return value of the dipole integral
 */
double Integrator::dipole_gto(const GTO& gto1,
                              const GTO& gto2,
                              unsigned int cc,
                              double ref) const {
    return this->dipole(gto1.get_alpha(), gto1.get_l(), gto1.get_m(), gto1.get_n(), gto1.get_position(),
                        gto2.get_alpha(), gto2.get_l(), gto2.get_m(), gto2.get_n(), gto2.get_position(),
                        cc, ref);
}

/**************************************************************************
 * KINETIC INTEGRALS
 **************************************************************************/

/**
 * @brief Calculate kinetic integral of two CGFs.
 *
 * @param cgf1  Contracted Gaussian Function 1
 * @param cgf2  Contracted Gaussian Function 2
 *
 * Calculates the value of < cgf1 | T | cgf2 >
 *
 * @return value of the kinetic integral
 */
double Integrator::kinetic(const CGF& cgf1, const CGF& cgf2) const {
    double sum = 0.0;

    // loop over all GTOs inside the CGF, calculate the kinetic integrals
    // and sum all the integral values
    for(unsigned int k = 0; k < cgf1.size(); k++) {
        for(unsigned int l = 0; l < cgf2.size(); l++) {
            sum += cgf1.get_norm_gto(k) *
                   cgf2.get_norm_gto(l) *
                   cgf1.get_coefficient_gto(k) *
                   cgf2.get_coefficient_gto(l) *
                   this->kinetic_gto(cgf1.get_gto(k), cgf2.get_gto(l) );
        }
    }

    const auto norm1 = cgf1.get_contraction_norm();
    const auto norm2 = cgf2.get_contraction_norm();

    return sum * norm1 * norm2;
}

/**
 * @brief Calculate kinetic integral of two GTOs.
 *
 * @param gto1  Gaussian Type Orbital 1
 * @param gto2  Gaussian Type Orbital 2
 *
 * Calculates the value of < gto1 | T | gto2 >
 *
 * @return value of the kinetic integral
 */
double Integrator::kinetic_gto(const GTO& gto1, const GTO& gto2) const {
    double term0 = gto2.get_alpha() *
                   (2.0 * ( gto2.get_l() + gto2.get_m() + gto2.get_n() ) + 3.0 ) *
                   this->overlap_gto(gto1, gto2);

    double term1 = -2.0 * pow(gto2.get_alpha(), 2.0) * (
        this->overlap(gto1.get_alpha(), gto1.get_l(), gto1.get_m(), gto1.get_n(), gto1.get_position(), gto2.get_alpha(), gto2.get_l()+2, gto2.get_m(), gto2.get_n(), gto2.get_position()) +
        this->overlap(gto1.get_alpha(), gto1.get_l(), gto1.get_m(), gto1.get_n(), gto1.get_position(), gto2.get_alpha(), gto2.get_l(), gto2.get_m()+2, gto2.get_n(), gto2.get_position()) +
        this->overlap(gto1.get_alpha(), gto1.get_l(), gto1.get_m(), gto1.get_n(), gto1.get_position(), gto2.get_alpha(), gto2.get_l(), gto2.get_m(), gto2.get_n()+2, gto2.get_position())
    );
    double term2 = -0.5 * (gto2.get_l() * (gto2.get_l() - 1) *
        this->overlap(gto1.get_alpha(), gto1.get_l(), gto1.get_m(), gto1.get_n(), gto1.get_position(), gto2.get_alpha(), gto2.get_l()-2, gto2.get_m(), gto2.get_n(), gto2.get_position()) +
                                                 gto2.get_m() * (gto2.get_m() - 1) *
        this->overlap(gto1.get_alpha(), gto1.get_l(), gto1.get_m(), gto1.get_n(), gto1.get_position(), gto2.get_alpha(), gto2.get_l(), gto2.get_m()-2, gto2.get_n(), gto2.get_position()) +
                                                 gto2.get_n() * (gto2.get_n() - 1) *
        this->overlap(gto1.get_alpha(), gto1.get_l(), gto1.get_m(), gto1.get_n(), gto1.get_position(), gto2.get_alpha(), gto2.get_l(), gto2.get_m(), gto2.get_n()-2, gto2.get_position()) );

    return term0 + term1 + term2;
}

/**
 * @brief Calculate derivative of kinetic integral of two CGFs.
 *
 * @param cgf1     Contracted Gaussian Function 1
 * @param cgf2     Contracted Gaussian Function 2
 * @param nucleus  nucleus position
 * @param coord    Cartesian direction (0=x,1=y,2=z)
 *
 * Calculates the value of d/dcx < cgf1 | -1/2 nabla^2 | cgf2 >
 *
 * @return value of the kinetic integral derivative
 */
double Integrator::kinetic_deriv(const CGF& cgf1,
                                 const CGF& cgf2,
                                 const Vec3& nucleus,
                                 unsigned int coord) const {
    return this->two_center_deriv(prepare_cgf_derivative(cgf1, nucleus), prepare_cgf_derivative(cgf2, nucleus),
                                  coord, &Integrator::kinetic, &Integrator::kinetic_deriv_gto);
}

/**
 * @brief Calculate derivative of kinetic integral of two GTOs.
 *
 * @param gto1   Gaussian Type Orbital 1
 * @param gto2   Gaussian Type Orbital 2
 * @param coord  Cartesian direction (0=x,1=y,2=z)
 *
 * Calculates the value of < d/dx gto1 |-1/2 nabla^2 | gto2 >
 *
 * @return value of the kinetic integral derivative
 */
double Integrator::kinetic_deriv_gto(const GTO& gto1, const GTO& gto2, unsigned int coord) const {
    std::array<unsigned int, 3> gto_ang = {gto1.get_l(), gto1.get_m(), gto1.get_n()};
    if(gto_ang[coord] != 0) {
        gto_ang[coord] += 1; // calculate l+1 term
        double term_plus = this->kinetic_gto(GTO(gto1.get_coefficient(),
                                             gto1.get_position()[0],
                                             gto1.get_position()[1],
                                             gto1.get_position()[2],
                                             gto1.get_alpha(),
                                             gto_ang[0],
                                             gto_ang[1],
                                             gto_ang[2]),
                                         gto2);

        gto_ang[coord] -= 2; // calculate l-1 term
        double term_min = this->kinetic_gto(GTO(gto1.get_coefficient(),
                                             gto1.get_position()[0],
                                             gto1.get_position()[1],
                                             gto1.get_position()[2],
                                             gto1.get_alpha(),
                                             gto_ang[0],
                                             gto_ang[1],
                                             gto_ang[2]),
                                         gto2);
        gto_ang[coord] += 1; // recover l

        return 2.0 * gto1.get_alpha() * term_plus - gto_ang[coord] * term_min;
    } else { // s-type GTO
        gto_ang[coord] += 1;
        double term1 = this->kinetic_gto(GTO(gto1.get_coefficient(),
                                         gto1.get_position()[0],
                                         gto1.get_position()[1],
                                         gto1.get_position()[2],
                                         gto1.get_alpha(),
                                         gto_ang[0],
                                         gto_ang[1],
                                         gto_ang[2]),
                                         gto2);
        return 2.0 * gto1.get_alpha() * term1;
    }
}

/**************************************************************************
 * NUCLEAR INTEGRALS
 **************************************************************************/

/**
 * @brief Calculate nuclear attraction integral of two CGFs.
 *
 * @param cgf1     Contracted Gaussian Function 1
 * @param cgf2     Contracted Gaussian Function 2
 * @param nucleus  nucleus position
 * @param charge   charge of the nucleus
 *
 * Calculates the value of < cgf1 | V | cgf2 >
 *
 * @return value of the nuclear integral
 */
double Integrator::nuclear(const CGF& cgf1,
                           const CGF& cgf2,
                           const Vec3 &nucleus,
                           unsigned int charge) const {
    double sum = 0.0;

    for(unsigned int k = 0; k < cgf1.size(); k++) {
        for(unsigned int l = 0; l < cgf2.size(); l++) {
            sum += cgf1.get_norm_gto(k) *
                   cgf2.get_norm_gto(l) *
                   cgf1.get_coefficient_gto(k) *
                   cgf2.get_coefficient_gto(l) *
                   this->nuclear_gto(cgf1.get_gto(k), cgf2.get_gto(l), nucleus);
        }
    }

    const auto norm1 = cgf1.get_contraction_norm();
    const auto norm2 = cgf2.get_contraction_norm();

    return sum * static_cast<double>(charge) * norm1 * norm2;
}

/**
 * @brief Calculate nuclear attraction integral of two GTOs.
 *
 * @param gto1     Gaussian Type Orbital 1
 * @param gto2     Gaussian Type Orbital 2
 * @param nucleus  nucleus position
 *
 * Calculates the value of < gto1 | V | gto2 >
 *
 * @return value of the nuclear integral
 */
double Integrator::nuclear_gto(const GTO& gto1,
                               const GTO& gto2,
                               const Vec3 &nucleus) const {
    return nuclear(gto1.get_position(), gto1.get_l(), gto1.get_m(), gto1.get_n(), gto1.get_alpha(),
                   gto2.get_position(), gto2.get_l(), gto2.get_m(), gto2.get_n(), gto2.get_alpha(), nucleus);
}

/**
 * @brief Calculate derivative of the nuclear attraction integral for two CGFs.
 *
 * @param cgf1     Contracted Gaussian Function
 * @param cgf2     Contracted Gaussian Function
 * @param nucleus  Position of the nucleus generating the potential
 * @param charge   Charge of the nucleus in a.u.
 * @param nucderiv Center with respect to which the derivative is taken
 * @param coord    Cartesian direction of the derivative (0=x,1=y,2=z)
 *
 * Calculates the value of
 *
 *     d/dR ⟨ cgf1 | V_nuc | cgf2 ⟩
 *
 * where R is the selected Cartesian coordinate of the given center.
 *
 * @return value of the nuclear attraction gradient
 */
double Integrator::nuclear_deriv(const CGF& cgf1, const CGF& cgf2, const Vec3& nucleus,
    unsigned int charge, const Vec3& nucderiv, unsigned int coord) const {
    return this->nuclear_deriv(prepare_cgf_derivative(cgf1, nucderiv), prepare_cgf_derivative(cgf2, nucderiv),
                               nucleus, charge, same_center(nucleus, nucderiv), coord);
}

double Integrator::nuclear_deriv(const CgfDerivative& bra, const CgfDerivative& ket, const Vec3& nucleus,
                                 unsigned int charge, bool nucleus_moves, unsigned int coord) const {
    if ((!nucleus_moves && bra.none_move() && ket.none_move()) ||
        (nucleus_moves && bra.all_move() && ket.all_move())) {
        return 0.0;
    }

    double sum = 0.0;
    for (unsigned int i = 0; i < bra.cgf.size(); ++i) {
        const GTO& g1 = bra.cgf.get_gto(i);
        const bool bra_moves_relative = bra.moves(i) != nucleus_moves;
        const double weight1 = g1.get_norm() * g1.get_coefficient();

        for (unsigned int j = 0; j < ket.cgf.size(); ++j) {
            const GTO& g2 = ket.cgf.get_gto(j);
            double deriv = 0.0;
            if (bra_moves_relative) {
                deriv += this->nuclear_deriv(g1, g2, nucleus, coord);
            }
            if (ket.moves(j) != nucleus_moves) {
                deriv += this->nuclear_deriv(g2, g1, nucleus, coord);
            }
            sum += weight1 * g2.get_norm() * g2.get_coefficient() * deriv;
        }
    }

    // In the potential nucleus's frame, stationary primitives move backwards
    // when that nucleus is displaced.
    const double signed_charge = nucleus_moves ? -static_cast<double>(charge) : static_cast<double>(charge);
    double result = signed_charge * sum * bra.cgf.get_contraction_norm() * ket.cgf.get_contraction_norm();
    const double log_norm_deriv = bra.log_norm_deriv[coord] + ket.log_norm_deriv[coord];
    if (log_norm_deriv != 0.0) {
        result += log_norm_deriv * this->nuclear(bra.cgf, ket.cgf, nucleus, charge);
    }
    return result;
}

/**
 * @brief Calculate derivative of the nuclear attraction integral for two GTOs.
 *
 * @param gto1     Gaussian Type Orbital
 * @param gto2     Gaussian Type Orbital
 * @param nucleus  Position of the nucleus generating the potential
 * @param coord    Cartesian direction of the derivative (0=x,1=y,2=z)
 *
 * Calculates the value of
 *
 *     d/dR ⟨ gto1 | V_nuc | gto2 ⟩
 *
 * where R is the selected Cartesian coordinate of the nuclear position.
 *
 * This routine evaluates the derivative of the nuclear attraction integral and
 * returns the corresponding energy gradient component.
 *
 * @return value of the nuclear attraction gradient
 */
double Integrator::nuclear_deriv(const GTO& gto1, const GTO& gto2,
                                 const Vec3& nucleus, unsigned int coord) const {
    std::array<unsigned int,3> ang = {
        gto1.get_l(),
        gto1.get_m(),
        gto1.get_n()
    };

    const unsigned int l = ang[coord];

    // (l + 1) term
    ang[coord] += 1;
    GTO gto_plus(
        gto1.get_coefficient(),
        gto1.get_position()[0],
        gto1.get_position()[1],
        gto1.get_position()[2],
        gto1.get_alpha(),
        ang[0], ang[1], ang[2]
    );

    const double term_plus =
        this->nuclear_gto(gto_plus, gto2, nucleus);

    // restore l+1 -> l
    ang[coord] -= 1;

    if (l > 0) {
        // (l - 1) term
        ang[coord] -= 1;
        GTO gto_min(
            gto1.get_coefficient(),
            gto1.get_position()[0],
            gto1.get_position()[1],
            gto1.get_position()[2],
            gto1.get_alpha(),
            ang[0], ang[1], ang[2]
        );

        const double term_min = this->nuclear_gto(gto_min, gto2, nucleus);

        // restore l
        ang[coord] += 1;

        return 2.0 * gto1.get_alpha() * term_plus - static_cast<double>(l) * term_min;
    }
    else {
        // s-type along this axis
        return 2.0 * gto1.get_alpha() * term_plus;
    }
}

/**
 * @brief Perform overlap integral evaluation for two primitive GTOs.
 *
 * @param alpha1  Gaussian exponent of the first GTO
 * @param l1      Power of x component of the polynomial of the first GTO
 * @param m1      Power of y component of the polynomial of the first GTO
 * @param n1      Power of z component of the polynomial of the first GTO
 * @param a       Center of the first GTO
 * @param alpha2  Gaussian exponent of the second GTO
 * @param l2      Power of x component of the polynomial of the second GTO
 * @param m2      Power of y component of the polynomial of the second GTO
 * @param n2      Power of z component of the polynomial of the second GTO
 * @param b       Center of the second GTO
 *
 * Calculates the value of < gto1 | gto2 >
 *
 * @return value of the overlap integral
 */
double Integrator::overlap(double alpha1, unsigned int l1, unsigned int m1, unsigned int n1, const Vec3 &a,
                           double alpha2, unsigned int l2, unsigned int m2, unsigned int n2, const Vec3 &b) const {
    return integrals::gaussian::overlap_gto(alpha1, l1, m1, n1, a, alpha2, l2, m2, n2, b);
}

/**
 * @brief Perform dipole integral evaluation for two primitive GTOs.
 *
 * @param alpha1  Gaussian exponent of the first GTO
 * @param l1      Power of x component of the polynomial of the first GTO
 * @param m1      Power of y component of the polynomial of the first GTO
 * @param n1      Power of z component of the polynomial of the first GTO
 * @param a       Center of the first GTO
 * @param alpha2  Gaussian exponent of the second GTO
 * @param l2      Power of x component of the polynomial of the second GTO
 * @param m2      Power of y component of the polynomial of the second GTO
 * @param n2      Power of z component of the polynomial of the second GTO
 * @param b       Center of the second GTO
 * @param cc      Cartesian direction (0=x,1=y,2=z)
 * @param cref    Reference position
 *
 * @return value of the dipole integral
 */
double Integrator::dipole(double alpha1, unsigned int l1, unsigned int m1, unsigned int n1, const Vec3 &a,
                          double alpha2, unsigned int l2, unsigned int m2, unsigned int n2, const Vec3 &b,
                          unsigned int cc, double cref) const {

    double rab2 = (a-b).norm2();
    double gamma = alpha1 + alpha2;

    // determine new product center
    Vec3 p = this->gaussian_product_center(alpha1, a, alpha2, b);

    // determine correcting pre-factor
    double pre = std::pow(M_PI / gamma, 1.5) * std::exp(-alpha1 * alpha2 * rab2 / gamma);

    // construct adjusted triple product
    Vec3 w;
    const std::array<unsigned int,3> o1({l1,m1,n1});
    const std::array<unsigned int,3> o2({l2,m2,n2});

    for(unsigned int i=0; i<3; i++) {
        w[i] = this->overlap_1D(o1[i], o2[i], p[i]-a[i], p[i]-b[i], gamma);
    }

    // create copy and adjust the "cc" value
    Vec3 wd = w;
    wd[cc] = this->overlap_1D(o1[cc], o2[cc]+1,
                             p[cc]-a[cc],
                             p[cc]-b[cc], gamma);

    return pre * (wd.prod() + (b[cc] - cref) * w.prod());
}

/**
 * @brief Calculate one-dimensional overlap integral.
 *
 * @param l1     Power of 'x' component of the polynomial of the first GTO
 * @param l2     Power of 'x' component of the polynomial of the second GTO
 * @param x1     'x' component of the position of the first GTO
 * @param x2     'x' component of the position of the second GTO
 * @param gamma  Sum of the two Gaussian exponents
 *
 * NOTE: in contrast to other places, here int has to be used rather than unsigned int
 *       because sometimes negative numbers are parsed
 *
 * @return value of the one dimensional overlap integral
 */
double Integrator::overlap_1D(int l1, int l2, double x1, double x2, double gamma) const {
    return integrals::gaussian::overlap_1d(l1, l2, x1, x2, gamma);
}

/**
 * @brief Calculates the Gaussian product center of two GTOs
 *
 * @param alpha1  Gaussian exponent of the first GTO
 * @param a       Center of the first GTO
 * @param alpha2  Gaussian exponent of the second GTO
 * @param b       Center of the second GTO
 *
 *
 * @return new gaussian product center
 */
Vec3 Integrator::gaussian_product_center(double alpha1, const Vec3& a,
                                         double alpha2, const Vec3& b) const {
    return integrals::gaussian::gaussian_product_center(alpha1, a, alpha2, b);
}

/**
 * @brief Compute the binomial prefactor for Hermite Gaussian recursion.
 *
 * @param s    Summation index
 * @param ia   Angular momentum component of first Gaussian
 * @param ib   Angular momentum component of second Gaussian
 * @param xpa  Distance P_x - A_x
 * @param xpb  Distance P_x - B_x
 *
 * @return binomial prefactor value
 */
double Integrator::binomial_prefactor(int s, int ia, int ib,
                                      double xpa, double xpb) const {
    return integrals::gaussian::binomial_prefactor(s, ia, ib, xpa, xpb);
}

/**
 * @brief Compute the binomial coefficient.
 *
 * @param a  Upper index
 * @param b  Lower index
 *
 * @return binomial coefficient
 */
double Integrator::binomial(int a, int b) const {
    return integrals::gaussian::binomial(a, b);
}

/**
 * @brief Perform nuclear attraction integral evaluation for primitive GTOs.
 *
 * @param a       Center of the first GTO
 * @param l1      Power of x component of the polynomial of the first GTO
 * @param m1      Power of y component of the polynomial of the first GTO
 * @param n1      Power of z component of the polynomial of the first GTO
 * @param alpha1  Gaussian exponent of the first GTO
 * @param b       Center of the second GTO
 * @param l2      Power of x component of the polynomial of the second GTO
 * @param m2      Power of y component of the polynomial of the second GTO
 * @param n2      Power of z component of the polynomial of the second GTO
 * @param alpha2  Gaussian exponent of the second GTO
 * @param c       Nuclear position
 *
 * @return value of the nuclear integral
 */
double Integrator::nuclear(const Vec3& a, int l1, int m1, int n1, double alpha1,
                           const Vec3& b, int l2, int m2, int n2,
                           double alpha2, const Vec3& c) const {


    double gamma = alpha1 + alpha2;

    Vec3 p = gaussian_product_center(alpha1, a, alpha2, b);
    double rab2 = (a-b).norm2();
    double rcp2 = (c-p).norm2();

    std::vector<double> ax = A_array(l1, l2, p[0]-a[0], p[0]-b[0], p[0]-c[0], gamma);
    std::vector<double> ay = A_array(m1, m2, p[1]-a[1], p[1]-b[1], p[1]-c[1], gamma);
    std::vector<double> az = A_array(n1, n2, p[2]-a[2], p[2]-b[2], p[2]-c[2], gamma);

    double sum = 0.0;

    for(int i=0; i<=l1+l2;i++) {
        for(int j=0; j<=m1+m2;j++) {
            for(int k=0; k<=n1+n2;k++) {
                sum += ax[i] * ay[j] * az[k] * this->boys_function(i+j+k,rcp2*gamma);
            }
        }
    }

    return -2.0 * M_PI / gamma * std::exp(-alpha1*alpha2*rab2/gamma) * sum;
}

/**
 * @brief Build the Hermite A array.
 *
 * @param l1  Angular momentum component for center A
 * @param l2  Angular momentum component for center B
 * @param pa  P-A distance component
 * @param pb  P-B distance component
 * @param cp  C-P distance component
 * @param g   Gaussian exponent sum
 *
 * @return Hermite A array
 */
std::vector<double> Integrator::A_array(const int l1, const int l2, const double pa, const double pb, const double cp, const double g) const {
    int imax = l1 + l2 + 1;
    std::vector<double> arrA(imax, 0);

    for(int i=0; i<imax; i++) {
        for(int r=0; r<=i/2; r++) {
            for(int u=0; u<=(i-2*r)/2; u++) {
                int iI = i - 2 * r - u;
                arrA[iI] += A_term(i, r, u, l1, l2, pa, pb, cp, g);
            }
        }
    }

    return arrA;
}

/**
 * @brief Evaluate a single Hermite A term.
 *
 * @param i      Hermite index
 * @param r      recursion index
 * @param u      Boys index
 * @param l1     Angular momentum component for center A
 * @param l2     Angular momentum component for center B
 * @param pax    P-A distance component
 * @param pbx    P-B distance component
 * @param cpx    C-P distance component
 * @param gamma  Gaussian exponent sum
 *
 * @return Hermite A term value
 */
double Integrator::A_term(const int i, const int r, const int u, const int l1, const int l2, const double pax, const double pbx, const double cpx, const double gamma) const {
    return  std::pow(-1,i) * this->binomial_prefactor(i,l1,l2,pax,pbx)*
            std::pow(-1,u) * factorial(i)*std::pow(cpx,i-2*r-2*u)*
            std::pow(0.25/gamma,r+u)/factorial(r)/factorial(u)/factorial(i-2*r-2*u);
}
