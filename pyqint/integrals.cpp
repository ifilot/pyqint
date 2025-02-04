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

/**
 * @brief Integrator constructor method
 *
 * @return Integrator class
 */
Integrator::Integrator(){}

/**
 * @brief      Evaluate all integrals for cgfs in buffer
 */
std::vector<double> Integrator::evaluate_cgfs(const std::vector<CGF>& cgfs,
                                              const std::vector<int>& charges,
                                              const std::vector<double>& px,
                                              const std::vector<double>& py,
                                              const std::vector<double>& pz) const {
    std::vector<double> results;

    size_t sz = cgfs.size();

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

    // calculate all two-electron integrals
    std::vector<double> tedouble(this->teindex(sz-1,sz-1,sz-1,sz-1) + 1, -1.0);

    // it is more efficient to first 'unroll' the fourfold nested loop
    // into a single vector of jobs to execute
    std::vector<std::array<size_t, 5>> jobs;
    for(size_t i=0; i<sz; i++) {
        for(size_t j=0; j<sz; j++) {
            size_t ij = i*(i+1)/2 + j;
            for(size_t k=0; k<sz; k++) {
                for(size_t l=0; l<sz; l++) {
                    size_t kl = k * (k+1)/2 + l;
                    if(ij <= kl) {
                        const size_t idx = this->teindex(i,j,k,l);

                        if(idx >= tedouble.size()) {
                            throw std::runtime_error("Process tried to access illegal array position");
                        }

                        if(tedouble[idx] < 0.0) {
                            tedouble[idx] = 1.0;
                            jobs.emplace_back(std::array<size_t,5>({idx, i, j, k, l}));
                        }
                    }
                }
            }
        }
    }

    // evaluate jobs
    #pragma omp parallel for schedule(dynamic)
    for(int s=0; s<(int)jobs.size(); s++) {  // have to use signed int for MSVC OpenMP here
        const size_t idx = jobs[s][0];
        const size_t i = jobs[s][1];
        const size_t j = jobs[s][2];
        const size_t k = jobs[s][3];
        const size_t l = jobs[s][4];
        tedouble[idx] = this->repulsion(cgfs[i], cgfs[j], cgfs[k], cgfs[l]);
    }

    // reorganize everyting into a tensor object
    std::vector<double> tetensor(std::pow(sz,4));
    for(size_t i=0; i<sz; i++) {
        for(size_t j=0; j<sz; j++) {
            for(size_t k=0; k<sz; k++) {
                for(size_t l=0; l<sz; l++) {
                    const size_t idx = this->teindex(i,j,k,l);
                    tetensor[i*sz*sz*sz + j*sz*sz + k*sz + l] = tedouble[idx];
                }
            }
        }
    }

    // package everything into results vector, will be unpacked in
    // connected Python class
    std::vector<double> Svec(S.data(), S.data()+sz*sz);
    results.insert(results.end(), Svec.begin(), Svec.end());

    // kinetic integrals
    std::vector<double> Tvec(T.data(), T.data() + sz*sz);
    results.insert(results.end(), Tvec.begin(), Tvec.end());

    // nuclear integrals
    std::vector<double> Vvec(V.data(), V.data() + sz*sz);
    results.insert(results.end(), Vvec.begin(), Vvec.end());

    // two electron in tegrals
    results.insert(results.end(), tetensor.begin(), tetensor.end());

    return results;
}

/**
 * @brief      Evaluate the geometric derivates for all cgfs in buffer
 */
std::vector<double> Integrator::evaluate_geometric_derivatives(const std::vector<CGF>& cgfs,
                                                               const std::vector<int>& charges,
                                                               const std::vector<double>& px,
                                                               const std::vector<double>& py,
                                                               const std::vector<double>& pz) const {
    const size_t sz = cgfs.size();

    // Construct 2x2 matrices to hold values for the overlap,
    // kinetic and two nuclear integral values, respectively.
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
                    S[idx] = this->overlap_deriv(cgfs[i], cgfs[j], Vec3(px[n], py[n], pz[n]), k);
                    T[idx] = this->kinetic_deriv(cgfs[i], cgfs[j], Vec3(px[n], py[n], pz[n]), k);


                    for(unsigned int l=0; l<charges.size(); l++) {
                        Vn[l*sz*sz + i*sz + j] = this->nuclear_deriv(cgfs[i], cgfs[j],
                                                         Vec3(px[l], py[l], pz[l]),
                                                         charges[l],
                                                         Vec3(px[n], py[n], pz[n]),
                                                         k);
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
    const unsigned int max_teints = this->teindex(sz-1,sz-1,sz-1,sz-1) + 1;
    std::vector<double> tedouble(max_teints * charges.size() * 3, -1.0);

    // it is more efficient to first 'unroll' the fourfold nested loop
    // into a single vector of jobs to execute
    std::vector<std::array<size_t, 7>> jobs;
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

    // // evaluate jobs
    #pragma omp parallel for schedule(dynamic)
    for(int s=0; s<(int)jobs.size(); s++) {  // have to use signed int for MSVC OpenMP here
        const size_t idx = jobs[s][0];
        const size_t i = jobs[s][1];
        const size_t j = jobs[s][2];
        const size_t k = jobs[s][3];
        const size_t l = jobs[s][4];
        const size_t n = jobs[s][5];
        const size_t d = jobs[s][6];
        tedouble[idx] = this->repulsion_deriv(cgfs[i], cgfs[j], cgfs[k], cgfs[l], Vec3(px[n], py[n], pz[n]), d);
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
 * @brief Calculates overlap integral of two CGF
 *
 * @param const CGF& cgf1   Contracted Gaussian Function
 * @param const CGF& cgf2   Contracted Gaussian Function
 *
 * Calculates the value of < cgf1 | cgf2 >
 *
 * @return double value of the overlap integral
 */
double Integrator::overlap(const CGF& cgf1, const CGF& cgf2) const {
    double sum = 0.0;

    // loop over all GTOs inside the CGF, calculate the overlap integrals
    // and sum all the integral values
    for(unsigned int k = 0; k < cgf1.size(); k++) {
        for(unsigned int l = 0; l < cgf2.size(); l++) {
            sum += cgf1.get_norm_gto(k) *
                   cgf2.get_norm_gto(l) *
                   cgf1.get_coefficient_gto(k) *
                   cgf2.get_coefficient_gto(l) *
                   this->overlap_gto(cgf1.get_gto(k), cgf2.get_gto(l) );
        }
    }

    return sum;
}

/**
 * @brief      Calculates overlap integral of two GTO
 *
 * @param[in]  gto1  Gaussian Type Orbital 1
 * @param[in]  gto2  Gaussian Type Orbital 2
 *
 * Calculates the value of < gto1 | gto2 >
 *
 * @return     double value of the overlap integral
 */
double Integrator::overlap_gto(const GTO& gto1, const GTO& gto2) const {
    return this->overlap(gto1.get_alpha(), gto1.get_l(), gto1.get_m(), gto1.get_n(), gto1.get_position(),
                         gto2.get_alpha(), gto2.get_l(), gto2.get_m(), gto2.get_n(), gto2.get_position());
}

/**
 * @brief      Calculates the geometric derivative of overlap integral of
 *             two CGF
 *
 * @param[in]  cgf1     Gaussian Contracted Functional 1
 * @param[in]  cgf2     Gaussian Contracted Functional 2
 * @param[in]  nucleus  Nucleus position
 * @param[in]  coord    Cartesian direction
 *
 * Calculates the value of d/dx < cgf1 | cgf2 >
 *
 * @return     double value of the nuclear integral
 */
double Integrator::overlap_deriv(const CGF& cgf1,
                                 const CGF& cgf2,
                                 const Vec3& nucleus,
                                 unsigned int coord) const {
    double sum = 0.0;

    // check if cgf originates from nucleus
    bool cgf1_nuc = (cgf1.get_r() - nucleus).norm2() < 0.0001;
    bool cgf2_nuc = (cgf2.get_r() - nucleus).norm2() < 0.0001;

    if(cgf1_nuc == cgf2_nuc) { // early exit
        return 0.0;
    }

    for(unsigned int k = 0; k < cgf1.size(); k++) {
        for(unsigned int l = 0; l < cgf2.size(); l++) {

            double t1 = cgf1_nuc ? this->overlap_deriv_gto(cgf1.get_gto(k),
                                                           cgf2.get_gto(l),
                                                           coord) : 0.0;
            double t2 = cgf2_nuc ? this->overlap_deriv_gto(cgf2.get_gto(l),
                                                           cgf1.get_gto(k),
                                                           coord) : 0.0;

            sum += cgf1.get_norm_gto(k) *
                   cgf2.get_norm_gto(l) *
                   cgf1.get_coefficient_gto(k) *
                   cgf2.get_coefficient_gto(l) *
                   (t1 + t2);
        }
    }

    return sum;
}

/**
 * @brief Calculates overlap integral of two GTO
 *
 * @param const GTO& gto1   Gaussian Type Orbital
 * @param const GTO& gto2   Gaussian Type Orbital
 *
 * Calculates the value of < gto1 | gto2 >
 *
 * @return double value of the overlap integral
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

/**
 * @brief      Calculates dipole integral of two CGF
 *
 * @param[in]  cgf1  Contracted Gaussian Function
 * @param[in]  cgf2  Contracted Gaussian Function
 * @param[in]  cc    Cartesian direction
 * @param[in]  ref   Reference position
 *
 * Calculates the value of < cgf1 | (cc - ref) | cgf2 >
 *
 * @return     value of the dipole integral
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

    return sum;
}

/**
 * @brief      Calculates dipole integral of two GTO
 *
 * @param[in]  gto1   Gaussian Type Orbital
 * @param[in]  gto2   Gaussian Type Orbital
 * @param      cc     Cartesian direction
 * @param[in]  ref    The reference
 *
 * Calculates the value of < gto1 | (cc - ref) | gto2 >
 *
 * @return     value of the dipole integral
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
 * @brief      Calculates kinetic integral of two CGF
 *
 * @param[in]  cgf1  Contracted Gaussian Functional 1
 * @param[in]  cgf2  Contracted Gaussian Functional 2
 *
 * Calculates the value of < cgf1 | T | cgf2 >
 *
 * @return     double value of the kinetic integral
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

    return sum;
}

/**
 * @brief      Calculates kinetic integral of two GTO
 *
 * @param[in]  gto1  Gaussian Type Orbital 1
 * @param[in]  gto2  Gaussian Type Orbital 2
 *
 * Calculates the value of < gto1 | T | gto2 >
 *
 * @return     double value of the kinetic integral
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
 * @brief      Calculates derivative of kinetic integral of two CGF
 *
 * @param[in]  cgf1     Contracted Gaussian Functional 1
 * @param[in]  cgf2     Contracted Gaussian Functional 2
 * @param[in]  nucleus  nucleus position
 * @param[in]  coord    derivative coordinate
 *
 * Calculates the value of d/dcx < cgf1 | -1/2 nabla^2 | cgf2 >
 *
 * @return     double value of the kinetic integral
 */
double Integrator::kinetic_deriv(const CGF& cgf1,
                                 const CGF& cgf2,
                                 const Vec3& nucleus,
                                 unsigned int coord) const {
    double sum = 0.0;

    // check if cgf originates from nucleus
    bool cgf1_nuc = (cgf1.get_r() - nucleus).norm2() < 0.0001;
    bool cgf2_nuc = (cgf2.get_r() - nucleus).norm2() < 0.0001;

    if(cgf1_nuc == cgf2_nuc) { // early exit
        return 0.0;
    }

    for(unsigned int k = 0; k < cgf1.size(); k++) {
        for(unsigned int l = 0; l < cgf2.size(); l++) {

            double t1 = cgf1_nuc ? this->kinetic_deriv_gto(cgf1.get_gto(k), cgf2.get_gto(l), coord) : 0.0;
            double t2 = cgf2_nuc ? this->kinetic_deriv_gto(cgf2.get_gto(l), cgf1.get_gto(k), coord) : 0.0;

            sum += cgf1.get_norm_gto(k) *
                   cgf2.get_norm_gto(l) *
                   cgf1.get_coefficient_gto(k) *
                   cgf2.get_coefficient_gto(l) *
                   (t1 + t2);
        }
    }

    return sum;
}

/**
 * @brief      Calculates derivative of kinetic integral of two GTOs
 *
 * @param[in]  gto1   The gto 1
 * @param[in]  gto2   The gto 2
 * @param      unsigned  int coord    Derivative direction
 *
 * Calculates the value of < d/dx gto1 |-1/2 nabla^2 | gto2 >
 *
 * @return     double value of the derivative of the kinetic integral
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
 * @brief      Calculates nuclear integral of two CGF
 *
 * @param[in]  cgf1     Contracted Gaussian Functional 1
 * @param[in]  cgf2     Contracted Gaussian Functional 2
 * @param[in]  nucleus  nucleus position
 * @param[in]  charge   charge of the nucleus
 *
 * Calculates the value of < cgf1 | V | cgf2 >
 *
 * @return     double value of the nuclear integral
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

    return sum * (double)charge;
}

/**
 * @brief      Calculates nuclear integral of two GTOs
 *
 * @param[in]  gto1     Gaussian Type Orbital 1
 * @param[in]  gto2     Gaussian Type Orbital 2
 * @param[in]  nucleus  nucleus position
 *
 * Calculates the value of < gto1 | V | gto2 >
 *
 * @return     double value of the nuclear integral
 */
double Integrator::nuclear_gto(const GTO& gto1,
                               const GTO& gto2,
                               const Vec3 &nucleus) const {
    return nuclear(gto1.get_position(), gto1.get_l(), gto1.get_m(), gto1.get_n(), gto1.get_alpha(),
                   gto2.get_position(), gto2.get_l(), gto2.get_m(), gto2.get_n(), gto2.get_alpha(), nucleus);
}

/**
 * @brief Calculates nuclear integral of two CGF
 *
 * @param const CGF& cgf1       Contracted Gaussian Function
 * @param const CGF& cgf2       Contracted Gaussian Function
 * @param unsigned int charge   charge of the nucleus in a.u.
 *
 * Calculates the value of < cgf1 | V | cgf2 >
 *
 * In contrast to the other derivatives, this will be done using a finite
 * difference procedure.
 *
 * @return double value of the nuclear integral
 */
double Integrator::nuclear_deriv(const CGF& cgf1, const CGF& cgf2, const Vec3 &nucleus, unsigned int charge,
                                 const Vec3& nucderiv, unsigned int coord) const {

    // check if cgf originates from nucleus
    bool n1 = (cgf1.get_r() - nucderiv).norm2() < 0.0001;
    bool n2 = (cgf2.get_r() - nucderiv).norm2() < 0.0001;
    bool n3 = (nucleus - nucderiv).norm2() < 0.0001;

    if(n1 == n2 && n2 == n3) {
        return 0.0;
    }

    static const double delta = 1e-5;

    // create copies of objects
    CGF cgf1m = cgf1;
    CGF cgf1p = cgf1;
    CGF cgf2m = cgf2;
    CGF cgf2p = cgf2;
    Vec3 nucm = nucleus;
    Vec3 nucp = nucleus;

    if(n1) {
        Vec3 rm = cgf1.get_r();
        rm[coord] -= delta;
        cgf1m.set_position(rm);
        Vec3 rp = cgf1.get_r();
        rp[coord] += delta;
        cgf1p.set_position(rp);
    }

    if(n2) {
        Vec3 rm = cgf2.get_r();
        rm[coord] -= delta;
        cgf2m.set_position(rm);
        Vec3 rp = cgf2.get_r();
        rp[coord] += delta;
        cgf2p.set_position(rp);
    }

    if(n3) {
        nucm[coord] -= delta;
        nucp[coord] += delta;
    }

    double vm = this->nuclear(cgf1m, cgf2m, nucm, charge);
    double vp = this->nuclear(cgf1p, cgf2p, nucp, charge);

    return (vp - vm) / (2.0 * delta);
}

/**************************************************************************
 * TWO-ELECTRON INTEGRALS
 **************************************************************************/

/**
 * @brief Calculates two-electron repulsion integral of four CGF
 *
 * @param const CGF& cgf1       Contracted Gaussian Function
 * @param const CGF& cgf2       Contracted Gaussian Function
 * @param const CGF& cgf3       Contracted Gaussian Function
 * @param const CGF& cgf4       Contracted Gaussian Function
 *
 * Calculates the value of < cgf1 | cgf2 | cgf3 | cgf4 >
 *
 * @return double value of the repulsion integral
 */
double Integrator::repulsion(const CGF &cgf1,const CGF &cgf2,const CGF &cgf3,const CGF &cgf4) const {
    double sum = 0;

    for(unsigned int i=0; i< cgf1.size(); i++) {
        for(unsigned int j=0; j< cgf2.size(); j++) {
            for(unsigned int k=0; k < cgf3.size(); k++) {
                for(unsigned int l=0; l < cgf4.size(); l++) {
                    const double n1 = cgf1.get_norm_gto(i);
                    const double n2 = cgf2.get_norm_gto(j);
                    const double n3 = cgf3.get_norm_gto(k);
                    const double n4 = cgf4.get_norm_gto(l);

                    double pre = cgf1.get_coefficient_gto(i) * cgf2.get_coefficient_gto(j) * cgf3.get_coefficient_gto(k) * cgf4.get_coefficient_gto(l);
                    sum += n1 * n2 * n3 * n4 * pre * this->repulsion(cgf1.get_gto(i), cgf2.get_gto(j), cgf3.get_gto(k), cgf4.get_gto(l));
                }
            }
        }
    }

    return sum;
}

/**
 * @brief Calculates derivative of the two-electron repulsion integral of four CGF
 *
 * @param const CGF& cgf1       Contracted Gaussian Function
 * @param const CGF& cgf2       Contracted Gaussian Function
 * @param const CGF& cgf3       Contracted Gaussian Function
 * @param const CGF& cgf4       Contracted Gaussian Function
 * @param const Vec3& nucleus   Nucleus coordinates
 * @param unsigned int coord    Derivative direction
 *
 * Calculates the value of d/dcx < cgf1 | cgf2 | cgf3 | cgf4 >
 *
 * @return double value of the repulsion integral
 */
double Integrator::repulsion_deriv(const CGF &cgf1, const CGF &cgf2, const CGF &cgf3, const CGF &cgf4,
    const Vec3& nucleus, unsigned int coord) const {
    double sum = 0;

    // check if cgf originates from nucleus
    bool cgf1_nuc = (cgf1.get_r() - nucleus).norm2() < 0.0001;
    bool cgf2_nuc = (cgf2.get_r() - nucleus).norm2() < 0.0001;
    bool cgf3_nuc = (cgf3.get_r() - nucleus).norm2() < 0.0001;
    bool cgf4_nuc = (cgf4.get_r() - nucleus).norm2() < 0.0001;

    // early exit
    if(cgf1_nuc == cgf2_nuc && cgf2_nuc == cgf3_nuc && cgf3_nuc == cgf4_nuc) {
        return 0.0;
    }

    for(unsigned int i=0; i< cgf1.size(); i++) {
        for(unsigned int j=0; j< cgf2.size(); j++) {
            for(unsigned int k=0; k < cgf3.size(); k++) {
                for(unsigned int l=0; l < cgf4.size(); l++) {

                    // calculate product of coefficients
                    double pre = cgf1.get_coefficient_gto(i) * cgf2.get_coefficient_gto(j) * cgf3.get_coefficient_gto(k) * cgf4.get_coefficient_gto(l);

                    // get normalization factors
                    const double n1 = cgf1.get_norm_gto(i);
                    const double n2 = cgf2.get_norm_gto(j);
                    const double n3 = cgf3.get_norm_gto(k);
                    const double n4 = cgf4.get_norm_gto(l);

                    // take the derivative towards the basis functions
                    double t1 = cgf1_nuc ? this->repulsion_deriv(cgf1.get_gto(i), cgf2.get_gto(j), cgf3.get_gto(k), cgf4.get_gto(l), coord) : 0.0;
                    double t2 = cgf2_nuc ? this->repulsion_deriv(cgf2.get_gto(j), cgf1.get_gto(i), cgf3.get_gto(k), cgf4.get_gto(l), coord) : 0.0;
                    double t3 = cgf3_nuc ? this->repulsion_deriv(cgf3.get_gto(k), cgf4.get_gto(l), cgf1.get_gto(i), cgf2.get_gto(j), coord) : 0.0;
                    double t4 = cgf4_nuc ? this->repulsion_deriv(cgf4.get_gto(l), cgf3.get_gto(k), cgf1.get_gto(i), cgf2.get_gto(j), coord) : 0.0;

                    sum += pre * n1 * n2 * n3 * n4 * (t1 + t2 + t3 + t4);
                }
            }
        }
    }

    return sum;
}

/**
 * @brief Calculates two-electron repulsion integral of four CGF
 *
 * @param const GTO& gto1       Contracted Gaussian Function
 * @param const GTO& gto2       Contracted Gaussian Function
 * @param const GTO& gto3       Contracted Gaussian Function
 * @param const GTO& gto4       Contracted Gaussian Function
 *
 * Calculates the value of < gto1 | gto2 | gto3 | gto4 >
 *
 * @return double value of the repulsion integral
 */
double Integrator::repulsion(const GTO &gto1, const GTO &gto2, const GTO &gto3, const GTO &gto4) const {

    double rep = this->repulsion(gto1.get_position(), gto1.get_l(), gto1.get_m(), gto1.get_n(), gto1.get_alpha(),
                                 gto2.get_position(), gto2.get_l(), gto2.get_m(), gto2.get_n(), gto2.get_alpha(),
                                 gto3.get_position(), gto3.get_l(), gto3.get_m(), gto3.get_n(), gto3.get_alpha(),
                                 gto4.get_position(), gto4.get_l(), gto4.get_m(), gto4.get_n(), gto4.get_alpha());

    return rep;
}

/**
 * @brief Calculates overlap integral of two GTO
 *
 * @param const GTO& gto1       Gaussian Type Orbital
 * @param const GTO& gto2       Gaussian Type Orbital
 * @param const GTO& gto3       Gaussian Type Orbital
 * @param const GTO& gto4       Gaussian Type Orbital
 * @param unsigned int coord    Derivative direction
 *
 * Calculates the value of < d/dx gto1 | gto2 | gto3 | gto4 >
 *
 * @return double value of the overlap integral
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
 * @brief Performs overlap integral evaluation
 *
 * @param double alpha1     Gaussian exponent of the first GTO
 * @param unsigned int l1   Power of x component of the polynomial of the first GTO
 * @param unsigned int m1   Power of y component of the polynomial of the first GTO
 * @param unsigned int n1   Power of z component of the polynomial of the first GTO
 * @param Vec3 a            Center of the Gaussian orbital of the first GTO
 * @param double alpha2     Gaussian exponent of the second GTO
 * @param unsigned int l2   Power of x component of the polynomial of the second GTO
 * @param unsigned int m2   Power of y component of the polynomial of the second GTO
 * @param unsigned int n2   Power of z component of the polynomial of the second GTO
 * @param Vec3 b            Center of the Gaussian orbital of the second GTO
 *
 * Calculates the value of < gto1 | gto2 >
 *
 * @return double value of the overlap integral
 */
double Integrator::overlap(double alpha1, unsigned int l1, unsigned int m1, unsigned int n1, const Vec3 &a,
                           double alpha2, unsigned int l2, unsigned int m2, unsigned int n2, const Vec3 &b) const {

    double rab2 = (a-b).norm2();
    double gamma = alpha1 + alpha2;
    Vec3 p = this->gaussian_product_center(alpha1, a, alpha2, b);

    double pre = std::pow(M_PI / gamma, 1.5) * std::exp(-alpha1 * alpha2 * rab2 / gamma);
    double wx = this->overlap_1D(l1, l2, p[0]-a[0], p[0]-b[0], gamma);
    double wy = this->overlap_1D(m1, m2, p[1]-a[1], p[1]-b[1], gamma);
    double wz = this->overlap_1D(n1, n2, p[2]-a[2], p[2]-b[2], gamma);

    return pre * wx * wy * wz;
}

/**
 * @brief Performs overlap integral evaluation
 *
 * @param double alpha1     Gaussian exponent of the first GTO
 * @param unsigned int l1   Power of x component of the polynomial of the first GTO
 * @param unsigned int m1   Power of y component of the polynomial of the first GTO
 * @param unsigned int n1   Power of z component of the polynomial of the first GTO
 * @param Vec3 a            Center of the Gaussian orbital of the first GTO
 * @param double alpha2     Gaussian exponent of the second GTO
 * @param unsigned int l2   Power of x component of the polynomial of the second GTO
 * @param unsigned int m2   Power of y component of the polynomial of the second GTO
 * @param unsigned int n2   Power of z component of the polynomial of the second GTO
 * @param Vec3 b            Center of the Gaussian orbital of the second GTO
 *
 * @return double value of the overlap integral
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
 * @brief Calculates one dimensional overlap integral
 *
 * @param int l1        Power of 'x' component of the polynomial of the first GTO
 * @param int l2        Power of 'x' component of the polynomial of the second GTO
 * @param double x1     'x' component of the position of the first GTO
 * @param double x2     'x' component of the position of the second GTO
 * @param double gamma  Sum of the two Gaussian exponents
 *
 * NOTE: in contrast to other places, here int has to be used rather than unsigned int
 *       because sometimes negative numbers are parsed
 *
 * @return double value of the one dimensional overlap integral
 */
double Integrator::overlap_1D(int l1, int l2, double x1, double x2, double gamma) const {
    double sum = 0.0;

    for(int i=0; i < (1 + std::floor(0.5 * (l1 + l2))); i++) {
        sum += this->binomial_prefactor(2*i, l1, l2, x1, x2) *
                     (i == 0 ? 1 : double_factorial(2 * i - 1) ) /
                     std::pow(2 * gamma, i);
    }

    return sum;
}

/**
 * @brief Calculates the Gaussian product center of two GTOs
 *
 * @param double alpha1     Gaussian exponent of the first GTO
 * @param double alpha2     Gaussian exponent of the second GTO
 * @param const Vec3 a      Center of the first GTO
 * @param const Vec3 b      Center of the second GTO
 *
 *
 * @return new gaussian product center
 */
Vec3 Integrator::gaussian_product_center(double alpha1, const Vec3& a,
                                         double alpha2, const Vec3& b) const {
    return (alpha1 * a + alpha2 * b) / (alpha1 + alpha2);
}

double Integrator::binomial_prefactor(int s, int ia, int ib,
                                      double xpa, double xpb) const {
    double sum = 0.0;

    for (int t=0; t < s+1; t++) {
        if ((s-ia <= t) && (t <= ib)) {
            sum += this->binomial(ia, s-t)   *
                   this->binomial(ib, t)     *
                   std::pow(xpa, ia - s + t) *
                   std::pow(xpb, ib - t);
        }
    }

    return sum;
}

double Integrator::binomial(int a, int b) const {
    if( (a < 0) || (b < 0) || (a-b < 0) ) {
        return 1.0;
    }
    return factorial(a) / (factorial(b) * factorial(a-b));
}

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
                sum += ax[i] * ay[j] * az[k] * this->gamma_inc.Fgamma(i+j+k,rcp2*gamma);
            }
        }
    }

    return -2.0 * M_PI / gamma * std::exp(-alpha1*alpha2*rab2/gamma) * sum;
}

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

double Integrator::A_term(const int i, const int r, const int u, const int l1, const int l2, const double pax, const double pbx, const double cpx, const double gamma) const {
    return  std::pow(-1,i) * this->binomial_prefactor(i,l1,l2,pax,pbx)*
            std::pow(-1,u) * factorial(i)*std::pow(cpx,i-2*r-2*u)*
            std::pow(0.25/gamma,r+u)/factorial(r)/factorial(u)/factorial(i-2*r-2*u);
}

/**
 * @brief Performs nuclear integral evaluation
 *
 * @param Vec3 a            Center of the Gaussian orbital of the first GTO
 * @param unsigned int l1   Power of x component of the polynomial of the first GTO
 * @param unsigned int m1   Power of y component of the polynomial of the first GTO
 * @param unsigned int n1   Power of z component of the polynomial of the first GTO
 * @param double alpha1     Gaussian exponent of the first GTO
 * @param Vec3 b            Center of the Gaussian orbital of the second GTO
 * @param unsigned int l2   Power of x component of the polynomial of the second GTO
 * @param unsigned int m2   Power of y component of the polynomial of the second GTO
 * @param unsigned int n2   Power of z component of the polynomial of the second GTO
 * @param double alpha2     Gaussian exponent of the second GTO
 * @param Vec3 c
 *
 * @return double value of the nuclear integral
 */
double Integrator::repulsion(const Vec3 &a, const int la, const int ma, const int na, const double alphaa,
                             const Vec3 &b, const int lb, const int mb, const int nb, const double alphab,
                             const Vec3 &c, const int lc, const int mc, const int nc, const double alphac,
                             const Vec3 &d, const int ld, const int md, const int nd, const double alphad) const {

    double rab2 = (a-b).norm2();
    double rcd2 = (c-d).norm2();

    Vec3 p = gaussian_product_center(alphaa, a, alphab, b);
    Vec3 q = gaussian_product_center(alphac, c, alphad, d);

    double rpq2 = (p-q).norm2();

    double gamma1 = alphaa + alphab;
    double gamma2 = alphac + alphad;
    double delta = 0.25 * (1.0 / gamma1 + 1.0 / gamma2);

    std::vector<double> bx = B_array(la, lb, lc, ld, p[0], a[0], b[0], q[0], c[0], d[0], gamma1, gamma2, delta);
    std::vector<double> by = B_array(ma, mb, mc, md, p[1], a[1], b[1], q[1], c[1], d[1], gamma1, gamma2, delta);
    std::vector<double> bz = B_array(na, nb, nc, nd, p[2], a[2], b[2], q[2], c[2], d[2], gamma1, gamma2, delta);

    double sum = 0.0;
    for(int i=0; i<=(la+lb+lc+ld); i++) {
        for(int j=0; j<=(ma+mb+mc+md); j++) {
            for(int k=0; k<=(na+nb+nc+nd); k++) {
                sum += bx[i]*by[j]*bz[k]*this->gamma_inc.Fgamma(i+j+k,0.25*rpq2/delta);
            }
        }
    }

    return 2.0 * std::pow(M_PI,2.5)/(gamma1*gamma2*std::sqrt(gamma1+gamma2))*
                 std::exp(-alphaa*alphab*rab2/gamma1)*
                 std::exp(-alphac*alphad*rcd2/gamma2)*sum;
}

/**
 * @brief Performs nuclear integral evaluation
 *
 * This function uses function-level caching of the Fgamma function
 *
 * @param Vec3 a            Center of the Gaussian orbital of the first GTO
 * @param unsigned int l1   Power of x component of the polynomial of the first GTO
 * @param unsigned int m1   Power of y component of the polynomial of the first GTO
 * @param unsigned int n1   Power of z component of the polynomial of the first GTO
 * @param double alpha1     Gaussian exponent of the first GTO
 * @param Vec3 b            Center of the Gaussian orbital of the second GTO
 * @param unsigned int l2   Power of x component of the polynomial of the second GTO
 * @param unsigned int m2   Power of y component of the polynomial of the second GTO
 * @param unsigned int n2   Power of z component of the polynomial of the second GTO
 * @param double alpha2     Gaussian exponent of the second GTO
 * @param Vec3 c
 *
 * @return double value of the nuclear integral
 */
double Integrator::repulsion_fgamma_cached(const Vec3 &a, const int la, const int ma, const int na, const double alphaa,
                                           const Vec3 &b, const int lb, const int mb, const int nb, const double alphab,
                                           const Vec3 &c, const int lc, const int mc, const int nc, const double alphac,
                                           const Vec3 &d, const int ld, const int md, const int nd, const double alphad) const {

    double rab2 = (a-b).norm2();
    double rcd2 = (c-d).norm2();

    Vec3 p = gaussian_product_center(alphaa, a, alphab, b);
    Vec3 q = gaussian_product_center(alphac, c, alphad, d);

    double rpq2 = (p-q).norm2();

    double gamma1 = alphaa + alphab;
    double gamma2 = alphac + alphad;
    double delta = 0.25 * (1.0 / gamma1 + 1.0 / gamma2);

    std::vector<double> bx = B_array(la, lb, lc, ld, p[0], a[0], b[0], q[0], c[0], d[0], gamma1, gamma2, delta);
    std::vector<double> by = B_array(ma, mb, mc, md, p[1], a[1], b[1], q[1], c[1], d[1], gamma1, gamma2, delta);
    std::vector<double> bz = B_array(na, nb, nc, nd, p[2], a[2], b[2], q[2], c[2], d[2], gamma1, gamma2, delta);

    // pre-calculate Fgamma values and cache them inside a vector
    const unsigned int imax = la+lb+lc+ld;
    const unsigned int jmax = ma+mb+mc+md;
    const unsigned int kmax = na+nb+nc+nd;
    std::vector<double> Fgamma_lookup(imax + jmax + kmax + 1, 0.0);
    for(unsigned int i=0; i<Fgamma_lookup.size(); i++) {
        Fgamma_lookup[i] = this->gamma_inc.Fgamma(i,0.25*rpq2/delta);
    }

    double sum = 0.0;
    for(int i=0; i<=(la+lb+lc+ld); i++) {
        for(int j=0; j<=(ma+mb+mc+md); j++) {
            for(int k=0; k<=(na+nb+nc+nd); k++) {
                sum += bx[i]*by[j]*bz[k]*Fgamma_lookup[i+j+k];
            }
        }
    }

    return 2.0 * std::pow(M_PI,2.5)/(gamma1*gamma2*std::sqrt(gamma1+gamma2))*
                 std::exp(-alphaa*alphab*rab2/gamma1)*
                 std::exp(-alphac*alphad*rcd2/gamma2)*sum;
}

std::vector<double> Integrator::B_array(const int l1,const int l2,const int l3,const int l4,
        const double p, const double a, const double b, const double q, const double c, const double d,
        const double g1, const double g2, const double delta) const {

    int imax = l1 + l2 + l3 + l4 + 1;
    std::vector<double> arrB(imax,0);

    for(int i1=0; i1<l1+l2+1; i1++) {
        for(int i2=0; i2<l3+l4+1; i2++) {
            for(int r1=0; r1 < i1/2+1; r1++) {
                for(int r2=0; r2 < i2/2+1; r2++) {
                    for(int u=0; u<(i1+i2)/2-r1-r2+1; u++) {
                        int i = i1+i2-2*(r1+r2)-u;
                        arrB[i] += this->B_term(i1,i2,r1,r2,u,l1,l2,l3,l4,
                                                            p,a,b,q,c,d,g1,g2,delta);
                    }
                }
            }
        }
    }
    return arrB;
}

double Integrator::B_term(const int i1, const int i2, const int r1, const int r2, const int u, const int l1, const int l2, const int l3, const int l4,
        const double px, const double ax, const double bx, const double qx, const double cx, const double dx, const double gamma1,
        const double gamma2, const double delta) const {
    return fB(i1,l1,l2,px,ax,bx,r1,gamma1)*
        pow(-1,i2) * fB(i2,l3,l4,qx,cx,dx,r2,gamma2)*
        pow(-1,u)*fact_ratio2(i1+i2-2*(r1+r2),u)*
        pow(qx-px,i1+i2-2*(r1+r2)-2*u)/
        pow(delta,i1+i2-2*(r1+r2)-u);
}

double Integrator::fB(const int i, const int l1, const int l2, const double p, const double a, const double b, const int r, const double g) const {
    return binomial_prefactor(i, l1, l2, p-a, p-b) * BB0(i, r, g);
}

double Integrator::BB0(int i, int r, double g) const {
    return fact_ratio2(i,r) * pow(4*g,r-i);
}

double Integrator::fact_ratio2(unsigned int a, unsigned int b) const {
    return factorial(a) / factorial(b) / factorial(a - 2*b);
}

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
