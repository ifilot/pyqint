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

#define _USE_MATH_DEFINES
#include <cmath>

#include <string>
#include <unordered_map>
#include <cstring>
#include <vector>
#include <array>
#include <memory>
#include <algorithm>

#include "fgamma.h"
#include "cgf.h"
#include "factorials.h"
#include "mathfuncs.h"
#include "eri/hellsing/hellsing_cache.h"
#include "eri/orchestrator.h"

#ifdef _OPENMP
#include <omp.h>
#endif

class Integrator {
private:
    std::string compile_date;
    std::string compile_time;
    std::string openmp_version;
    std::string compiler_version;
    std::string compiler_type;

    BoysFunction boys_function;
    EriEngine eri_engine;
    mutable integrals::eri::hellsing::HellsingCacheTable1D hellsing_cache;

public:
    /**
     * @brief Construct an Integrator instance.
     *
     * @param lmax        Maximum angular momentum to cache initially
     * @param nu_max      Maximum Boys function order to tabulate
     * @param eri_engine  Primitive electron repulsion engine
     */
    Integrator(int lmax=4, int nu_max=12, EriEngine eri_engine=EriEngine::Hellsing);

    /**
     * @brief Return the number of OpenMP threads available.
     *
     * @return number of threads, or -1 when OpenMP is unavailable
     */
    inline int get_num_threads() const {
        #ifdef _OPENMP
            int numthreads = 1;
            #pragma omp parallel
            {
                #pragma omp single
                numthreads = omp_get_num_threads();
            }
            return numthreads;
        #else
            return -1;
        #endif
    }

    /**
     * @brief      Gets the compiler version.
     *
     * @return     The compiler version.
     */
    inline const char* get_compiler_version() const {
        return this->compiler_version.c_str();
    }

    /**
     * @brief      Gets the compile time.
     *
     * @return     The compile time.
     */
    inline const char* get_compile_time() const {
        return this->compile_time.c_str();
    }

    /**
     * @brief      Gets the compile date.
     *
     * @return     The compile date.
     */
    inline const char* get_compile_date() const {
        return this->compile_date.c_str();
    }

    /**
     * @brief      Gets the openmp version.
     *
     * @return     The openmp version.
     */
    inline const char* get_openmp_version() const {
        return this->openmp_version.c_str();
    }

    /**
     * @brief      Gets the compiler type.
     *
     * @return     The compiler type.
     */
    inline const char* get_compiler_type() const {
        return this->compiler_type.c_str();
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
     * @return flat vector of overlap, kinetic, nuclear, and dense TEI values
     */
    std::vector<double> evaluate_cgfs(const std::vector<CGF>& cgfs,
                                      const std::vector<int>& charges,
                                      const std::vector<double>& px,
                                      const std::vector<double>& py,
                                      const std::vector<double>& pz,
                                      size_t planner_budget) const;

    /** @brief Evaluate overlap, kinetic, and nuclear attraction matrices in sequence. */
    std::vector<double> evaluate_one_electron_integrals(const std::vector<CGF>& cgfs,
                                                        const std::vector<int>& charges,
                                                        const std::vector<double>& px,
                                                        const std::vector<double>& py,
                                                        const std::vector<double>& pz) const;

    /** @brief Evaluate Cartesian position integrals relative to an origin. */
    std::vector<double> evaluate_dipoles(const std::vector<CGF>& cgfs, double x, double y, double z) const;

    /**
     * @brief Evaluate all two-electron integrals within the basis set.
     *
     * @param cgfs  Contracted Gaussian functions
     * @param planner_budget  Bytes available for ERI quartet planning
     *
     * @return dense vector of two-electron integrals
     */
    std::vector<double> evaluate_tei(const std::vector<CGF>& cgfs, size_t planner_budget) const;

    /**
     * @brief Evaluate the ERI tensor and record separate setup and kernel timings.
     */
    std::vector<double> evaluate_tei(const std::vector<CGF>& cgfs, integrals::eri::EriBuildProfile& profile,
                                     size_t planner_budget) const;

    /**
     * @brief Evaluate each unique ERI in teindex order.
     * @param cgfs    Contracted Gaussian functions
     * @param values  Zero-initialized array with one element per unique ERI
     * @param profile Optional setup and evaluation timings
     * @param planner_budget  Bytes available for ERI quartet planning
     */
    void evaluate_packed_tei(const std::vector<CGF>& cgfs, double* values,
                             integrals::eri::EriBuildProfile* profile, size_t planner_budget) const;

    /**
     * @brief Contract packed ERIs into dense Coulomb and exchange matrices.
     * @param n          Number of basis functions
     * @param values     Unique ERIs in teindex order
     * @param density_a  Restricted or alpha density matrix in row-major order
     * @param density_b  Beta density matrix, or nullptr for restricted contraction
     * @param coulomb    Coulomb matrix from the total density
     * @param exchange_a Restricted or alpha exchange matrix
     * @param exchange_b Beta exchange matrix, or nullptr for restricted contraction
     */
    void contract_packed_tei(size_t n, const double* values, const double* density_a,
                             const double* density_b, double* coulomb, double* exchange_a,
                             double* exchange_b) const;

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
    std::vector<double> evaluate_geometric_derivatives(const std::vector<CGF>& cgfs,
                                                       const std::vector<int>& charges,
                                                       const std::vector<double>& px,
                                                       const std::vector<double>& py,
                                                       const std::vector<double>& pz) const;

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
    double overlap(const CGF& cgf1, const CGF& cgf2) const;

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
    double overlap_gto(const GTO& gto1, const GTO& gto2) const;

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
    double overlap_deriv(const CGF& cgf1,
                         const CGF& cgf2,
                         const Vec3& nucleus,
                         unsigned int coord) const;

    /**
     * @brief Calculate the geometric derivative of the overlap integral for two GTOs.
     *
     * @param gto1   Gaussian Type Orbital 1
     * @param gto2   Gaussian Type Orbital 2
     * @param coord  Cartesian direction (0=x,1=y,2=z)
     *
     * Calculates the value of d/dx < gto1 | gto2 >
     *
     * @return value of the overlap derivative
     */
    double overlap_deriv_gto(const GTO& gto1,
                             const GTO& gto2,
                             unsigned int coord) const;

/**************************************************************************
 * DIPOLE INTEGRALS
 **************************************************************************/

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
    double dipole(const CGF& cgf1,
                  const CGF& cgf2,
                  unsigned int cc,
                  double ref = 0.0) const;
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
    double dipole_gto(const GTO& gto1,
                      const GTO& gto2,
                      unsigned int cc,
                      double ref = 0.0) const;

    // expanded interface for Cython
    /**
     * @brief Calculate the geometric derivative of the overlap integral.
     *
     * @param cgf1   Contracted Gaussian Function 1
     * @param cgf2   Contracted Gaussian Function 2
     * @param cx     nucleus x coordinate
     * @param cy     nucleus y coordinate
     * @param cz     nucleus z coordinate
     * @param coord  Cartesian direction (0=x,1=y,2=z)
     *
     * @return value of the overlap derivative
     */
    inline double overlap_deriv(const CGF& cgf1, const CGF& cgf2, double cx, double cy, double cz, unsigned int coord) const {
        return this->overlap_deriv(cgf1, cgf2, Vec3(cx, cy, cz), coord);
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
    double kinetic(const CGF& cgf1, const CGF& cgf2) const;

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
    double kinetic_gto(const GTO& gto1, const GTO& gto2) const;

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
    double kinetic_deriv(const CGF& cgf1,
                         const CGF& cgf2,
                         const Vec3& nucleus,
                         unsigned int coord) const;


    /**
     * @brief Calculate kinetic derivative (Cython interface).
     *
     * Expanded interface for Cython
     *
     * @param cgf1   Contracted Gaussian Function 1
     * @param cgf2   Contracted Gaussian Function 2
     * @param cx     nucleus position x
     * @param cy     nucleus position y
     * @param cz     nucleus position z
     * @param coord  Cartesian direction (0=x,1=y,2=z)
     *
     * @return value of the kinetic integral derivative
     */
    inline double kinetic_deriv(const CGF& cgf1,
                                const CGF& cgf2,
                                double cx, double cy, double cz,
                                unsigned int coord) const {
        return this->kinetic_deriv(cgf1, cgf2, Vec3(cx, cy, cz), coord);
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
    double kinetic_deriv_gto(const GTO& gto1,
                             const GTO& gto2,
                             unsigned int coord) const;

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
    double nuclear(const CGF &cgf1,
                   const CGF &cgf2,
                   const Vec3& nucleus,
                   unsigned int charge) const;

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
    double nuclear_gto(const GTO &gto1,
                       const GTO &gto2,
                       const Vec3& nucleus) const;


    /**
     * @brief Calculate nuclear attraction integral of two CGFs.
     *
     * @param cgf1    Contracted Gaussian Function 1
     * @param cgf2    Contracted Gaussian Function 2
     * @param cx      nuclear coordinate x
     * @param cy      nuclear coordinate y
     * @param cz      nuclear coordinate z
     * @param charge  charge of the nucleus
     *
     * Calculates the value of < cgf1 | V | cgf2 >
     *
     * Expanded notation for Cython interface
     *
     * @return     double value of the nuclear integral
     */
    inline double nuclear(const CGF &cgf1,
                          const CGF &cgf2,
                          double cx, double cy, double cz,
                          unsigned int charge) const {
        return this->nuclear(cgf1, cgf2, Vec3(cx, cy, cz), charge);
    }

    /**
     * @brief Calculates derivative of the nuclear attraction integral of two CGF
     *
     * @param cgf1      Contracted Gaussian Function
     * @param cgf2      Contracted Gaussian Function
     * @param nucleus   Position of the nucleus generating the potential
     * @param charge    Charge of the nucleus in a.u.
     * @param nucderiv  Center with respect to which the derivative is taken
     * @param coord     Cartesian direction of the derivative (0=x,1=y,2=z)
     *
     * Calculates the value of
     *
     *     d/dR ⟨ cgf1 | V_nuc | cgf2 ⟩
     *
     * where R is the selected Cartesian coordinate of the given center.
     *
     * In contrast to the other integrals, this routine evaluates the
     * nuclear attraction **derivative** directly and returns the
     * corresponding energy gradient component.
     *
     * @return double value of the nuclear attraction gradient
     */
    double nuclear_deriv(const CGF &cgf1, const CGF &cgf2, const Vec3& nucleus, 
                         unsigned int charge, const Vec3& nucderiv, unsigned int coord) const;

    // expanded notation for Cython interface
    /**
     * @brief Calculate nuclear attraction integral derivatives (Cython interface).
     *
     * @param cgf1   Contracted Gaussian Function
     * @param cgf2   Contracted Gaussian Function
     * @param cx     nucleus x coordinate
     * @param cy     nucleus y coordinate
     * @param cz     nucleus z coordinate
     * @param charge Charge of the nucleus in a.u.
     * @param dx     derivative center x coordinate
     * @param dy     derivative center y coordinate
     * @param dz     derivative center z coordinate
     * @param coord  Cartesian direction (0=x,1=y,2=z)
     *
     * @return value of the nuclear attraction gradient
     */
    inline double nuclear_deriv(const CGF &cgf1, const CGF &cgf2, double cx, double cy, double cz, unsigned int charge,
                                double dx, double dy, double dz, unsigned int coord) const {
        return this->nuclear_deriv(cgf1, cgf2, Vec3(cx, cy, cz), charge, Vec3(dx, dy, dz), coord);
    }

    /**
     * @brief Calculate nuclear attraction integral of two GTOs (Cython interface).
     *
     * @param gto1  Gaussian Type Orbital
     * @param gto2  Gaussian Type Orbital
     * @param cx    nucleus x coordinate
     * @param cy    nucleus y coordinate
     * @param cz    nucleus z coordinate
     *
     * Calculates the value of < gto1 | V | gto2 >
     *
     * @return double value of the nuclear integral
     */
    inline double nuclear_gto(const GTO &gto1, const GTO &gto2, double cx, double cy, double cz) const {
        return this->nuclear_gto(gto1, gto2, Vec3(cx, cy, cz));
    }

    /**
     * @brief Calculates derivative of the nuclear attraction integral of two GTO
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
     * This routine evaluates the derivative of the nuclear attraction
     * integral and returns the corresponding energy gradient component.
     *
     * @return double value of the nuclear attraction gradient
     */
    double nuclear_deriv(const GTO& gto1, const GTO& gto2, const Vec3& nucleus,
                         unsigned int coord) const;

/**************************************************************************
 * TWO-ELECTRON INTEGRALS
 **************************************************************************/

    /**
     * @brief Calculates two-electron repulsion integral of four CGF
     *
     * @param cgf1  Contracted Gaussian Function
     * @param cgf2  Contracted Gaussian Function
     * @param cgf3  Contracted Gaussian Function
     * @param cgf4  Contracted Gaussian Function
     *
     * Calculates the value of < cgf1 | cgf2 | cgf3 | cgf4 >
     *
     * @return double value of the repulsion integral
     */
    double repulsion(const CGF &cgf1, const CGF &cgf2, const CGF &cgf3, const CGF &cgf4) const;

    /**
     * @brief Calculates two-electron repulsion integral of four CGF
     *
     * @param gto1  Gaussian Type Orbital
     * @param gto2  Gaussian Type Orbital
     * @param gto3  Gaussian Type Orbital
     * @param gto4  Gaussian Type Orbital
     *
     * Calculates the value of < gto1 | gto2 | gto3 | gto4 >
     *
     * @return double value of the repulsion integral
     */
    double repulsion(const GTO &gto1, const GTO &gto2, const GTO &gto3, const GTO &gto4) const;

    /**
     * @brief Calculates derivative of the two-electron repulsion integral of four CGF
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
     * @return double value of the repulsion integral
     */
    double repulsion_deriv(const CGF &cgf1,const CGF &cgf2,const CGF &cgf3,const CGF &cgf4,
        const Vec3& nucleus, unsigned int coord) const;

    /**
     * @brief Calculate derivative of the two-electron repulsion integral (Cython interface).
     *
     * @param cgf1   Contracted Gaussian Function
     * @param cgf2   Contracted Gaussian Function
     * @param cgf3   Contracted Gaussian Function
     * @param cgf4   Contracted Gaussian Function
     * @param cx     nucleus x coordinate
     * @param cy     nucleus y coordinate
     * @param cz     nucleus z coordinate
     * @param coord  Cartesian direction (0=x,1=y,2=z)
     *
     * @return value of the repulsion integral derivative
     */
    inline double repulsion_deriv(const CGF &cgf1,const CGF &cgf2,const CGF &cgf3,const CGF &cgf4,
        double cx, double cy, double cz, unsigned int coord) const {
        return this->repulsion_deriv(cgf1, cgf2, cgf3, cgf4, Vec3(cx, cy, cz), coord);
    }

    /**
     * @brief Calculates overlap integral of two GTO
     *
     * @param gto1   Gaussian Type Orbital
     * @param gto2   Gaussian Type Orbital
     * @param gto3   Gaussian Type Orbital
     * @param gto4   Gaussian Type Orbital
     * @param coord  Cartesian direction (0=x,1=y,2=z)
     *
     * Calculates the value of < d/dx gto1 | gto2 | gto3 | gto4 >
     *
     * @return double value of the overlap integral
     */
    double repulsion_deriv(const GTO& gto1, const GTO& gto2, const GTO &gto3, const GTO &gto4, unsigned int coord) const;

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
    size_t teindex(size_t i, size_t j, size_t k, size_t l) const;

    /**
     * @brief Ensure the Hellsing cache supports the basis set.
     *
     * @param cgf1  Contracted Gaussian Function
     * @param cgf2  Contracted Gaussian Function
     * @param cgf3  Contracted Gaussian Function
     * @param cgf4  Contracted Gaussian Function
     */
    void ensure_hellsing_cache(const CGF &cgf1, const CGF &cgf2, const CGF &cgf3, const CGF &cgf4);

private:
    /**
     * @brief Describe how a contracted function responds to moving one nucleus.
     *
     * Uniform motion uses one flag for the full contraction. A mixed contraction
     * uses per-primitive flags and has a geometry-dependent normalization.
     */
    struct CgfDerivative {
        const CGF& cgf;
        bool first_moves;
        std::vector<bool> primitive_moves;
        std::array<double, 3> log_norm_deriv;

        bool moves(unsigned int i) const {
            return primitive_moves.empty() ? first_moves : primitive_moves[i];
        }

        bool all_move() const { return primitive_moves.empty() && first_moves; }
        bool none_move() const { return primitive_moves.empty() && !first_moves; }
    };

    /**
     * @brief Identify the moving primitives and differentiate their contraction norm.
     *
     * @param cgf      Contracted Gaussian function
     * @param nucleus  Position of the nucleus being displaced
     * @param shared_center  Verified common center when available
     *
     * @return Primitive motion and the derivatives of the logarithmic normalization
     */
    CgfDerivative prepare_cgf_derivative(const CGF& cgf, const Vec3& nucleus,
                                         const std::optional<Vec3>& shared_center = std::nullopt) const;

    /**
     * @brief Match a Gaussian or operator center to a nuclear position.
     */
    static bool same_center(const Vec3& center, const Vec3& nucleus) {
        return (center - nucleus).norm2() < 0.0001;
    }

    using ContractedIntegral = double (Integrator::*)(const CGF&, const CGF&) const;
    using PrimitiveDerivative = double (Integrator::*)(const GTO&, const GTO&, unsigned int) const;

    /**
     * @brief Differentiate a contracted overlap or kinetic integral.
     */
    double two_center_deriv(const CgfDerivative& bra, const CgfDerivative& ket, unsigned int coord,
                            ContractedIntegral integral, PrimitiveDerivative derivative) const;

    /**
     * @brief Differentiate nuclear attraction using prepared contraction derivatives.
     */
    double nuclear_deriv(const CgfDerivative& bra, const CgfDerivative& ket, const Vec3& nucleus,
                         unsigned int charge, bool nucleus_moves, unsigned int coord) const;

    /**
     * @brief Differentiate electron repulsion using prepared contraction derivatives.
     */
    double repulsion_deriv(const CgfDerivative& a, const CgfDerivative& b,
                           const CgfDerivative& c, const CgfDerivative& d, unsigned int coord) const;

    /**
     * @brief Evaluate a contracted quartet with its selected kernel and optional prepared pairs.
     */
    double repulsion(const CGF& a, const CGF& b, const CGF& c, const CGF& d, EriEngine engine,
                     const integrals::eri::EriOrchestrator::Pair* bra = nullptr,
                     const integrals::eri::EriOrchestrator::Pair* ket = nullptr) const;

    /**
     * @brief Evaluate a primitive quartet using the selected ERI kernel.
     */
    double repulsion(const GTO& a, const GTO& b, const GTO& c, const GTO& d, EriEngine engine) const;

    /**
     * @brief Performs overlap integral evaluation
     *
     * @param alpha1  Gaussian exponent of the first GTO
     * @param l1      Power of x component of the polynomial of the first GTO
     * @param m1      Power of y component of the polynomial of the first GTO
     * @param n1      Power of z component of the polynomial of the first GTO
     * @param a       Center of the Gaussian orbital of the first GTO
     * @param alpha2  Gaussian exponent of the second GTO
     * @param l2      Power of x component of the polynomial of the second GTO
     * @param m2      Power of y component of the polynomial of the second GTO
     * @param n2      Power of z component of the polynomial of the second GTO
     * @param b       Center of the Gaussian orbital of the second GTO
     *
     * @return double value of the overlap integral
     */
    double overlap(double alpha1, unsigned int l1, unsigned int m1, unsigned int n1, const Vec3 &a,
                   double alpha2, unsigned int l2, unsigned int m2, unsigned int n2, const Vec3 &b) const;

    /**
     * @brief Performs dipole integral evaluation
     *
     * @param alpha1  Gaussian exponent of the first GTO
     * @param l1      Power of x component of the polynomial of the first GTO
     * @param m1      Power of y component of the polynomial of the first GTO
     * @param n1      Power of z component of the polynomial of the first GTO
     * @param a       Center of the Gaussian orbital of the first GTO
     * @param alpha2  Gaussian exponent of the second GTO
     * @param l2      Power of x component of the polynomial of the second GTO
     * @param m2      Power of y component of the polynomial of the second GTO
     * @param n2      Power of z component of the polynomial of the second GTO
     * @param b       Center of the Gaussian orbital of the second GTO
     * @param cc      Cartesian direction (0=x,1=y,2=z)
     * @param cref    Reference position
     *
     * @return double value of the overlap integral
     */
    double dipole(double alpha1, unsigned int l1, unsigned int m1, unsigned int n1, const Vec3 &a,
                  double alpha2, unsigned int l2, unsigned int m2, unsigned int n2, const Vec3 &b,
                  unsigned int cc, double cref = 0.0) const;

    /**
     * @brief Performs nuclear integral evaluation
     *
     * @param a       Center of the Gaussian orbital of the first GTO
     * @param l1      Power of x component of the polynomial of the first GTO
     * @param m1      Power of y component of the polynomial of the first GTO
     * @param n1      Power of z component of the polynomial of the first GTO
     * @param alpha1  Gaussian exponent of the first GTO
     * @param b       Center of the Gaussian orbital of the second GTO
     * @param l2      Power of x component of the polynomial of the second GTO
     * @param m2      Power of y component of the polynomial of the second GTO
     * @param n2      Power of z component of the polynomial of the second GTO
     * @param alpha2  Gaussian exponent of the second GTO
     * @param c       Nuclear position
     *
     * @return double value of the nuclear integral
     */
    double nuclear(const Vec3& a,
                   int l1, int m1, int n1,
                   double alpha1,
                   const Vec3& b,
                   int l2, int m2, int n2,
                   double alpha2,
                   const Vec3& c) const;

    /**
     * @brief Performs nuclear integral evaluation
     *
     * @param a       Center of the Gaussian orbital of the first GTO
     * @param l1      Power of x component of the polynomial of the first GTO
     * @param m1      Power of y component of the polynomial of the first GTO
     * @param n1      Power of z component of the polynomial of the first GTO
     * @param alpha1  Gaussian exponent of the first GTO
     * @param b       Center of the Gaussian orbital of the second GTO
     * @param l2      Power of x component of the polynomial of the second GTO
     * @param m2      Power of y component of the polynomial of the second GTO
     * @param n2      Power of z component of the polynomial of the second GTO
     * @param alpha2  Gaussian exponent of the second GTO
     * @param c       Nuclear position
     * @param coord   Cartesian direction to derive nuclear coordinate towards
     *
     * @return double value of the nuclear integral derived towards nuclear coordinate
     */
    double nuclear_deriv_op(const Vec3& a, int l1, int m1, int n1, double alpha1,
                            const Vec3& b, int l2, int m2, int n2,
                            double alpha2, const Vec3& c, unsigned int coord) const;

    /**
     * @brief Calculates one dimensional overlap integral
     *
     * @param l1     Power of 'x' component of the polynomial of the first GTO
     * @param l2     Power of 'x' component of the polynomial of the second GTO
     * @param x1     'x' component of the position of the first GTO
     * @param x2     'x' component of the position of the second GTO
     * @param gamma  Sum of the two Gaussian exponents
     *
     * @return double value of the one dimensional overlap integral
     */
    double overlap_1D(int l1, int l2, double x1, double x2, double gamma) const;

    /************************
     *
     * AUXILIARY FUNCTIONS
     *
     ************************/

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
    Vec3 gaussian_product_center(double alpha1, const Vec3 &a,
                                 double alpha2, const Vec3 &b) const;

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
    double binomial_prefactor(int s, int ia, int ib, double xpa, double xpb) const;

    /**
     * @brief Compute the binomial coefficient.
     *
     * @param a  Upper index
     * @param b  Lower index
     *
     * @return binomial coefficient
     */
    double binomial(int a, int b) const;

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
    std::vector<double> A_array(const int l1, const int l2,
                                const double pa, const double pb,
                                const double cp, const double g) const;

    /**
     * @brief Build the derivative of the Hermite A array.
     *
     * @param l1  Angular momentum component for center A
     * @param l2  Angular momentum component for center B
     * @param pa  P-A distance component
     * @param pb  P-B distance component
     * @param cp  C-P distance component
     * @param g   Gaussian exponent sum
     *
     * @return Hermite A array derivative
     */
    std::vector<double> A_array_deriv(const int l1, const int l2,
                                      const double pa, const double pb,
                                      const double cp, const double g) const;

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
    double A_term(const int i, const int r, const int u,
                  const int l1, const int l2,
                  const double pax, const double pbx,
                  const double cpx, const double gamma) const;

    /**
     * @brief Evaluate the Boys function in gamma form.
     *
     * @param m  Boys order
     * @param x  argument
     *
     * @return Boys function value
     */
    double gamma(const double m, double x) const;

    /**
     * @brief Compute all two-electron integrals and expand them into a dense tensor.
     *
     * @param cgfs      Contracted Gaussian functions
     * @param tetensor  Output tensor of two-electron integrals
     */
    void calculate_two_electron_integrals(const std::vector<CGF>& cgfs, std::vector<double>& tetensor,
                                         integrals::eri::EriBuildProfile* profile, size_t planner_budget) const;
};
