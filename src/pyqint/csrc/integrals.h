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
#include "hellsing_cache.h"

#ifdef _OPENMP
#include <omp.h>
#endif

struct HellsingBTerm {
    double c;  // coefficient
    int mu;    // mu
    int u;     // u
};

class Integrator {
private:
    std::string compile_date;
    std::string compile_time;
    std::string openmp_version;
    std::string compiler_version;
    std::string compiler_type;

    BoysFunction boys_function;
    mutable HellsingCacheTable1D hellsing_cache;

public:
    /**
     * @brief Construct an Integrator instance.
     *
     * @param lmax    Maximum angular momentum
     * @param nu_max  Maximum Boys function order
     */
    Integrator(int lmax=4, int nu_max=12);

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
     *
     * @return packed vector of overlap, kinetic, nuclear, and TEI values
     */
    std::vector<double> evaluate_cgfs(const std::vector<CGF>& cgfs,
                                      const std::vector<int>& charges,
                                      const std::vector<double>& px,
                                      const std::vector<double>& py,
                                      const std::vector<double>& pz) const;

    /**
     * @brief Evaluate all two-electron integrals within the basis set.
     *
     * @param cgfs  Contracted Gaussian functions
     *
     * @return packed vector of two-electron integrals
     */
    std::vector<double> evaluate_tei(const std::vector<CGF>& cgfs) const;

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
     * @brief Performs electron repulsion integral evaluation.
     *
     * @param a       Center of the first GTO
     * @param la      Power of x component of the polynomial of the first GTO
     * @param ma      Power of y component of the polynomial of the first GTO
     * @param na      Power of z component of the polynomial of the first GTO
     * @param alphaa  Gaussian exponent of the first GTO
     * @param b       Center of the second GTO
     * @param lb      Power of x component of the polynomial of the second GTO
     * @param mb      Power of y component of the polynomial of the second GTO
     * @param nb      Power of z component of the polynomial of the second GTO
     * @param alphab  Gaussian exponent of the second GTO
     * @param c       Center of the third GTO
     * @param lc      Power of x component of the polynomial of the third GTO
     * @param mc      Power of y component of the polynomial of the third GTO
     * @param nc      Power of z component of the polynomial of the third GTO
     * @param alphac  Gaussian exponent of the third GTO
     * @param d       Center of the fourth GTO
     * @param ld      Power of x component of the polynomial of the fourth GTO
     * @param md      Power of y component of the polynomial of the fourth GTO
     * @param nd      Power of z component of the polynomial of the fourth GTO
     * @param alphad  Gaussian exponent of the fourth GTO
     *
     * @return Value of the electron repulsion integral
     */
    double repulsion(const Vec3 &a, const int la, const int ma, const int na, const double alphaa,
                     const Vec3 &b, const int lb, const int mb, const int nb, const double alphab,
                     const Vec3 &c, const int lc, const int mc, const int nc, const double alphac,
                     const Vec3 &d, const int ld, const int md, const int nd, const double alphad) const;

    /**
     * @brief Performs electron repulsion integral (ERI) evaluation with cached Boys function
     *
     * @param a        Center of the first Gaussian-type orbital (GTO)
     * @param la       Power of x component of the polynomial of the first GTO
     * @param ma       Power of y component of the polynomial of the first GTO
     * @param na       Power of z component of the polynomial of the first GTO
     * @param alphaa  Gaussian exponent of the first GTO
     *
     * @param b        Center of the second Gaussian-type orbital (GTO)
     * @param lb       Power of x component of the polynomial of the second GTO
     * @param mb       Power of y component of the polynomial of the second GTO
     * @param nb       Power of z component of the polynomial of the second GTO
     * @param alphab  Gaussian exponent of the second GTO
     *
     * @param c        Center of the third Gaussian-type orbital (GTO)
     * @param lc       Power of x component of the polynomial of the third GTO
     * @param mc       Power of y component of the polynomial of the third GTO
     * @param nc       Power of z component of the polynomial of the third GTO
     * @param alphac  Gaussian exponent of the third GTO
     *
     * @param d        Center of the fourth Gaussian-type orbital (GTO)
     * @param ld       Power of x component of the polynomial of the fourth GTO
     * @param md       Power of y component of the polynomial of the fourth GTO
     * @param nd       Power of z component of the polynomial of the fourth GTO
     * @param alphad  Gaussian exponent of the fourth GTO
     *
     * @return Value of the electron repulsion integral
     */
    double repulsion_boys_cached(const Vec3 &a, const int la, const int ma, const int na, const double alphaa,
                                 const Vec3 &b, const int lb, const int mb, const int nb, const double alphab,
                                 const Vec3 &c, const int lc, const int mc, const int nc, const double alphac,
                                 const Vec3 &d, const int ld, const int md, const int nd, const double alphad) const;

    /**
     * @brief Performs electron repulsion integral (ERI) evaluation with in-situ Hellsing engine
     *
     * @param a        Center of the first Gaussian-type orbital (GTO)
     * @param la       Power of x component of the polynomial of the first GTO
     * @param ma       Power of y component of the polynomial of the first GTO
     * @param na       Power of z component of the polynomial of the first GTO
     * @param alphaa  Gaussian exponent of the first GTO
     *
     * @param b        Center of the second Gaussian-type orbital (GTO)
     * @param lb       Power of x component of the polynomial of the second GTO
     * @param mb       Power of y component of the polynomial of the second GTO
     * @param nb       Power of z component of the polynomial of the second GTO
     * @param alphab  Gaussian exponent of the second GTO
     *
     * @param c        Center of the third Gaussian-type orbital (GTO)
     * @param lc       Power of x component of the polynomial of the third GTO
     * @param mc       Power of y component of the polynomial of the third GTO
     * @param nc       Power of z component of the polynomial of the third GTO
     * @param alphac  Gaussian exponent of the third GTO
     *
     * @param d        Center of the fourth Gaussian-type orbital (GTO)
     * @param ld       Power of x component of the polynomial of the fourth GTO
     * @param md       Power of y component of the polynomial of the fourth GTO
     * @param nd       Power of z component of the polynomial of the fourth GTO
     * @param alphad  Gaussian exponent of the fourth GTO
     *
     * @return Value of the electron repulsion integral
     */
    double repulsion_hellsing (const Vec3 &a, const int la, const int ma, const int na, const double alphaa,
                               const Vec3 &b, const int lb, const int mb, const int nb, const double alphab,
                               const Vec3 &c, const int lc, const int mc, const int nc, const double alphac,
                               const Vec3 &d, const int ld, const int md, const int nd, const double alphad) const;

    /**
     * @brief Performs electron repulsion integral (ERI) evaluation with cached Hellsing engine
     *
     * @param a        Center of the first Gaussian-type orbital (GTO)
     * @param la       Power of x component of the polynomial of the first GTO
     * @param ma       Power of y component of the polynomial of the first GTO
     * @param na       Power of z component of the polynomial of the first GTO
     * @param alphaa  Gaussian exponent of the first GTO
     *
     * @param b        Center of the second Gaussian-type orbital (GTO)
     * @param lb       Power of x component of the polynomial of the second GTO
     * @param mb       Power of y component of the polynomial of the second GTO
     * @param nb       Power of z component of the polynomial of the second GTO
     * @param alphab  Gaussian exponent of the second GTO
     *
     * @param c        Center of the third Gaussian-type orbital (GTO)
     * @param lc       Power of x component of the polynomial of the third GTO
     * @param mc       Power of y component of the polynomial of the third GTO
     * @param nc       Power of z component of the polynomial of the third GTO
     * @param alphac  Gaussian exponent of the third GTO
     *
     * @param d        Center of the fourth Gaussian-type orbital (GTO)
     * @param ld       Power of x component of the polynomial of the fourth GTO
     * @param md       Power of y component of the polynomial of the fourth GTO
     * @param nd       Power of z component of the polynomial of the fourth GTO
     * @param alphad  Gaussian exponent of the fourth GTO
     *
     * @return Value of the electron repulsion integral
     */
    double repulsion_hellsing_cached (const Vec3 &a, const int la, const int ma, const int na, const double alphaa,
                                      const Vec3 &b, const int lb, const int mb, const int nb, const double alphab,
                                      const Vec3 &c, const int lc, const int mc, const int nc, const double alphac,
                                      const Vec3 &d, const int ld, const int md, const int nd, const double alphad) const;

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
     * @brief Build the B array for the Obara–Saika recurrence.
     *
     * @param l1     Angular momentum component for center A
     * @param l2     Angular momentum component for center B
     * @param l3     Angular momentum component for center C
     * @param l4     Angular momentum component for center D
     * @param p      Gaussian exponent sum for the first pair
     * @param a      Center A coordinate
     * @param b      Center B coordinate
     * @param q      Gaussian exponent sum for the second pair
     * @param c      Center C coordinate
     * @param d      Center D coordinate
     * @param g1     Gaussian exponent sum gamma1
     * @param g2     Gaussian exponent sum gamma2
     * @param delta  delta parameter
     *
     * @return B array values
     */
    std::vector<double> B_array(const int l1,const int l2,const int l3,const int l4,
                                const double p, const double a, const double b, const double q, const double c, const double d,
                                const double g1, const double g2, const double delta) const;

    /**
     * @brief Build the Hellsing B array terms.
     *
     * @param l1  Angular momentum component for center A
     * @param l2  Angular momentum component for center B
     * @param l3  Angular momentum component for center C
     * @param l4  Angular momentum component for center D
     * @param a1  Gaussian exponent of center A
     * @param a2  Gaussian exponent of center B
     * @param a3  Gaussian exponent of center C
     * @param a4  Gaussian exponent of center D
     * @param ax  A coordinate component
     * @param bx  B coordinate component
     * @param cx  C coordinate component
     * @param dx  D coordinate component
     * @param px  P coordinate component
     * @param qx  Q coordinate component
     * @param g1  Gaussian exponent sum gamma1
     * @param g2  Gaussian exponent sum gamma2
     *
     * @return Hellsing B-term array
     */
    std::vector<HellsingBTerm> B_array_hellsing(
        int l1, int l2, int l3, int l4,
        double a1, double a2, double a3, double a4,
        double ax, double bx, double cx, double dx,
        double px, double qx,
        double g1, double g2) const;

    /**
     * @brief Evaluate a single B-term for the Obara–Saika recurrence.
     *
     * @param i1     Hermite index for the first pair
     * @param i2     Hermite index for the second pair
     * @param r1     Recursion index for the first pair
     * @param r2     Recursion index for the second pair
     * @param u      Boys index
     * @param l1     Angular momentum component for center A
     * @param l2     Angular momentum component for center B
     * @param l3     Angular momentum component for center C
     * @param l4     Angular momentum component for center D
     * @param px     P coordinate component
     * @param ax     A coordinate component
     * @param bx     B coordinate component
     * @param qx     Q coordinate component
     * @param cx     C coordinate component
     * @param dx     D coordinate component
     * @param gamma1 Gaussian exponent sum gamma1
     * @param gamma2 Gaussian exponent sum gamma2
     * @param delta  delta parameter
     *
     * @return B-term value
     */
    double B_term(const int i1, const int i2, const int r1, const int r2, const int u, const int l1, const int l2, const int l3, const int l4,
    const double px, const double ax, const double bx, const double qx, const double cx, const double dx, const double gamma1,
    const double gamma2, const double delta) const;

    /**
     * @brief Evaluate helper term used in B-array construction.
     *
     * @param i   Hermite index
     * @param l1  Angular momentum component for center A
     * @param l2  Angular momentum component for center B
     * @param p   Gaussian exponent sum
     * @param a   Center A coordinate
     * @param b   Center B coordinate
     * @param r   Recursion index
     * @param q   Gaussian exponent sum for the second pair
     *
     * @return fB term value
     */
    double fB(const int i, const int l1, const int l2, const double p, const double a, const double b, const int r, const double q) const;

    /**
     * @brief Evaluate the BB0 helper term.
     *
     * @param i  Hermite index
     * @param r  Recursion index
     * @param q  Gaussian exponent sum
     *
     * @return BB0 term value
     */
    double BB0(int i, int r, double q) const;

    /**
     * @brief Compute factorial ratio used in recursion relations.
     *
     * @param a  first index
     * @param b  second index
     *
     * @return factorial ratio
     */
    double fact_ratio2(unsigned int a, unsigned int b) const;

    /**
     * @brief Compute and pack all two-electron integrals.
     *
     * @param cgfs      Contracted Gaussian functions
     * @param tetensor  Output tensor of two-electron integrals
     */
    void calculate_two_electron_integrals(const std::vector<CGF>& cgfs, std::vector<double>& tetensor) const;
};
