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

#include <string>
#include <unordered_map>
#include <cstring>
#include <vector>
#include <array>

#include "gamma.h"
#include "cgf.h"
#include "factorials.h"

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

public:
    /**
     * @brief Integrator constructor method
     *
     * @return Integrator class
     */
    Integrator();

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
     * @brief      Evaluate all integrals for cgfs in buffer
     */
    std::vector<double> evaluate_cgfs(const std::vector<CGF>& cgfs,
                                      const std::vector<int>& charges,
                                      const std::vector<double>& px,
                                      const std::vector<double>& py,
                                      const std::vector<double>& pz) const;

    /**
     * @brief      Evaluate all integrals for cgfs in buffer
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
     * @brief      Calculates overlap integral of two CGF
     *
     * @param[in]  cgf1  Contracted Gaussian Functional 1
     * @param[in]  cgf2  Contracted Gaussian Functional 2
     *
     * Calculates the value of < cgf1 | cgf2 >
     *
     * @return     double value of the overlap integral
     */
    double overlap(const CGF& cgf1, const CGF& cgf2) const;

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
    double overlap_gto(const GTO& gto1, const GTO& gto2) const;

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
    double overlap_deriv(const CGF& cgf1,
                         const CGF& cgf2,
                         const Vec3& nucleus,
                         unsigned int coord) const;

    /**
     * @brief      Calculates the geometric derivative of overlap integral of
     *             two GTOs
     *
     * @param[in]  gto1   Gaussian Type Orbital 1
     * @param[in]  gto2   Gaussian Type Orbital 2
     * @param[in]  coord  Cartesian direction
     * @param[in]  nucleus  Nucleus position
     *
     * Calculates the value of d/dx < gto1 | gto2 >
     *
     * @return     double value of the nuclear integral
     */
    double overlap_deriv_gto(const GTO& gto1,
                             const GTO& gto2,
                             unsigned int coord) const;

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
    double dipole(const CGF& cgf1,
                  const CGF& cgf2,
                  unsigned int cc,
                  double ref = 0.0) const;
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
    double dipole_gto(const GTO& gto1,
                      const GTO& gto2,
                      unsigned int cc,
                      double ref = 0.0) const;

    // expanded interface for Cython
    inline double overlap_deriv(const CGF& cgf1, const CGF& cgf2, double cx, double cy, double cz, unsigned int coord) const {
        return this->overlap_deriv(cgf1, cgf2, Vec3(cx, cy, cz), coord);
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
    double kinetic(const CGF& cgf1, const CGF& cgf2) const;

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
    double kinetic_gto(const GTO& gto1, const GTO& gto2) const;

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
    double kinetic_deriv(const CGF& cgf1,
                         const CGF& cgf2,
                         const Vec3& nucleus,
                         unsigned int coord) const;


    /**
     * @brief      Calculate kinetic derivative
     *
     * Expanded interface for Cython
     *
     * @param[in]  cgf1   The cgf 1
     * @param[in]  cgf2   The cgf 2
     * @param[in]  cx     nucleus position x
     * @param[in]  cy     nucleus position y
     * @param[in]  cz     nucleus position z
     * @param[in]  coord  The coordinate
     *
     * @return     double value of the kinetic integral
     */
    inline double kinetic_deriv(const CGF& cgf1,
                                const CGF& cgf2,
                                double cx, double cy, double cz,
                                unsigned int coord) const {
        return this->kinetic_deriv(cgf1, cgf2, Vec3(cx, cy, cz), coord);
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
    double kinetic_deriv_gto(const GTO& gto1,
                             const GTO& gto2,
                             unsigned int coord) const;

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
    double nuclear(const CGF &cgf1,
                   const CGF &cgf2,
                   const Vec3& nucleus,
                   unsigned int charge) const;

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
    double nuclear_gto(const GTO &gto1,
                       const GTO &gto2,
                       const Vec3& nucleus) const;


    /**
     * @brief      Calculates nuclear integral of two CGF
     *
     * @param[in]  cgf1    Contracted Gaussian Functional 1
     * @param[in]  cgf2    Contracted Gaussian Functional 2
     * @param[in]  cx      nuclear coordinate x
     * @param[in]  cy      nuclear coordinate y
     * @param[in]  cz      nuclear coordinate z
     * @param[in]  charge  charge of the nucleus
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
     * @brief Calculates nuclear integral of two CGF
     *
     * @param const CGF& cgf1       Contracted Gaussian Function
     * @param const CGF& cgf2       Contracted Gaussian Function
     * @param const Vec3 nucleus    Position of the nucleus
     * @param unsigned int charge   charge of the nucleus in a.u.
     *
     * Calculates the value of < cgf1 | V | cgf2 >
     *
     * @return double value of the nuclear integral
     */
    double nuclear_deriv(const CGF &cgf1, const CGF &cgf2, const Vec3& nucleus, unsigned int charge,
                         const Vec3& nucderiv, unsigned int coord) const;

    // expanded notation for Cython interface
    inline double nuclear_deriv(const CGF &cgf1, const CGF &cgf2, double cx, double cy, double cz, unsigned int charge,
                                double dx, double dy, double dz, unsigned int coord) const {
        return this->nuclear_deriv(cgf1, cgf2, Vec3(cx, cy, cz), charge, Vec3(dx, dy, dz), coord);
    }

    /**
     * @brief Calculates nuclear integral of two CGF
     *
     * @param const GTO& gto1       Contracted Gaussian Function
     * @param const GTO& gto2       Contracted Gaussian Function
     * @param unsigned int charge   charge of the nucleus in a.u.
     *
     * Calculates the value of < gto1 | V | gto2 >
     *
     * @return double value of the nuclear integral
     */
    inline double nuclear_gto(const GTO &gto1, const GTO &gto2, double cx, double cy, double cz) const {
        return this->nuclear_gto(gto1, gto2, Vec3(cx, cy, cz));
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
    double repulsion(const CGF &cgf1, const CGF &cgf2, const CGF &cgf3, const CGF &cgf4) const;

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
    double repulsion(const GTO &gto1, const GTO &gto2, const GTO &gto3, const GTO &gto4) const;

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
    double repulsion_deriv(const CGF &cgf1,const CGF &cgf2,const CGF &cgf3,const CGF &cgf4,
        const Vec3& nucleus, unsigned int coord) const;

    inline double repulsion_deriv(const CGF &cgf1,const CGF &cgf2,const CGF &cgf3,const CGF &cgf4,
        double cx, double cy, double cz, unsigned int coord) const {
        return this->repulsion_deriv(cgf1, cgf2, cgf3, cgf4, Vec3(cx, cy, cz), coord);
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
    double repulsion_deriv(const GTO& gto1, const GTO& gto2, const GTO &gto3, const GTO &gto4, unsigned int coord) const;

    size_t teindex(size_t i, size_t j, size_t k, size_t l) const;

private:
    /*
     * @var     gamma_inc
     * @brief   class that handles the evaluation of the Gamma function
     */
    GammaInc gamma_inc;

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
    double overlap(double alpha1, unsigned int l1, unsigned int m1, unsigned int n1, const Vec3 &a,
                   double alpha2, unsigned int l2, unsigned int m2, unsigned int n2, const Vec3 &b) const;

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
    double dipole(double alpha1, unsigned int l1, unsigned int m1, unsigned int n1, const Vec3 &a,
                  double alpha2, unsigned int l2, unsigned int m2, unsigned int n2, const Vec3 &b,
                  unsigned int cc, double cref = 0.0) const;

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
     * @param Vec3 c            Nuclear position
     * @param coord             Cartesian direction to derive nuclear coordinate towards
     *
     * @return double value of the nuclear integral derived towards nuclear coordinate
     */
    double nuclear_deriv_op(const Vec3& a, int l1, int m1, int n1, double alpha1,
                            const Vec3& b, int l2, int m2, int n2,
                            double alpha2, const Vec3& c, unsigned int coord) const;

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
    double repulsion(const Vec3 &a, const int la, const int ma, const int na, const double alphaa,
                     const Vec3 &b, const int lb, const int mb, const int nb, const double alphab,
                     const Vec3 &c, const int lc, const int mc, const int nc, const double alphac,
                     const Vec3 &d, const int ld, const int md, const int nd, const double alphad) const;

    /**
     * @brief Performs nuclear integral evaluation including caching of Fgamma
     *
     * This function uses function-level caching of the Fgamma function; this implementation
     * was suggested in https://github.com/ifilot/hfcxx/issues/8, but explicit unit testing
     * actually shows not appreciable difference in speed.
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
    double repulsion_fgamma_cached(const Vec3 &a, const int la, const int ma, const int na, const double alphaa,
                                   const Vec3 &b, const int lb, const int mb, const int nb, const double alphab,
                                   const Vec3 &c, const int lc, const int mc, const int nc, const double alphac,
                                   const Vec3 &d, const int ld, const int md, const int nd, const double alphad) const;

    /**
     * @brief Calculates one dimensional overlap integral
     *
     * @param int l1        Power of 'x' component of the polynomial of the first GTO
     * @param int l2        Power of 'x' component of the polynomial of the second GTO
     * @param double x1     'x' component of the position of the first GTO
     * @param double x2     'x' component of the position of the second GTO
     * @param double gamma  Sum of the two Gaussian exponents
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
     * @param double alpha1     Gaussian exponent of the first GTO
     * @param double alpha2     Gaussian exponent of the second GTO
     * @param const Vec3 a      Center of the first GTO
     * @param const Vec3 b      Center of the second GTO
     *
     *
     * @return new gaussian product center
     */
    Vec3 gaussian_product_center(double alpha1, const Vec3 &a,
                                 double alpha2, const Vec3 &b) const;

    double binomial_prefactor(int s, int ia, int ib, double xpa, double xpb) const;

    double binomial(int a, int b) const;

    std::vector<double> A_array(const int l1, const int l2,
                                const double pa, const double pb,
                                const double cp, const double g) const;

    std::vector<double> A_array_deriv(const int l1, const int l2,
                                      const double pa, const double pb,
                                      const double cp, const double g) const;

    double A_term(const int i, const int r, const int u,
                  const int l1, const int l2,
                  const double pax, const double pbx,
                  const double cpx, const double gamma) const;

    double gamma(const double m, double x) const;

    std::vector<double> B_array(const int l1,const int l2,const int l3,const int l4,
                                const double p, const double a, const double b, const double q, const double c, const double d,
                                const double g1, const double g2, const double delta) const;

    double B_term(const int i1, const int i2, const int r1, const int r2, const int u, const int l1, const int l2, const int l3, const int l4,
    const double px, const double ax, const double bx, const double qx, const double cx, const double dx, const double gamma1,
    const double gamma2, const double delta) const;

    double fB(const int i, const int l1, const int l2, const double p, const double a, const double b, const int r, const double q) const;
    double BB0(int i, int r, double q) const;
    double fact_ratio2(unsigned int a, unsigned int b) const;
};
