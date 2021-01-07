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

#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <unordered_map>

#include "gamma.h"
#include "cgf.h"

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

    inline const char* get_compiler_version() const {
        return this->compiler_version.c_str();
    }

    inline const char* get_compile_time() const {
        return this->compile_time.c_str();
    }

    inline const char* get_compile_date() const {
        return this->compile_date.c_str();
    }

    inline const char* get_openmp_version() const {
        return this->openmp_version.c_str();
    }

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
     * @brief Calculates overlap integral of two CGF
     *
     * @param const CGF& cgf1   Contracted Gaussian Function
     * @param const CGF& cgf2   Contracted Gaussian Function
     *********************
     * Calculates the value of < cgf1 | cgf2 >
     *
     * @return double value of the overlap integral
     */
    double overlap(const CGF& cgf1, const CGF& cgf2) const;

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
    double overlap(const GTO& gto1, const GTO& gto2) const;

    /**
     * @brief Calculates derivative of overlap integral of two CGF
     *
     * @param const CGF& cgf1       Contracted Gaussian Function
     * @param const CGF& cgf2       Contracted Gaussian Function
     * @param unsigned int coord    Derivative direction
     *
     * Calculates the value of d/dx < cgf1 | cgf2 >
     *
     * @return double value of the nuclear integral
     */
    double overlap_deriv(const CGF& cgf1, const CGF& cgf2, const vec3& nucleus, unsigned int coord) const;

    // expanded interface for Cython
    inline double overlap_deriv(const CGF& cgf1, const CGF& cgf2, double cx, double cy, double cz, unsigned int coord) const {
        return this->overlap_deriv(cgf1, cgf2, vec3(cx, cy, cz), coord);
    }

    /**
     * @brief Calculates derivative of the overlap
     * integral (of the first BF) towards coordinate
     *
     * @param const GTO& gto1       Gaussian Type Orbital
     * @param const GTO& gto2       Gaussian Type Orbital
     * @param unsigned int coord    Derivative direction
     *
     * Calculates the value of < d/d1x gto1 | gto2 >
     *
     * @return double value of the overlap integral
     */
    double overlap_deriv(const GTO& gto1, const GTO& gto2, unsigned int coord) const;

    /**
     * @brief Calculates kinetic integral of two CGF
     *
     * @param const CGF& cgf1   Contracted Gaussian Function
     * @param const CGF& cgf2   Contracted Gaussian Function
     *
     * Calculates the value of < cgf1 | T | cgf2 >
     *
     * @return double value of the kinetic integral
     */
    double kinetic(const CGF& cgf1, const CGF& cgf2) const;

    /**
     * @brief Calculates kinetic integral of two GTO
     *
     * @param const GTO& gto1   Gaussian Type Orbital
     * @param const GTO& gto2   Gaussian Type Orbital
     *
     * Calculates the value of < gto1 | H | gto2 >
     *
     * @return double value of the kinetic integral
     */
    double kinetic(const GTO& gto1, const GTO& gto2) const;

    /**
     * @brief Calculates derivative of overlap integral of two CGF
     *
     * @param const CGF& cgf1       Contracted Gaussian Function
     * @param const CGF& cgf2       Contracted Gaussian Function
     * @param const vec3& nucleus   Nucleus coordinates
     * @param unsigned int coord    Derivative direction
     *
     * Calculates the value of d/dcx < cgf1 | -1/2 nabla^2 | cgf2 >
     *
     * @return double value of the nuclear integral
     */
    double kinetic_deriv(const CGF& cgf1, const CGF& cgf2, const vec3& nucleus, unsigned int coord) const;

    // expanded interface for Cython
    inline double kinetic_deriv(const CGF& cgf1, const CGF& cgf2, double cx, double cy, double cz, unsigned int coord) const {
        return this->kinetic_deriv(cgf1, cgf2, vec3(cx, cy, cz), coord);
    }

    /**
     * @brief Calculates derivative of kinetic integral of two GTO
     *
     * @param const GTO& gto1       Gaussian Type Orbital
     * @param const GTO& gto2       Gaussian Type Orbital
     * @param unsigned int coord    Derivative direction
     *
     * Calculates the value of < d/dx gto1 |-1/2 nabla^2 | gto2 >
     *
     * @return double value of the overlap integral
     */
    double kinetic_deriv(const GTO& gto1, const GTO& gto2, unsigned int coord) const;

    /**
     * @brief Calculates nuclear integral of two CGF
     *
     * @param const CGF& cgf1       Contracted Gaussian Function
     * @param const CGF& cgf2       Contracted Gaussian Function
     * @param const vec3 nucleus    Position of the nucleus
     * @param unsigned int charge   charge of the nucleus in a.u.
     *
     * Calculates the value of < cgf1 | V | cgf2 >
     *
     * @return double value of the nuclear integral
     */
    double nuclear(const CGF &cgf1, const CGF &cgf2, const vec3& nucleus, unsigned int charge) const;

    // expanded notation for Cython interface
    inline double nuclear(const CGF &cgf1, const CGF &cgf2, double cx, double cy, double cz, unsigned int charge) const {
        return this->nuclear(cgf1, cgf2, vec3(cx, cy, cz), charge);
    }

    /**
     * @brief Calculates nuclear integral of two CGF
     *
     * @param const CGF& cgf1       Contracted Gaussian Function
     * @param const CGF& cgf2       Contracted Gaussian Function
     * @param const vec3 nucleus    Position of the nucleus
     * @param unsigned int charge   charge of the nucleus in a.u.
     *
     * Calculates the value of < cgf1 | V | cgf2 >
     *
     * @return double value of the nuclear integral
     */
    double nuclear_deriv(const CGF &cgf1, const CGF &cgf2, const vec3& nucleus,
        unsigned int charge, unsigned int coord) const;

    // expanded notation for Cython interface
    inline double nuclear_deriv(const CGF &cgf1, const CGF &cgf2, double cx, double cy, double cz,
        unsigned int charge, unsigned int coord) const {
        return this->nuclear_deriv(cgf1, cgf2, vec3(cx, cy, cz), charge, coord);
    }

    /**
     * @brief Calculates nuclear integral of two CGF
     *
     * @param const GTO& gto1       Contracted Gaussian Function
     * @param const GTO& gto2       Contracted Gaussian Function
     * @param const vec3 nucleus    Position of the nucleus
     * @param unsigned int charge   charge of the nucleus in a.u.
     *
     * Calculates the value of < gto1 | V | gto2 >
     *
     * @return double value of the nuclear integral
     */
    double nuclear(const GTO &gto1, const GTO &gto2, const vec3& nucleus) const;

    /**
     * @brief Calculates derivative towards one of the basis functions in the
     * nuclear integral of two GTOs
     *
     * @param const GTO& gto1       Contracted Gaussian Function
     * @param const GTO& gto2       Contracted Gaussian Function
     * @param const vec3 nucleus    Position of the nucleus
     * @param unsigned int charge   charge of the nucleus in a.u.
     * @param unsigned int coord    Cartesian direction to take derivative towards
     *
     * Calculates the value of < d/dx1 gto1 | V | gto2 >
     *
     * @return double value of the nuclear integral
     */
    double nuclear_deriv_bf(const GTO &gto1, const GTO &gto2, const vec3& nucleus, unsigned int coord) const;

    /**
     * @brief Calculates derivative of the nuclear attraction operator
     * in a nuclear attraction integral of two CGF
     *
     * @param const GTO& gto1       Contracted Gaussian Function
     * @param const GTO& gto2       Contracted Gaussian Function
     * @param const vec3 nucleus    Position of the nucleus
     * @param unsigned int charge   charge of the nucleus in a.u.
     * @param unsigned int coord    Cartesian direction to take derivative towards
     *
     * Calculates the value of < gto1 | d/dcx V | gto2 >
     *
     * @return double value of the nuclear integral
     */
    inline double nuclear_deriv_bf(const GTO &gto1, const GTO &gto2, double rx, double ry, double rz, unsigned int coord) const {
        return this->nuclear_deriv_bf(gto1, gto2, vec3(rx, ry, rz), coord);
    }

    /**
     * @brief Calculates nuclear integral of two CGF
     *
     * @param const GTO& gto1       Contracted Gaussian Function
     * @param const GTO& gto2       Contracted Gaussian Function
     * @param const vec3 nucleus    Position of the nucleus
     * @param unsigned int charge   charge of the nucleus in a.u.
     * @param unsigned int coord    Cartesian direction to take derivative towards
     *
     * Calculates the value of < gto1 | V | gto2 >
     *
     * @return double value of the nuclear integral
     */
    inline double nuclear_deriv_op(const GTO &gto1, const GTO &gto2, const vec3& nucleus, unsigned int coord) const;

    /**
     * @brief Calculates derivative of the nuclear attraction operator
     * in a nuclear attraction integral of two CGF
     *
     * @param const GTO& gto1       Contracted Gaussian Function
     * @param const GTO& gto2       Contracted Gaussian Function
     * @param const vec3 nucleus    Position of the nucleus
     * @param unsigned int charge   charge of the nucleus in a.u.
     * @param unsigned int coord    Cartesian direction to take derivative towards
     *
     * Calculates the value of < gto1 | d/dcx V | gto2 >
     *
     * @return double value of the nuclear integral
     */
    inline double nuclear_deriv_op(const GTO &gto1, const GTO &gto2, double rx, double ry, double rz, unsigned int coord) const {
        return this->nuclear_deriv_op(gto1, gto2, vec3(rx, ry, rz), coord);
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
    inline double nuclear(const GTO &gto1, const GTO &gto2, double cx, double cy, double cz) const {
        return this->nuclear(gto1, gto2, vec3(cx, cy, cz));
    }

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
     * @param const vec3& nucleus   Nucleus coordinates
     * @param unsigned int coord    Derivative direction
     *
     * Calculates the value of d/dcx < cgf1 | cgf2 | cgf3 | cgf4 >
     *
     * @return double value of the repulsion integral
     */
    double repulsion_deriv(const CGF &cgf1,const CGF &cgf2,const CGF &cgf3,const CGF &cgf4,
        const vec3& nucleus, unsigned int coord) const;

    inline double repulsion_deriv(const CGF &cgf1,const CGF &cgf2,const CGF &cgf3,const CGF &cgf4,
        double cx, double cy, double cz, unsigned int coord) const {
        return this->repulsion_deriv(cgf1, cgf2, cgf3, cgf4, vec3(cx, cy, cz), coord);
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

    const unsigned int teindex(unsigned int i, unsigned int j, unsigned int k, unsigned int l) const;

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
     * @param vec3 a            Center of the Gaussian orbital of the first GTO
     * @param double alpha2     Gaussian exponent of the second GTO
     * @param unsigned int l2   Power of x component of the polynomial of the second GTO
     * @param unsigned int m2   Power of y component of the polynomial of the second GTO
     * @param unsigned int n2   Power of z component of the polynomial of the second GTO
     * @param vec3 b            Center of the Gaussian orbital of the second GTO
     *
     * @return double value of the overlap integral
     */
    double overlap(double alpha1, unsigned int l1, unsigned int m1, unsigned int n1, const vec3 &a,
                   double alpha2, unsigned int l2, unsigned int m2, unsigned int n2, const vec3 &b) const;

    /**
     * @brief Performs nuclear integral evaluation
     *
     * @param vec3 a            Center of the Gaussian orbital of the first GTO
     * @param unsigned int l1   Power of x component of the polynomial of the first GTO
     * @param unsigned int m1   Power of y component of the polynomial of the first GTO
     * @param unsigned int n1   Power of z component of the polynomial of the first GTO
     * @param double alpha1     Gaussian exponent of the first GTO
     * @param vec3 b            Center of the Gaussian orbital of the second GTO
     * @param unsigned int l2   Power of x component of the polynomial of the second GTO
     * @param unsigned int m2   Power of y component of the polynomial of the second GTO
     * @param unsigned int n2   Power of z component of the polynomial of the second GTO
     * @param double alpha2     Gaussian exponent of the second GTO
     * @param vec3 c
     *
     * @return double value of the nuclear integral
     */
    double nuclear(const vec3& a,
                   int l1, int m1, int n1,
                   double alpha1,
                   const vec3& b,
                   int l2, int m2, int n2,
                   double alpha2,
                   const vec3& c) const;

    /**
     * @brief Performs nuclear integral evaluation
     *
     * @param vec3 a            Center of the Gaussian orbital of the first GTO
     * @param unsigned int l1   Power of x component of the polynomial of the first GTO
     * @param unsigned int m1   Power of y component of the polynomial of the first GTO
     * @param unsigned int n1   Power of z component of the polynomial of the first GTO
     * @param double alpha1     Gaussian exponent of the first GTO
     * @param vec3 b            Center of the Gaussian orbital of the second GTO
     * @param unsigned int l2   Power of x component of the polynomial of the second GTO
     * @param unsigned int m2   Power of y component of the polynomial of the second GTO
     * @param unsigned int n2   Power of z component of the polynomial of the second GTO
     * @param double alpha2     Gaussian exponent of the second GTO
     * @param vec3 c            Nuclear position
     * @param coord             Cartesian direction to derive nuclear coordinate towards
     *
     * @return double value of the nuclear integral derived towards nuclear coordinate
     */
    double nuclear_deriv_op(const vec3& a, int l1, int m1, int n1, double alpha1,
                            const vec3& b, int l2, int m2, int n2,
                            double alpha2, const vec3& c, unsigned int coord) const;

    /**
     * @brief Performs nuclear integral evaluation
     *
     * @param vec3 a            Center of the Gaussian orbital of the first GTO
     * @param unsigned int l1   Power of x component of the polynomial of the first GTO
     * @param unsigned int m1   Power of y component of the polynomial of the first GTO
     * @param unsigned int n1   Power of z component of the polynomial of the first GTO
     * @param double alpha1     Gaussian exponent of the first GTO
     * @param vec3 b            Center of the Gaussian orbital of the second GTO
     * @param unsigned int l2   Power of x component of the polynomial of the second GTO
     * @param unsigned int m2   Power of y component of the polynomial of the second GTO
     * @param unsigned int n2   Power of z component of the polynomial of the second GTO
     * @param double alpha2     Gaussian exponent of the second GTO
     * @param vec3 c
     *
     * @return double value of the nuclear integral
     */
    double repulsion(const vec3 &a, const double norma, const int la, const int ma, const int na, const double alphaa,
                     const vec3 &b, const double normb, const int lb, const int mb, const int nb, const double alphab,
                     const vec3 &c, const double normc, const int lc, const int mc, const int nc, const double alphac,
                     const vec3 &d, const double normd, const int ld, const int md, const int nd, const double alphad) const;

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
     * @param const vec3 a      Center of the first GTO
     * @param const vec3 b      Center of the second GTO
     *
     *
     * @return new gaussian product center
     */
    vec3 gaussian_product_center(double alpha1, const vec3 &a,
                                 double alpha2, const vec3 &b) const;

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
    double B0(int i, int r, double q) const;
    double fact_ratio2(const int a, const int b) const;

    void init();
};
