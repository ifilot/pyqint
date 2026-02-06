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

#include <iostream>
#include <vector>
#include <cstddef>

#include "factorials.h"
#include "mathfuncs.h"
#include "vec3.h"

/*
 * Gaussian Type Orbital
 *
 * N * (x-X)^l * (y-Y)^m * (z-Z)^n * exp(-alpha * r^2)
 *
 * where r = sqrt(x^2 + y^2 + z^2)
 * and N a normalization constant such that <GTO|GTO> = 1
 *
 */

class GTO { // Gaussian Type Orbital
private:
    double c;               // coefficient
    double alpha;           // alpha value in the exponent
    unsigned int l,m,n;     // powers of the polynomial
    Vec3 position;          // position vector (unit = Bohr)
    double norm;            // normalization constant

public:
    /**
     * @brief Default constructor.
     */
    GTO(){}

    /**
     * @brief Construct Gaussian Type Orbital.
     *
     * @param _c        coefficient
     * @param _x        position of the Gaussian
     * @param _y        position of the Gaussian
     * @param _z        position of the Gaussian
     * @param _alpha    alpha value in the exponent
     * @param _l        power of x in the polynomial
     * @param _m        power of y in the polynomial
     * @param _n        power of z in the polynomial
     */
    GTO(double _c,
        double _x,
        double _y,
        double _z,
        double _alpha,
        unsigned int _l,
        unsigned int _m,
        unsigned int _n);

    /**
     * @brief Construct Gaussian Type Orbital.
     *
     * @param _c        coefficient
     * @param _position position of the Gaussian
     * @param _alpha    alpha value in the exponent
     * @param _l        power of x in the polynomial
     * @param _m        power of y in the polynomial
     * @param _n        power of z in the polynomial
     */
    GTO(double _c,
        const Vec3& _position,
        double _alpha,
        unsigned int _l,
        unsigned int _m,
        unsigned int _n);

    /*
     * INLINE GETTERS
     */

    /**
     * @brief Get the coefficient.
     *
     * @return coefficient value
     */
    inline const double get_coefficient() const noexcept {
        return this->c;
    }

    /**
     * @brief Get the alpha value in the polynomial.
     *
     * @return alpha value
     */
    inline const double get_alpha() const noexcept {
        return this->alpha;
    }

    /**
     * @brief Get power of the x component.
     *
     * @return x power
     */
    inline const unsigned int get_l() const noexcept {
        return this->l;
    }

    /**
     * @brief Get power of the y component.
     *
     * @return y power
     */
    inline const unsigned int get_m() const noexcept {
        return this->m;
    }

    /**
     * @brief Get power of the z component.
     *
     * @return z power
     */
    inline const unsigned int get_n() const noexcept {
        return this->n;
    }

    /**
     * @brief Return the normalization constant.
     *
     * @return normalization constant
     */
    inline const double get_norm() const noexcept {
        return this->norm;
    }

    /**
     * @brief Get the center of the Gaussian.
     *
     * @return position vector
     */
    inline const Vec3& get_position() const noexcept {
        return this->position;
    }

    /**
     * @brief Get the amplitude of the GTO.
     *
     * @param r    coordinates
     *
     * @return amplitude value
     */
    const double get_amp(const Vec3& r) const noexcept;

    /**
     * @brief Get the amplitude of the GTO.
     *
     * @param x    x coordinate
     * @param y    y coordinate
     * @param z    z coordinate
     *
     * @return amplitude value
     */
    inline double get_amp(double x, double y, double z) const noexcept {
        return this->get_amp(Vec3(x,y,z));
    }

    /**
     * @brief Get the gradient of the GTO.
     *
     * @param r    coordinates
     *
     * @return gradient vector
     */
    Vec3 get_grad(const Vec3& r) const noexcept;

    /**
     * @brief Set a (new) position of the GTO.
     *
     * @param _position  new center position
     */
    inline void set_position(const Vec3& _position) {
        this->position = _position;
    }

private:
    /**
     * @brief Calculate the normalization constant so that <GTO|GTO>=1.
     */
    void calculate_normalization_constant();
};


class CGF { // Contracted Gaussian Function
private:
    std::vector<GTO> gtos;  // vector holding all gtos
    Vec3 r;                 // position of the CGF

public:
    /**
     * @brief Construct an empty CGF at the origin.
     */
    CGF();

    /**
     * @brief Construct a CGF at the provided coordinates.
     *
     * @param x  x coordinate
     * @param y  y coordinate
     * @param z  z coordinate
     */
    CGF(double x, double y, double z);

    /**
     * @brief Construct a CGF at the provided position.
     *
     * @param _r  center position
     */
    CGF(const Vec3& _r);

    // type of GTOs to add
    enum{
        GTO_S,
        GTO_PX,
        GTO_PY,
        GTO_PZ,
        GTO_DX2,
        GTO_DXY,
        GTO_DXZ,
        GTO_DY2,
        GTO_DYZ,
        GTO_DZ2,

        NUM_GTO
    };

    /**
     * @brief Get the vector position.
     *
     * @return position vector
     */
    inline const Vec3& get_r() const noexcept {
        return r;
    }

    /**
     * @brief Return the length of the contraction.
     *
     * @return length of the contraction
     */
    inline size_t size() const noexcept {
        return this->gtos.size();
    }

    /**
     * @brief Return the normalization constant of a GTO.
     *
     * @param i  ith GTO in the CGF
     *
     * @return normalization constant
     */
    inline const double get_norm_gto(const unsigned int i) const noexcept {
        return this->gtos[i].get_norm();
    }

    /**
     * @brief Return the coefficient of the GTO.
     *
     * @param i  ith GTO in the CGF
     *
     * @return GTO coefficient
     */
    inline const double get_coefficient_gto(const unsigned int i) const noexcept {
        return this->gtos[i].get_coefficient();
    }

    /**
     * @brief Return the GTO.
     *
     * @param i  ith GTO in the CGF
     *
     * @return GTO reference
     */
    inline const GTO& get_gto(const unsigned int i) const noexcept {
        return this->gtos[i];
    }

    /**
     * @brief Get the amplitude of the CGF.
     *
     * @param r    coordinates
     *
     * @return amplitude value
     */
    const double get_amp(const Vec3& r) const noexcept;

    /**
     * @brief Get the amplitude of the CGF.
     *
     * @param x  x coordinate
     * @param y  y coordinate
     * @param z  z coordinate
     *
     * @return amplitude value
     */
    inline double get_amp(double x, double y, double z) const noexcept {
        return this->get_amp(Vec3(x,y,z));
    }

    /**
     * @brief Get the gradient of the CGF.
     *
     * @param r    coordinates
     *
     * @return gradient components
     */
    std::vector<double> get_grad(const Vec3& r) const noexcept;

    /**
     * @brief Get the gradient of the CGF.
     *
     * @param x  x coordinate
     * @param y  y coordinate
     * @param z  z coordinate
     *
     * @return gradient components
     */
    inline std::vector<double> get_grad(double x, double y, double z) const noexcept {
        return this->get_grad(Vec3(x,y,z));
    }

    /**
     * @brief Add a GTO to the CGF.
     *
     * @param type   type of the orbital (see above for the list)
     * @param alpha  alpha value
     * @param c      coefficient
     */
    void add_gto(unsigned int type,
                 double alpha,
                 double c);

    /**
     * @brief Add a GTO to the CGF.
     *
     * @param c      coefficient
     * @param alpha  alpha value
     * @param l      l angular momentum x
     * @param m      m angular momentum y
     * @param n      n angular momentum z
     */
    void add_gto(double c,
                 double alpha,
                 unsigned int l,
                 unsigned int m,
                 unsigned int n);

    /**
     * @brief Add a GTO to the CGF with an explicit position.
     *
     * @param c      coefficient
     * @param px     px value
     * @param py     py value
     * @param pz     pz value
     * @param alpha  alpha value
     * @param l      l angular momentum x
     * @param m      m angular momentum y
     * @param n      n angular momentum z
     */
    void add_gto_with_position(double c,
                 double px,
                 double py,
                 double pz,
                 double alpha,
                 unsigned int l,
                 unsigned int m,
                 unsigned int n);

    /**
     * @brief Set a (new) center for the CGF.
     *
     * @param pos  center of the CGF
     */
    void set_position(const Vec3 &pos);

    /**
     * @brief Get maximum l value among GTOs.
     *
     * @return maximum l value among GTOs
     */
    unsigned int max_primitive_l() const noexcept;

    /**
     * @brief Get the normalization constant for the pair of CGFs.
     *
     * N < φ_i | φ_i > = 1 => N = ...
     * for φ_i is a CGF with angular momentum shell pair (l,m,n)
     * see: https://arxiv.org/pdf/2007.12057 page 10 for more details
     *
     * @return normalization constant
     */
    double get_contraction_norm() const;
};
