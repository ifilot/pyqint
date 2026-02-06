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

#include "cgf.h"
#include "gaussian_integrals.h"
#include <cassert>

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
GTO::GTO(double _c,
         const Vec3& _position,     // position (unit = Bohr)
         double _alpha,
         unsigned int _l,
         unsigned int _m,
         unsigned int _n):
    c(_c),
    alpha(_alpha),
    l(_l),
    m(_m),
    n(_n),
    position(_position) {

    // calculate the normalization constant
    this->calculate_normalization_constant();
}

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
GTO::GTO(double _c,
         double _x,
         double _y,
         double _z,
         double _alpha,
         unsigned int _l,
         unsigned int _m,
         unsigned int _n):
    c(_c),
    alpha(_alpha),
    l(_l),
    m(_m),
    n(_n),
    position(Vec3(_x,_y,_z)) {

    // calculate the normalization constant
    this->calculate_normalization_constant();
}

/**
 * @brief Get the amplitude of the GTO.
 *
 * @param r    coordinates
 *
 * @return amplitude value
 */
const double GTO::get_amp(const Vec3& r) const noexcept {
    const double dx = r[0] - this->position[0];
    const double dy = r[1] - this->position[1];
    const double dz = r[2] - this->position[2];
    const double r2 = dx*dx + dy*dy + dz*dz;

    return this->norm *
        ipow(dx, this->l) *
        ipow(dy, this->m) *
        ipow(dz, this->n) *
        std::exp(-alpha * r2);
}

/**
 * @brief Get the gradient of the GTO.
 *
 * @param r    coordinates
 *
 * @return gradient vector
 */
Vec3 GTO::get_grad(const Vec3& r) const noexcept {
    // calculate exponential term and its product with the cartesian terms
    // for x,y,z components
    const double dx = r[0] - this->position[0];
    const double dy = r[1] - this->position[1];
    const double dz = r[2] - this->position[2];

    const double ex = std::exp(-this->alpha * dx * dx);
    const double ey = std::exp(-this->alpha * dy * dy);
    const double ez = std::exp(-this->alpha * dz * dz);

    const double fx = ipow(dx, this->l) * ex;
    const double fy = ipow(dy, this->m) * ey;
    const double fz = ipow(dz, this->n) * ez;

    // calculate first derivative of the exponential term
    double gx = -2.0 * this->alpha * (r[0]-this->position[0]) * fx;
    double gy = -2.0 * this->alpha * (r[1]-this->position[1]) * fy;
    double gz = -2.0 * this->alpha * (r[2]-this->position[2]) * fz;

    // if there is a Cartesian component (l,m,n > 0), apply the product rule
    // and add the contribution of this term
    if(this->l > 0) {
        gx += this->l * ipow(r[0] - this->position[0], this->l-1) * ex;
    }
    if(this->m > 0) {
        gy += this->m * ipow(r[1] - this->position[1], this->m-1) * ey;
    }
    if(this->n > 0) {
        gz += this->n * ipow(r[2] - this->position[2], this->n-1) * ez;
    }

    // return vector with derivative towards x,y and z
    return Vec3(this->norm * gx * fy * fz,
                this->norm * fx * gy * fz,
                this->norm * fx * fy * gz);
}

/**
 * @brief Calculate the normalization constant so that <GTO|GTO>=1.
 */
void GTO::calculate_normalization_constant() {
    double nom =   std::pow(2.0, 2.0 * (l + m + n) + 3.0 / 2.0) *
                   std::pow(alpha, (l + m + n) + 3.0 / 2.0);

    double denom = (l < 1 ? 1 : double_factorial(2 * l - 1) )*
                   (m < 1 ? 1 : double_factorial(2 * m - 1) )*
                   (n < 1 ? 1 : double_factorial(2 * n - 1) )*
                   std::pow(M_PI, 3.0 / 2.0);

    this->norm = std::sqrt(nom / denom);
}

/**
 * @brief Construct an empty CGF at the origin.
 */
CGF::CGF():
    r(Vec3(0,0,0)) {
        // do nothing
}

/**
 * @brief Construct a CGF at the provided coordinates.
 *
 * @param x  x coordinate
 * @param y  y coordinate
 * @param z  z coordinate
 */
CGF::CGF(double x, double y, double z) :
    r(Vec3(x,y,z)) {}

/**
 * @brief Construct a CGF at the provided position.
 *
 * @param _r  center position
 */
CGF::CGF(const Vec3& _r):
    r(_r) {
        // do nothing
}

/**
 * @brief Get the amplitude of the CGF.
 *
 * @param r    coordinates
 *
 * @return amplitude value
 */
const double CGF::get_amp(const Vec3& r) const noexcept {
    double sum = 0.0;

    for(const auto& gto : this->gtos) {
        sum += gto.get_coefficient() * gto.get_amp(r);
    }

    return sum;
}

/**
 * @brief Get the gradient of the CGF.
 *
 * @param r    coordinates
 *
 * @return gradient components
 */
std::vector<double> CGF::get_grad(const Vec3& r) const noexcept {
    Vec3 sum = Vec3(0,0,0);

    for(const auto& gto : this->gtos) {
        sum += gto.get_coefficient() * gto.get_grad(r);
    }

    return std::vector<double>({sum[0], sum[1], sum[2]});
}

/**
 * @brief Add a GTO to the CGF.
 *
 * @param type   type of the orbital (see above for the list)
 * @param alpha  alpha value
 * @param c      coefficient
 */
void CGF::add_gto(unsigned int type,  // type of the orbital (see above for the list)
                  double alpha,       // alpha value
                  double c) {         // coefficient

    switch(type) {
        // S ORBITAL
        case GTO_S:
            this->gtos.push_back(GTO(c, this->r, alpha, 0,0,0));
        break;

        // P ORBITALS
        case GTO_PX:
            this->gtos.push_back(GTO(c, this->r, alpha, 1,0,0));
        break;
        case GTO_PY:
            this->gtos.push_back(GTO(c, this->r, alpha, 0,1,0));
        break;
        case GTO_PZ:
            this->gtos.push_back(GTO(c, this->r, alpha, 0,0,1));
        break;

        // D ORBITALS
        case GTO_DX2:
            this->gtos.push_back(GTO(c, this->r, alpha, 2,0,0));
        break;
        case GTO_DXY:
            this->gtos.push_back(GTO(c, this->r, alpha, 1,1,0));
        break;
        case GTO_DXZ:
            this->gtos.push_back(GTO(c, this->r, alpha, 1,0,1));
        break;
        case GTO_DY2:
            this->gtos.push_back(GTO(c, this->r, alpha, 0,2,0));
        break;
        case GTO_DYZ:
            this->gtos.push_back(GTO(c, this->r, alpha, 0,1,1));
        break;
        case GTO_DZ2:
            this->gtos.push_back(GTO(c, this->r, alpha, 0,0,2));
        break;

        default:
            std::cerr << "Undefined orbital type. Exiting..." << std::endl;
            exit(-1);
        break;
    }
}

/**
 * @brief Add a GTO to the CGF.
 *
 * @param c      coefficient
 * @param alpha  alpha value
 * @param l      l angular momentum x
 * @param m      m angular momentum y
 * @param n      n angular momentum z
 */
void CGF::add_gto(double c,
                  double alpha,
                  unsigned int l,
                  unsigned int m,
                  unsigned int n) {
    this->gtos.push_back(GTO(c, this->r, alpha, l, m, n));
}

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
void CGF::add_gto_with_position(double c,
                  double px,
                  double py,
                  double pz,
                  double alpha,
                  unsigned int l,
                  unsigned int m,
                  unsigned int n) {
    this->gtos.push_back(GTO(c, Vec3(px, py, pz), alpha, l, m, n));
}

/**
 * @brief Set a (new) center for the CGF.
 *
 * @param pos  center of the CGF
 */
void CGF::set_position(const Vec3 &pos) {
    this->r = pos;

    for(unsigned int i=0; i<this->gtos.size(); i++) {
        this->gtos[i].set_position(pos);
    }
}

/**
 * @brief Get maximum l value among GTOs.
 *
 * @return maximum l value among GTOs
 */
unsigned int CGF::max_primitive_l() const noexcept {
    unsigned int max_l = 0;

    for (unsigned int i = 0; i < this->size(); ++i) {
        const GTO& g = this->get_gto(i);

        max_l = std::max(max_l, g.get_l());
        max_l = std::max(max_l, g.get_m());
        max_l = std::max(max_l, g.get_n());
    }

    return max_l;
}

/**
 * @brief Get the normalization constant for the pair of CGFs.
 *
 * N < φ_i | φ_i > = 1 => N = ...
 * for φ_i is a CGF with angular momentum shell pair (l,m,n)
 * see: https://arxiv.org/pdf/2007.12057 page 10 for more details
 *
 * @return normalization constant
 */
double CGF::get_contraction_norm() const {
    // NOTE: any GTO is fine for the same shell tuple (lmn)
    assert(gtos.size() != 0LU && "No GTOs found!");

    // Check if all GTOs have the same angular momentum
    const auto& first = this->get_gto(0);
    const auto l = first.get_l();
    const auto m = first.get_m();
    const auto n = first.get_n();

    // NOTE: only needed for the spherical harmonics test; otherwise not needed...
    for (size_t i = 1; i < this->size(); ++i) {
        const auto& gto = this->get_gto(i);
        if (gto.get_l() != l || gto.get_m() != m || gto.get_n() != n) {
            // Mixed angular momentum (e.g., spherical harmonics)
            // Coefficients are pre-normalised; no additional factor needed
            return 1.0;
        }
    }

    double sum = 0.0;
    for (size_t i = 0; i < this->size(); ++i) {
        for (size_t j = 0; j < this->size(); ++j) {
            // const double a_i = this->get_coefficient_gto(i);
            // const double a_j = this->get_coefficient_gto(j);
            const double a_i = this->get_norm_gto(i) * this->get_coefficient_gto(i);
            const double a_j = this->get_norm_gto(j) * this->get_coefficient_gto(j);
            const auto& gto_i = this->get_gto(i);
            const auto& gto_j = this->get_gto(j);

            sum += a_i * a_j *
                   integrals::gaussian::overlap_gto(gto_i.get_alpha(),
                                                    gto_i.get_l(),
                                                    gto_i.get_m(),
                                                    gto_i.get_n(),
                                                    gto_i.get_position(),
                                                    gto_j.get_alpha(),
                                                    gto_j.get_l(),
                                                    gto_j.get_m(),
                                                    gto_j.get_n(),
                                                    gto_j.get_position());
        }
    }

    return 1.0 / std::sqrt(sum);
}
