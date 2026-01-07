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
#include <cassert>

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

/*
 * @fn get_amp
 * @brief Gets the amplitude of the GTO
 *
 * @param Vec3 r    coordinates
 *
 * @return const double amplitude
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

/*
 * @fn get_gradient
 * @brief Gets the gradient of the GTO
 *
 * @param Vec3 r    coordinates
 *
 * @return gradient
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

/*
 * @fn calculate_normalization_constant
 * @brief Calculates the normalization constant so that <GTO|GTO>=1
 *
 * @return void
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

/*
 * @fn CGF
 * @brief Constructor
 *
 * @return CGF
 */
CGF::CGF():
    r(Vec3(0,0,0)) {
        // do nothing
}

/*
 * @fn CGF
 * @brief Default constructor
 *
 * @return CGF
 */
CGF::CGF(double x, double y, double z) :
    r(Vec3(x,y,z)) {}

/*
 * @fn CGF
 * @brief Default constructor
 *
 * @return CGF
 */
CGF::CGF(const Vec3& _r):
    r(_r) {
        // do nothing
}

/*
 * @fn get_amp
 * @brief Gets the amplitude of the CGF
 *
 * @param Vec3 r    coordinates
 *
 * @return const double amplitude
 */
const double CGF::get_amp(const Vec3& r) const noexcept {
    double sum = 0.0;

    for(const auto& gto : this->gtos) {
        sum += gto.get_coefficient() * gto.get_amp(r);
    }

    return sum;
}

/*
 * @fn get_grad
 * @brief Gets the gradient of the CGF
 *
 * @param Vec3 r    coordinates
 *
 * @return gradient
 */
std::vector<double> CGF::get_grad(const Vec3& r) const noexcept {
    Vec3 sum = Vec3(0,0,0);

    for(const auto& gto : this->gtos) {
        sum += gto.get_coefficient() * gto.get_grad(r);
    }

    return std::vector<double>({sum[0], sum[1], sum[2]});
}

/*
 * @fn add_GTO
 * @brief Add a GTO to the CGF
 *
 * @param unsigned int type     type of the orbital (see above for the list)
 * @param double alpha          alpha value
 * @param double c              coefficient
 * @param const Vec3& Vec3      position
 *
 * @return void
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

/*
 * @fn add_gto
 * @brief Add a GTO to the CGF
 *
 * @param double c              coefficient
 * @param double alpha          alpha value
 * @param unsigned int l        l angular momentum x
 * @param unsigned int m        m angular momentum y
 * @param unsigned int n        n angular momentum z
 *
 * @return void
 */
void CGF::add_gto(double c,
                  double alpha,
                  unsigned int l,
                  unsigned int m,
                  unsigned int n) {
    this->gtos.push_back(GTO(c, this->r, alpha, l, m, n));
}

/*
 * @fn add_gto_with_position
 * @brief Add a GTO to the CGF
 *
 * @param double c              coefficient
 * @param double px             px value
 * @param double py             px value
 * @param double pz             px value
 * @param double alpha          alpha value
 * @param unsigned int l        l angular momentum x
 * @param unsigned int m        m angular momentum y
 * @param unsigned int n        n angular momentum z
 *
 * @return void
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

/*
 * @fn set_position
 * @brief Set a (new) center for the CGF
 *
 * @param pos   center of the CGF
 *
 * @return void
 */
void CGF::set_position(const Vec3 &pos) {
    this->r = pos;

    for(unsigned int i=0; i<this->gtos.size(); i++) {
        this->gtos[i].set_position(pos);
    }
}

/*
 * @fn max_primitive_l
 * @brief Get maximum l value among GTOs
 *
 * @return unsigned int maximum l value among GTOs
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

    // Standard contraction: all primitives have same (l,m,n)
    const int L = l + m + n;

    static const double pi_three_half_pow = M_PI * std::sqrt(M_PI);
    const double prefactor = pi_three_half_pow *
                       (l < 1 ? 1.0 : double_factorial(2*l - 1)) *
                       (m < 1 ? 1.0 : double_factorial(2*m - 1)) *
                       (n < 1 ? 1.0 : double_factorial(2*n - 1)) /
                       static_cast<double>(1 << L);

    double sum = 0.0;
    for (size_t i = 0; i < this->size(); ++i) {
        for (size_t j = 0; j < this->size(); ++j) {
            // const double a_i = this->get_coefficient_gto(i);
            // const double a_j = this->get_coefficient_gto(j);
            const double a_i = this->get_norm_gto(i) * this->get_coefficient_gto(i);
            const double a_j = this->get_norm_gto(j) * this->get_coefficient_gto(j);
            const double alpha_i = this->get_gto(i).get_alpha();
            const double alpha_j = this->get_gto(j).get_alpha();

            sum += a_i * a_j / std::pow(alpha_i + alpha_j, L + 1.5);
        }
    }

    return 1.0 / std::sqrt(prefactor * sum);
}
