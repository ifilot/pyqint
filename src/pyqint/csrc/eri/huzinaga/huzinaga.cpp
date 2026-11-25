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

#include "huzinaga.h"

#include "factorials.h"
#include "gaussian_integrals.h"

namespace integrals::eri::huzinaga {

namespace {

/**
 * @brief Compute factorial ratio used in recursion relations.
 *
 * @param a  first index
 * @param b  second index
 *
 * @return factorial ratio
 */
double fact_ratio2(unsigned int a, unsigned int b) {
    return factorial(a) / factorial(b) / factorial(a - 2*b);
}

/**
 * @brief Evaluate the BB0 helper term.
 *
 * @param i  Hermite index
 * @param r  Recursion index
 * @param g  Gaussian exponent sum
 *
 * @return BB0 term value
 */
double BB0(int i, int r, double g) {
    return fact_ratio2(i,r) * std::pow(4*g,r-i);
}

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
 * @param g   Gaussian exponent sum for the second pair
 *
 * @return fB term value
 */
double fB(const int i, const int l1, const int l2, const double p, const double a, const double b, const int r, const double g) {
    return gaussian::binomial_prefactor(i, l1, l2, p-a, p-b) * BB0(i, r, g);
}

/**
 * @brief Evaluate a single B-term for the Huzinaga expansion.
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
    const double gamma2, const double delta) {
    return fB(i1,l1,l2,px,ax,bx,r1,gamma1)*
        std::pow(-1,i2) * fB(i2,l3,l4,qx,cx,dx,r2,gamma2)*
        std::pow(-1,u)*fact_ratio2(i1+i2-2*(r1+r2),u)*
        std::pow(qx-px,i1+i2-2*(r1+r2)-2*u)/
        std::pow(delta,i1+i2-2*(r1+r2)-u);
}

/**
 * @brief Build the B array for the Huzinaga expansion.
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
        const double g1, const double g2, const double delta) {

    int imax = l1 + l2 + l3 + l4 + 1;
    std::vector<double> arrB(imax,0);

    for(int i1=0; i1<l1+l2+1; i1++) {
        for(int i2=0; i2<l3+l4+1; i2++) {
            for(int r1=0; r1 < i1/2+1; r1++) {
                for(int r2=0; r2 < i2/2+1; r2++) {
                    for(int u=0; u<(i1+i2)/2-r1-r2+1; u++) {
                        int i = i1+i2-2*(r1+r2)-u;
                        arrB[i] += B_term(i1,i2,r1,r2,u,l1,l2,l3,l4,
                                                            p,a,b,q,c,d,g1,g2,delta);
                    }
                }
            }
        }
    }
    return arrB;
}

}  // namespace

double evaluate_scalar(
    const GTO& gto_a, const GTO& gto_b,
    const GTO& gto_c, const GTO& gto_d,
    const BoysFunction& boys_function) {
    const Vec3& a = gto_a.get_position();
    const Vec3& b = gto_b.get_position();
    const Vec3& c = gto_c.get_position();
    const Vec3& d = gto_d.get_position();

    const double alphaa = gto_a.get_alpha();
    const double alphab = gto_b.get_alpha();
    const double alphac = gto_c.get_alpha();
    const double alphad = gto_d.get_alpha();

    const int la = gto_a.get_l(), ma = gto_a.get_m(), na = gto_a.get_n();
    const int lb = gto_b.get_l(), mb = gto_b.get_m(), nb = gto_b.get_n();
    const int lc = gto_c.get_l(), mc = gto_c.get_m(), nc = gto_c.get_n();
    const int ld = gto_d.get_l(), md = gto_d.get_m(), nd = gto_d.get_n();

    double rab2 = (a-b).norm2();
    double rcd2 = (c-d).norm2();

    Vec3 p = gaussian::gaussian_product_center(alphaa, a, alphab, b);
    Vec3 q = gaussian::gaussian_product_center(alphac, c, alphad, d);

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
                sum += bx[i]*by[j]*bz[k] * boys_function(i+j+k,0.25*rpq2/delta);
            }
        }
    }

    return 2.0 * std::pow(M_PI,2.5)/(gamma1*gamma2*std::sqrt(gamma1+gamma2))*
                 std::exp(-alphaa*alphab*rab2/gamma1)*
                 std::exp(-alphac*alphad*rcd2/gamma2)*sum;
}

double evaluate_block(
    const GTO& gto_a, const GTO& gto_b,
    const GTO& gto_c, const GTO& gto_d,
    const BoysFunction& boys_function) {
    const Vec3& a = gto_a.get_position();
    const Vec3& b = gto_b.get_position();
    const Vec3& c = gto_c.get_position();
    const Vec3& d = gto_d.get_position();

    const double alphaa = gto_a.get_alpha();
    const double alphab = gto_b.get_alpha();
    const double alphac = gto_c.get_alpha();
    const double alphad = gto_d.get_alpha();

    const int la = gto_a.get_l(), ma = gto_a.get_m(), na = gto_a.get_n();
    const int lb = gto_b.get_l(), mb = gto_b.get_m(), nb = gto_b.get_n();
    const int lc = gto_c.get_l(), mc = gto_c.get_m(), nc = gto_c.get_n();
    const int ld = gto_d.get_l(), md = gto_d.get_m(), nd = gto_d.get_n();

    const double rab2 = (a-b).norm2();
    const double rcd2 = (c-d).norm2();

    const Vec3 p = gaussian::gaussian_product_center(alphaa, a, alphab, b);
    const Vec3 q = gaussian::gaussian_product_center(alphac, c, alphad, d);

    const double rpq2 = (p-q).norm2();

    const double gamma1 = alphaa + alphab;
    const double gamma2 = alphac + alphad;
    const double delta = 0.25 * (1.0 / gamma1 + 1.0 / gamma2);
    const double T = 0.25*rpq2/delta;

    std::vector<double> bx = B_array(la, lb, lc, ld, p[0], a[0], b[0], q[0], c[0], d[0], gamma1, gamma2, delta);
    std::vector<double> by = B_array(ma, mb, mc, md, p[1], a[1], b[1], q[1], c[1], d[1], gamma1, gamma2, delta);
    std::vector<double> bz = B_array(na, nb, nc, nd, p[2], a[2], b[2], q[2], c[2], d[2], gamma1, gamma2, delta);

    // The Cartesian terms reuse Boys values with the same total order.
    const int lmax = la + lb + lc + ld + ma + mb + mc + md + na + nb + nc + nd;
    std::vector<double> fg(la+lb+lc+ld+ma+mb+mc+md+na+nb+nc+nd+1);
    boys_function.compute_block(lmax, T, fg.data());

    double sum = 0.0;
    for(int i=0; i<=(la+lb+lc+ld); i++) {
        for(int j=0; j<=(ma+mb+mc+md); j++) {
            for(int k=0; k<=(na+nb+nc+nd); k++) {
                sum += bx[i]*by[j]*bz[k]*fg[i+j+k];
            }
        }
    }

    return 2.0 * std::pow(M_PI,2.5)/(gamma1*gamma2*std::sqrt(gamma1+gamma2))*
                 std::exp(-alphaa*alphab*rab2/gamma1)*
                 std::exp(-alphac*alphad*rcd2/gamma2)*sum;
}

}  // namespace integrals::eri::huzinaga
