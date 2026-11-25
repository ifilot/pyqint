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

#include "hellsing.h"

#include "factorials.h"
#include "gaussian_integrals.h"
#include "mathfuncs.h"

namespace integrals::eri::hellsing {

namespace {

struct HellsingBTerm {
    double c;  // coefficient
    int mu;    // mu
    int u;     // u
};

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
    double g1, double g2) {

    const double pre1 =
        sign_pow(l1 + l2) *
        factorial(l1) * factorial(l2) /
        std::pow(g1, l1 + l2);

    const double pre2 =
        factorial(l3) * factorial(l4) /
        std::pow(g2, l3 + l4);

    const double eta = g1 * g2 / (g1 + g2);

    std::vector<HellsingBTerm> arr;

    for (int i1 = 0; i1 <= l1 / 2; ++i1)
    for (int i2 = 0; i2 <= l2 / 2; ++i2)
    for (int o1 = 0; o1 <= l1 - 2 * i1; ++o1)
    for (int o2 = 0; o2 <= l2 - 2 * i2; ++o2)
    for (int r1 = 0; r1 <= (o1 + o2) / 2; ++r1)
    {
        const double t11 =
            sign_pow(o2 + r1) *
            factorial(o1 + o2) /
            ipow(4.0, i1 + i2 + r1) /
            factorial(i1) /
            factorial(i2) /
            factorial(o1) /
            factorial(o2) /
            factorial(r1);

        const double t12 =
            ipow(a1, o2 - i1 - r1) *
            ipow(a2, o1 - i2 - r1) *
            ipow(g1, 2 * (i1 + i2) + r1) *
            ipow(ax - bx, o1 + o2 - 2 * r1) /
            factorial(l1 - 2 * i1 - o1) /
            factorial(l2 - 2 * i2 - o2) /
            factorial(o1 + o2 - 2 * r1);

        for (int i3 = 0; i3 <= l3 / 2; ++i3)
        for (int i4 = 0; i4 <= l4 / 2; ++i4)
        for (int o3 = 0; o3 <= l3 - 2 * i3; ++o3)
        for (int o4 = 0; o4 <= l4 - 2 * i4; ++o4)
        for (int r2 = 0; r2 <= (o3 + o4) / 2; ++r2)
        {
            const double t21 =
                sign_pow(o3 + r2) *
                factorial(o3 + o4) /
                ipow(4.0, i3 + i4 + r2) /
                factorial(i3) /
                factorial(i4) /
                factorial(o3) /
                factorial(o4) /
                factorial(r2);

            const double t22 =
                ipow(a3, o4 - i3 - r2) *
                ipow(a4, o3 - i4 - r2) *
                ipow(g2, 2 * (i3 + i4) + r2) *
                ipow(cx - dx, o3 + o4 - 2 * r2) /
                factorial(l3 - 2 * i3 - o3) /
                factorial(l4 - 2 * i4 - o4) /
                factorial(o3 + o4 - 2 * r2);

            const int mu =
                l1 + l2 + l3 + l4
                - 2 * (i1 + i2 + i3 + i4)
                - (o1 + o2 + o3 + o4);

            for (int u = 0; u <= mu / 2; ++u)
            {
                const double t3 =
                    sign_pow(u) *
                    factorial(mu) *
                    ipow(eta, mu - u) *
                    ipow(px - qx, mu - 2 * u) /
                    ipow(4.0, u) /
                    factorial(u) /
                    factorial(mu - 2 * u);

                arr.push_back({
                    pre1 * pre2 * t11 * t12 * t21 * t22 * t3,
                    mu,
                    u
                });
            }
        }
    }

    return arr;
}

}  // namespace

double evaluate_direct(
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
    const double eta    = (gamma1 * gamma2) / (gamma1 + gamma2);
    const double T = eta * rpq2;

    const int nu_max =
        (la + ma + na) +
        (lb + mb + nb) +
        (lc + mc + nc) +
        (ld + md + nd);

    // universal pre-factor
    constexpr double PI25 = 17.493418327624862846262821679871;
    const double pref = 2.0 * PI25 / (gamma1 * gamma2 * std::sqrt(gamma1 + gamma2)) *
                        std::exp(-alphaa * alphab * rab2 / gamma1) *
                        std::exp(-alphac * alphad * rcd2 / gamma2);

    // An (ss|ss) quartet only needs the zeroth Boys order.
    if (nu_max == 0) {
        return pref * boys_function(0, T);
    }

    // Build B-arrays for each Cartesian component (x, y, z)
    const auto Bx = B_array_hellsing(
        la, lb, lc, ld,
        alphaa, alphab, alphac, alphad,
        a[0], b[0], c[0], d[0],
        p[0], q[0],
        gamma1, gamma2);

    const auto By = B_array_hellsing(
        ma, mb, mc, md,
        alphaa, alphab, alphac, alphad,
        a[1], b[1], c[1], d[1],
        p[1], q[1],
        gamma1, gamma2);

    const auto Bz = B_array_hellsing(
        na, nb, nc, nd,
        alphaa, alphab, alphac, alphad,
        a[2], b[2], c[2], d[2],
        p[2], q[2],
        gamma1, gamma2);

    // build Fgamma block
    std::vector<double> F(nu_max + 1);
    boys_function.compute_block(nu_max, T, F.data());

    // triple sum
    double s = 0.0;
    for (const auto& tx : Bx) {
        for (const auto& ty : By) {
            for (const auto& tz : Bz) {
                const int nu = (tx.mu + ty.mu + tz.mu) - (tx.u + ty.u + tz.u);
                s += tx.c * ty.c * tz.c * F[nu];
            }
        }
    }

    return pref * s;
}

}  // namespace integrals::eri::hellsing
