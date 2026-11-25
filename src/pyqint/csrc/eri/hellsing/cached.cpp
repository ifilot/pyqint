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

#include "mathfuncs.h"

namespace integrals::eri::hellsing {

double evaluate_cached(
    const GTO& gto_a, const GTO& gto_b,
    const GTO& gto_c, const GTO& gto_d,
    const BoysFunction& boys_function, const HellsingCacheTable1D& cache) {
    const Vec3 ab = gto_a.get_position() - gto_b.get_position();
    const Vec3 cd = gto_c.get_position() - gto_d.get_position();
    return evaluate_prepared(gto_a, gto_b, gto_c, gto_d,
                             prepare_pair(gto_a, gto_b, ab.norm2()),
                             prepare_pair(gto_c, gto_d, cd.norm2()),
                             ab, cd, boys_function, cache);
}

double evaluate_prepared(
    const GTO& gto_a, const GTO& gto_b, const GTO& gto_c, const GTO& gto_d,
    const PrimitivePair& bra, const PrimitivePair& ket, const Vec3& ab, const Vec3& cd,
    const BoysFunction& boys_function, const HellsingCacheTable1D& cache) {
    const double alphaa = gto_a.get_alpha();
    const double alphab = gto_b.get_alpha();
    const double alphac = gto_c.get_alpha();
    const double alphad = gto_d.get_alpha();

    const int la = gto_a.get_l(), ma = gto_a.get_m(), na = gto_a.get_n();
    const int lb = gto_b.get_l(), mb = gto_b.get_m(), nb = gto_b.get_n();
    const int lc = gto_c.get_l(), mc = gto_c.get_m(), nc = gto_c.get_n();
    const int ld = gto_d.get_l(), md = gto_d.get_m(), nd = gto_d.get_n();

    const Vec3 pq = bra.product_center - ket.product_center;
    const double gamma1 = bra.exponent_sum;
    const double gamma2 = ket.exponent_sum;
    const double eta = (gamma1 * gamma2) / (gamma1 + gamma2);
    const double T = eta * pq.norm2();

    const int nu_max =
        (la + ma + na) +
        (lb + mb + nb) +
        (lc + mc + nc) +
        (ld + md + nd);

    // ---------------- universal prefactor ----------------
    constexpr double PI25 = 17.493418327624862846262821679871;
    const double pref =
        2.0 * PI25 / (gamma1 * gamma2 * std::sqrt(gamma1 + gamma2)) *
        bra.attenuation * ket.attenuation;

    // ---------------- Boys block ----------------
    std::vector<double> F(nu_max + 1);
    boys_function.compute_block(nu_max, T, F.data());

    // ---------------- cached kernels ----------------
    const auto& Kx = cache.get(la, lb, lc, ld);
    const auto& Ky = cache.get(ma, mb, mc, md);
    const auto& Kz = cache.get(na, nb, nc, nd);

    // ---------------- build 1D polynomials ----------------
    auto build_poly =
        [&](const HellsingCache1D& K,
            double AB, double CD, double PQ,
            std::vector<double>& poly,
            std::vector<int>& nu)
    {
        const std::size_t n = K.scalar.size();
        poly.resize(n);
        nu.resize(n);

        for (std::size_t i = 0; i < n; ++i) {
            const auto& pwr = K.powers[i];

            poly[i] =
                K.scalar[i] *
                ipow(alphaa, pwr[0]) *
                ipow(alphab, pwr[1]) *
                ipow(alphac, pwr[2]) *
                ipow(alphad, pwr[3]) *
                ipow(gamma1, pwr[4]) *
                ipow(gamma2, pwr[5]) *
                ipow(AB,     pwr[6]) *
                ipow(CD,     pwr[7]) *
                ipow(eta,    pwr[8]) *
                ipow(PQ,     pwr[9]);

            // Each expansion term contributes to Boys order mu - u.
            nu[i] = static_cast<int>(K.mu[i]) - static_cast<int>(K.u[i]);
        }
    };

    std::vector<double> Bx, By, Bz;
    std::vector<int>    nux, nuy, nuz;

    build_poly(Kx, ab[0], cd[0], pq[0], Bx, nux);
    build_poly(Ky, ab[1], cd[1], pq[1], By, nuy);
    build_poly(Kz, ab[2], cd[2], pq[2], Bz, nuz);

    // ---------------- triple sum ----------------
    double s = 0.0;
    for (std::size_t ix = 0; ix < Bx.size(); ++ix) {
        const double cx = Bx[ix];
        const int    nx = nux[ix];

        for (std::size_t iy = 0; iy < By.size(); ++iy) {
            const double cxy = cx * By[iy];
            const int    nxy = nx + nuy[iy];

            for (std::size_t iz = 0; iz < Bz.size(); ++iz) {
                s += cxy * Bz[iz] * F[nxy + nuz[iz]];
            }
        }
    }

    return pref * s;
}

}  // namespace integrals::eri::hellsing
