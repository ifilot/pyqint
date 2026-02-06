#include "integrals.h"

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
double Integrator::repulsion_hellsing (
    const Vec3 &a, const int la, const int ma, const int na, const double alphaa,
    const Vec3 &b, const int lb, const int mb, const int nb, const double alphab,
    const Vec3 &c, const int lc, const int mc, const int nc, const double alphac,
    const Vec3 &d, const int ld, const int md, const int nd, const double alphad) const {

    const double rab2 = (a-b).norm2();
    const double rcd2 = (c-d).norm2();

    const Vec3 p = gaussian_product_center(alphaa, a, alphab, b);
    const Vec3 q = gaussian_product_center(alphac, c, alphad, d);
    const double rpq2 = (p-q).norm2();

    const double gamma1 = alphaa + alphab;
    const double gamma2 = alphac + alphad;
    const double eta    = (gamma1 * gamma2) / (gamma1 + gamma2);
    const double T = eta * rpq2;

    // evaluate nu_ma
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

    // early exit for (ss|ss); saves a couple of cycles
    if(nu_max == 0) {
        return pref * this->boys_function(0, T);
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
    this->boys_function.compute_block(nu_max, T, F.data());

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
double Integrator::repulsion_hellsing_cached(
    const Vec3 &a, const int la, const int ma, const int na, const double alphaa,
    const Vec3 &b, const int lb, const int mb, const int nb, const double alphab,
    const Vec3 &c, const int lc, const int mc, const int nc, const double alphac,
    const Vec3 &d, const int ld, const int md, const int nd, const double alphad) const
{
    // ---------------- geometry ----------------
    const double rab2 = (a - b).norm2();
    const double rcd2 = (c - d).norm2();

    const Vec3 p = gaussian_product_center(alphaa, a, alphab, b);
    const Vec3 q = gaussian_product_center(alphac, c, alphad, d);
    const double rpq2 = (p - q).norm2();

    // ---------------- scalars ----------------
    const double gamma1 = alphaa + alphab;
    const double gamma2 = alphac + alphad;
    const double eta    = (gamma1 * gamma2) / (gamma1 + gamma2);
    const double T      = eta * rpq2;

    const int nu_max =
        (la + ma + na) +
        (lb + mb + nb) +
        (lc + mc + nc) +
        (ld + md + nd);

    // ---------------- universal prefactor ----------------
    constexpr double PI25 = 17.493418327624862846262821679871;
    const double pref =
        2.0 * PI25 / (gamma1 * gamma2 * std::sqrt(gamma1 + gamma2)) *
        std::exp(-alphaa * alphab * rab2 / gamma1) *
        std::exp(-alphac * alphad * rcd2 / gamma2);

    // ---------------- Boys block ----------------
    std::vector<double> F(nu_max + 1);
    this->boys_function.compute_block(nu_max, T, F.data());

    // ---------------- cached kernels ----------------
    const auto& Kx = hellsing_cache.get(la, lb, lc, ld);
    const auto& Ky = hellsing_cache.get(ma, mb, mc, md);
    const auto& Kz = hellsing_cache.get(na, nb, nc, nd);

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

            // this is exactly: mu - u
            nu[i] = static_cast<int>(K.mu[i]) - static_cast<int>(K.u[i]);
        }
    };

    std::vector<double> Bx, By, Bz;
    std::vector<int>    nux, nuy, nuz;

    build_poly(Kx, a[0] - b[0], c[0] - d[0], p[0] - q[0], Bx, nux);
    build_poly(Ky, a[1] - b[1], c[1] - d[1], p[1] - q[1], By, nuy);
    build_poly(Kz, a[2] - b[2], c[2] - d[2], p[2] - q[2], Bz, nuz);

    // ---------------- triple sum (same logic as reference) ----------------
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
std::vector<HellsingBTerm> Integrator::B_array_hellsing(
    int l1, int l2, int l3, int l4,
    double a1, double a2, double a3, double a4,
    double ax, double bx, double cx, double dx,
    double px, double qx,
    double g1, double g2) const {

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
