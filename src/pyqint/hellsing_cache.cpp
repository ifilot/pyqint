#include "hellsing_cache.h"

std::size_t
HellsingCacheTable1D::index(int l1, int l2, int l3, int l4, int lmax)
{
    const std::size_t n = static_cast<std::size_t>(lmax) + 1;
    return ((static_cast<std::size_t>(l1) * n + static_cast<std::size_t>(l2)) * n
          + static_cast<std::size_t>(l3)) * n
          + static_cast<std::size_t>(l4);
}

HellsingCacheTable1D::HellsingCacheTable1D(int lmax)
: lmax_(lmax)
{
    const std::size_t n =
        static_cast<std::size_t>(lmax_ + 1) *
        static_cast<std::size_t>(lmax_ + 1) *
        static_cast<std::size_t>(lmax_ + 1) *
        static_cast<std::size_t>(lmax_ + 1);

    table_.resize(n);

    for (int l1 = 0; l1 <= lmax_; ++l1)
        for (int l2 = 0; l2 <= lmax_; ++l2)
            for (int l3 = 0; l3 <= lmax_; ++l3)
                for (int l4 = 0; l4 <= lmax_; ++l4)
                    table_[index(l1, l2, l3, l4, lmax_)] =
                        build_kernel_1d(l1, l2, l3, l4);
}

const HellsingCache1D&
HellsingCacheTable1D::get(int l1, int l2, int l3, int l4) const
{
    return table_[index(l1, l2, l3, l4, lmax_)];
}

HellsingCache1D
HellsingCacheTable1D::build_kernel_1d(int l1, int l2, int l3, int l4)
{
    HellsingCache1D K;

    // A small (very rough) reserve heuristic helps avoid tons of reallocations.
    // Not “optimization-critical” yet, just a sanity default.
    const int L = l1 + l2 + l3 + l4;
    K.reserve(static_cast<std::size_t>(10 + 20 * (L + 1)));

    // Overall prefactors (no g1/g2 here; those are encoded into exponents below)
    const double pre1 =
        sign_pow(l1 + l2) *
        factorial(l1) *
        factorial(l2);

    const double pre2 =
        factorial(l3) *
        factorial(l4);

    for (int i1 = 0; i1 <= l1 / 2; ++i1)
    for (int i2 = 0; i2 <= l2 / 2; ++i2)
    for (int o1 = 0; o1 <= l1 - 2 * i1; ++o1)
    for (int o2 = 0; o2 <= l2 - 2 * i2; ++o2)
    for (int r1 = 0; r1 <= (o1 + o2) / 2; ++r1)
    {
        const double t11 =
            sign_pow(o2 + r1) *
            factorial(o1 + o2) /
            (ipow(4.0, i1 + i2 + r1) *
             factorial(i1) *
             factorial(i2) *
             factorial(o1) *
             factorial(o2) *
             factorial(r1));

        const double d12 =
            factorial(l1 - 2 * i1 - o1) *
            factorial(l2 - 2 * i2 - o2) *
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
                (ipow(4.0, i3 + i4 + r2) *
                 factorial(i3) *
                 factorial(i4) *
                 factorial(o3) *
                 factorial(o4) *
                 factorial(r2));

            const double d22 =
                factorial(l3 - 2 * i3 - o3) *
                factorial(l4 - 2 * i4 - o4) *
                factorial(o3 + o4 - 2 * r2);

            const int mu =
                l1 + l2 + l3 + l4
                - 2 * (i1 + i2 + i3 + i4)
                - (o1 + o2 + o3 + o4);

            for (int u = 0; u <= mu / 2; ++u)
            {
                const double t3 =
                    sign_pow(u) *
                    factorial(mu) /
                    (ipow(4.0, u) *
                     factorial(u) *
                     factorial(mu - 2 * u));

                K.scalar.push_back(pre1 * pre2 * (t11 / d12) * (t21 / d22) * t3);

                // Exponent pattern:
                // a1,a2,a3,a4,g1,g2,AB,CD,eta,PQ
                K.powers.push_back({
                    static_cast<int16_t>(o2 - i1 - r1),
                    static_cast<int16_t>(o1 - i2 - r1),
                    static_cast<int16_t>(o4 - i3 - r2),
                    static_cast<int16_t>(o3 - i4 - r2),
                    static_cast<int16_t>(2 * (i1 + i2) + r1 - (l1 + l2)),
                    static_cast<int16_t>(2 * (i3 + i4) + r2 - (l3 + l4)),
                    static_cast<int16_t>(o1 + o2 - 2 * r1),
                    static_cast<int16_t>(o3 + o4 - 2 * r2),
                    static_cast<int16_t>(mu - u),
                    static_cast<int16_t>(mu - 2 * u)
                });

                K.mu.push_back(static_cast<uint16_t>(mu));
                K.u.push_back(static_cast<uint16_t>(u));
            }
        }
    }

    return K;
}