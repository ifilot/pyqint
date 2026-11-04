#include "hellsing_cache.h"

/**
 * @brief Compute the flat index into the kernel table
 *
 * Maps a 4-tuple of angular momenta (l1, l2, l3, l4) into a single
 * linear index assuming all values satisfy:
 *
 *   0 <= li <= lmax
 *
 * The layout corresponds to a row-major ordering over the four indices.
 *
 * @param l1   Angular momentum of the first GTO
 * @param l2   Angular momentum of the second GTO
 * @param l3   Angular momentum of the third GTO
 * @param l4   Angular momentum of the fourth GTO
 * @param lmax Maximum angular momentum supported by the table
 *
 * @return Flat index into the kernel table
 */
std::size_t HellsingCacheTable1D::index(int l1, int l2, int l3, int l4, int lmax) {
    const std::size_t n = static_cast<std::size_t>(lmax) + 1;

    return ((static_cast<std::size_t>(l1) * n
           + static_cast<std::size_t>(l2)) * n
           + static_cast<std::size_t>(l3)) * n
           + static_cast<std::size_t>(l4);
}

/**
 * @brief Construct and precompute all one-dimensional Hellsing kernels
 *
 * Builds a complete lookup table containing Hellsing kernels for all
 * combinations (l1, l2, l3, l4) such that:
 *
 *   0 <= li <= lmax
 *
 * All kernels are generated eagerly and stored in a flat array for
 * fast lookup during integral evaluation.
 *
 * @param lmax Maximum angular momentum supported by the cache
 */
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

/**
 * @brief Retrieve a precomputed one-dimensional Hellsing kernel
 *
 * This function performs a pure lookup into the kernel table.
 * The caller must ensure that the cache has already been expanded
 * to support the requested angular momenta.
 *
 * @param l1 Angular momentum of the first GTO
 * @param l2 Angular momentum of the second GTO
 * @param l3 Angular momentum of the third GTO
 * @param l4 Angular momentum of the fourth GTO
 *
 * @return Reference to the corresponding Hellsing kernel
 */
const HellsingCache1D& HellsingCacheTable1D::get(int l1, int l2, int l3, int l4) const {
    return table_[index(l1, l2, l3, l4, lmax_)];
}

/**
 * @brief Ensure the cache supports at least the requested angular momentum
 *
 * If the current cache already supports the required angular momentum,
 * this function does nothing. Otherwise, the cache is rebuilt from scratch
 * with the larger limit.
 *
 * This function must be called outside any parallel region.
 *
 * @param required_lmax Required maximum angular momentum
 */
void HellsingCacheTable1D::ensure_lmax(int required_lmax) {
    if (required_lmax <= lmax_) {
        return;
    }

    // Build a completely new cache with the larger angular momentum
    HellsingCacheTable1D new_cache(required_lmax);

    // Atomically replace the current cache
    *this = std::move(new_cache);
}

/**
 * @brief Build a one-dimensional Hellsing kernel for a specific quartet
 *
 * Constructs all combinatorial terms required to evaluate the Hellsing
 * expansion for a given set of angular momenta (l1, l2, l3, l4).
 *
 * The resulting kernel contains:
 *  - scalar prefactors
 *  - integer exponent patterns
 *  - Boys-function bookkeeping indices
 *
 * @param l1 Angular momentum of the first GTO
 * @param l2 Angular momentum of the second GTO
 * @param l3 Angular momentum of the third GTO
 * @param l4 Angular momentum of the fourth GTO
 *
 * @return Fully constructed one-dimensional Hellsing kernel
 */
HellsingCache1D HellsingCacheTable1D::build_kernel_1d(int l1, int l2, int l3, int l4) {
    HellsingCache1D K;

    // Rough reserve heuristic to reduce reallocations
    const int L = l1 + l2 + l3 + l4;
    K.reserve(static_cast<std::size_t>(10 + 20 * (L + 1)));

    // Overall prefactors (Gaussian exponents encoded later via powers)
    const double pre1 =
        sign_pow(l1 + l2) *
        factorial(l1) *
        factorial(l2);

    const double pre2 =
        factorial(l3) *
        factorial(l4);

    // First pair (A,B)
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

        // Second pair (C,D)
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

            // Boys-function summation
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
                std::array<int16_t, 10> p = {
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
                };

                K.powers.push_back(p);
                K.update_bounds(p);

                K.mu.push_back(static_cast<uint16_t>(mu));
                K.u.push_back(static_cast<uint16_t>(u));
            }
        }
    }

    return K;
}
