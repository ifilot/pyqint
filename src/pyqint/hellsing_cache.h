#pragma once

#include <array>
#include <cstdint>
#include <vector>
#include <cstddef>
#include <limits>

#include "mathfuncs.h"
#include "factorials.h"

/**
 * @brief One-dimensional Hellsing kernel cache
 *
 * This structure stores all combinatorial terms required to evaluate
 * one Cartesian component (x, y, or z) of the electron repulsion integral
 * using the Hellsing expansion.
 *
 * Each entry corresponds to a single term in the expansion and consists of:
 *  - a scalar prefactor
 *  - a set of integer exponents describing powers of geometric and
 *    Gaussian-related quantities
 *  - Boys-function bookkeeping indices
 *
 * The cache is designed to be evaluated efficiently by precomputing
 * exponent bounds for tight power tables.
 */
struct HellsingCache1D {
    /**
     * @brief Scalar combinatorial prefactor for each expansion term
     *
     * These values include factorials, signs, and normalization factors
     * arising from the Hellsing recurrence relations.
     */
    std::vector<double> scalar;

    /**
     * @brief Integer exponent pattern for each expansion term
     *
     * The entries correspond to the following quantities:
     *  - [0] a1   : power of alpha_a
     *  - [1] a2   : power of alpha_b
     *  - [2] a3   : power of alpha_c
     *  - [3] a4   : power of alpha_d
     *  - [4] g1   : power of gamma_1
     *  - [5] g2   : power of gamma_2
     *  - [6] AB   : power of (A - B)
     *  - [7] CD   : power of (C - D)
     *  - [8] eta  : power of eta
     *  - [9] PQ   : power of (P - Q)
     *
     * All powers are stored as signed 16-bit integers.
     */
    std::vector<std::array<int16_t, 10>> powers;

    /**
     * @brief Boys-function order parameter (mu)
     *
     * This is the total order of the Boys function required for the term.
     */
    std::vector<uint16_t> mu;

    /**
     * @brief Boys-function summation index (u)
     *
     * Each term contributes to Boys function values of order (mu - u).
     */
    std::vector<uint16_t> u;

    /**
     * @brief Minimum exponent encountered for each power component
     *
     * Used to construct tight power tables and avoid unnecessary evaluations.
     */
    std::array<int16_t, 10> pmin;

    /**
     * @brief Maximum exponent encountered for each power component
     *
     * Used to construct tight power tables and avoid unnecessary evaluations.
     */
    std::array<int16_t, 10> pmax;

    /**
     * @brief Default constructor
     *
     * Initializes exponent bounds to extreme values so they can be
     * tightened during kernel construction.
     */
    HellsingCache1D() {
        pmin.fill(std::numeric_limits<int16_t>::max());
        pmax.fill(std::numeric_limits<int16_t>::min());
    }

    /**
     * @brief Number of expansion terms stored
     *
     * @return Number of cached Hellsing expansion terms
     */
    std::size_t size() const noexcept { return scalar.size(); }

    /**
     * @brief Clear all cached data
     *
     * Resets all vectors and exponent bounds.
     */
    void clear() {
        scalar.clear();
        powers.clear();
        mu.clear();
        u.clear();
        pmin.fill(std::numeric_limits<int16_t>::max());
        pmax.fill(std::numeric_limits<int16_t>::min());
    }

    /**
     * @brief Reserve memory for expansion terms
     *
     * @param n Number of terms to reserve space for
     */
    void reserve(std::size_t n) {
        scalar.reserve(n);
        powers.reserve(n);
        mu.reserve(n);
        u.reserve(n);
    }

    /**
     * @brief Update exponent bounds for a newly added power pattern
     *
     * @param p Exponent array corresponding to a single expansion term
     */
    inline void update_bounds(const std::array<int16_t,10>& p) noexcept {
        for (int k = 0; k < 10; ++k) {
            if (p[k] < pmin[k]) pmin[k] = p[k];
            if (p[k] > pmax[k]) pmax[k] = p[k];
        }
    }
};

/**
 * @brief Flat lookup table for one-dimensional Hellsing kernels
 *
 * This class stores precomputed HellsingCache1D kernels for all combinations
 * of angular momenta (l1, l2, l3, l4) such that:
 *
 *   0 <= li <= lmax
 *
 * The kernels are indexed in a flat array for fast lookup during
 * electron repulsion integral evaluation.
 */
class HellsingCacheTable1D {
public:
    /**
     * @brief Construct and precompute all kernels up to a given angular momentum
     *
     * @param lmax Maximum angular momentum supported by the cache
     */
    explicit HellsingCacheTable1D(int lmax);

    /**
     * @brief Ensure the cache supports at least the requested angular momentum
     *
     * If the current cache is too small, it is rebuilt to accommodate
     * the larger value.
     *
     * @param required_lmax Required maximum angular momentum
     */
    void ensure_lmax(int required_lmax);

    /**
     * @brief Return the maximum angular momentum currently supported
     *
     * @return Maximum angular momentum supported by the cache
     */
    int lmax() const noexcept { return lmax_; }

    /**
     * @brief Retrieve a precomputed Hellsing kernel
     *
     * @param l1 Angular momentum of the first GTO
     * @param l2 Angular momentum of the second GTO
     * @param l3 Angular momentum of the third GTO
     * @param l4 Angular momentum of the fourth GTO
     *
     * @return Reference to the corresponding one-dimensional Hellsing kernel
     */
    const HellsingCache1D& get(int l1, int l2, int l3, int l4) const;

private:
    /**
     * @brief Compute flat array index for a kernel
     *
     * @param l1 Angular momentum of the first GTO
     * @param l2 Angular momentum of the second GTO
     * @param l3 Angular momentum of the third GTO
     * @param l4 Angular momentum of the fourth GTO
     * @param lmax Maximum angular momentum supported
     *
     * @return Flat index into the kernel table
     */
    static std::size_t index(int l1, int l2, int l3, int l4, int lmax);

    /**
     * @brief Build a one-dimensional Hellsing kernel for a specific quartet
     *
     * @param l1 Angular momentum of the first GTO
     * @param l2 Angular momentum of the second GTO
     * @param l3 Angular momentum of the third GTO
     * @param l4 Angular momentum of the fourth GTO
     *
     * @return Fully constructed one-dimensional Hellsing kernel
     */
    HellsingCache1D build_kernel_1d(int l1, int l2, int l3, int l4);

    /// Maximum angular momentum supported by the table
    int lmax_{0};

    /// Flat storage of all precomputed kernels
    std::vector<HellsingCache1D> table_;
};