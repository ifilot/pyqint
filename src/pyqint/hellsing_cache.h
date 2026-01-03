#pragma once

#include <array>
#include <cstdint>
#include <vector>
#include <cstddef>

#include "mathfuncs.h"
#include "factorials.h"

#include <array>
#include <cstdint>
#include <limits>
#include <vector>

struct HellsingCache1D {
    // Combinatorial scalar prefactor per term
    std::vector<double> scalar;

    // Integer exponent pattern
    // [0] a1, [1] a2, [2] a3, [3] a4,
    // [4] g1, [5] g2, [6] AB, [7] CD,
    // [8] eta, [9] PQ
    std::vector<std::array<int16_t, 10>> powers;

    // Boys bookkeeping
    std::vector<uint16_t> mu;
    std::vector<uint16_t> u;

    // NEW: exponent bounds for tight power tables
    std::array<int16_t, 10> pmin;
    std::array<int16_t, 10> pmax;

    HellsingCache1D() {
        pmin.fill(std::numeric_limits<int16_t>::max());
        pmax.fill(std::numeric_limits<int16_t>::min());
    }

    std::size_t size() const noexcept { return scalar.size(); }

    void clear() {
        scalar.clear();
        powers.clear();
        mu.clear();
        u.clear();
        pmin.fill(std::numeric_limits<int16_t>::max());
        pmax.fill(std::numeric_limits<int16_t>::min());
    }

    void reserve(std::size_t n) {
        scalar.reserve(n);
        powers.reserve(n);
        mu.reserve(n);
        u.reserve(n);
    }

    inline void update_bounds(const std::array<int16_t,10>& p) noexcept {
        for (int k = 0; k < 10; ++k) {
            if (p[k] < pmin[k]) pmin[k] = p[k];
            if (p[k] > pmax[k]) pmax[k] = p[k];
        }
    }
};

// Flat table of kernels for all 0 <= li <= lmax.
class HellsingCacheTable1D {
public:
    explicit HellsingCacheTable1D(int lmax);

    int lmax() const noexcept { return lmax_; }

    const HellsingCache1D& get(int l1, int l2, int l3, int l4) const;

private:
    static std::size_t index(int l1, int l2, int l3, int l4, int lmax);

    HellsingCache1D build_kernel_1d(int l1, int l2, int l3, int l4);

    int lmax_{0};
    std::vector<HellsingCache1D> table_;
};