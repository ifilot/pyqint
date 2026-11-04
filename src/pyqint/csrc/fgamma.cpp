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

#include "fgamma.h"

// -------------------- Series expansion --------------------

/**
 * @brief Evaluate F_n(T) by its series of positive terms,
 *
 *   F_n(T) = exp(-T) Σ_k (2T)^k / ((2n+1)(2n+3)…(2n+2k+1)),
 *
 * which converges for all T without cancellation. It needs about T + 40
 * terms, so it is used for moderate T only.
 *
 * @param n  Boys order
 * @param T  argument
 *
 * @return F_n(T)
 */
inline double boys_series(int n, double T) {
    double term = 1.0 / (2.0 * n + 1.0);
    double sum  = term;
    for (int k = 1; k < 10000; ++k) {
        term *= 2.0 * T / (2.0 * n + 2.0 * k + 1.0);
        sum += term;
        if (term <= sum * std::numeric_limits<double>::epsilon())
            break;
    }
    return std::exp(-T) * sum;
}

// -------------------- Large-T evaluation --------------------

/**
 * @brief Evaluate F_0(T) … F_nu(T) from F_0(T) = sqrt(π/T) erf(sqrt T) / 2 by
 *        upward recurrence. The recurrence is stable only for T well above nu.
 *
 * @param nu  maximum order
 * @param T   argument (T > 0)
 * @param F   output array (size ≥ nu+1)
 */
inline void boys_upward(int nu, double T, double* F) {
    const double expmT = std::exp(-T);
    const double inv2T = 0.5 / T;
    F[0] = 0.5 * std::sqrt(M_PI / T) * std::erf(std::sqrt(T));
    for (int k = 0; k < nu; ++k)
        F[k + 1] = ((2.0 * k + 1.0) * F[k] - expmT) * inv2T;
}

// ==================== Construction ====================

/**
 * @brief Construct a Boys function evaluator.
 *
 * @param nu_max  Maximum Boys order to precompute.
 */
BoysFunction::BoysFunction(int nu_max)
    : nu_tab_max_(-1),
      nrows_(0),
      npoints_(0) {
    init_table(nu_max);
}

// ==================== Public API ====================

/**
 * @brief Evaluate the Boys function F_ν(T).
 *
 * @param nu  Boys order
 * @param T   argument
 *
 * @return Boys function value
 */
double BoysFunction::operator()(int nu, double T) const
{
    if (T < 0.0)
        return std::numeric_limits<double>::quiet_NaN();
    if (T < T_TABLE && nu <= nu_tab_max_)
        return taylor(nu, T);
    return Fgamma_exact(nu, T);
}

/**
 * @brief Compute F_0(T) … F_νmax(T).
 *
 * @param nu_max  Maximum order
 * @param T       argument
 * @param F       output array (size ≥ nu_max+1)
 */
void BoysFunction::compute_block(int nu_max, double T, double* F) const
{
    if (T >= T_TABLE && T > 2.0 * nu_max) {
        // -------- Upward recurrence (stable for T >> nu) --------
        boys_upward(nu_max, T, F);
        return;
    }

    // -------- Downward recurrence (stable) --------
    F[nu_max] = (*this)(nu_max, T);

    const double expmT = std::exp(-T);
    for (int nu = nu_max; nu > 0; --nu)
        F[nu - 1] =
            (2.0 * T * F[nu] + expmT) / (2.0 * nu - 1.0);
}

/**
 * @brief Return the maximum order stored in the table.
 *
 * @return maximum Boys order
 */
int BoysFunction::max_order() const noexcept
{
    return nu_tab_max_;
}

// ==================== Table initialization ====================

/**
 * @brief Initialize the table of F_ν(T_i) for the orders needed by the
 *        Taylor expansion of F_0 … F_νmax.
 *
 * @param nu_max  maximum order
 */
void BoysFunction::init_table(int nu_max) {
    nu_tab_max_ = std::max(0, nu_max);
    nrows_      = nu_tab_max_ + NTAYLOR;
    npoints_    = static_cast<int>(T_TABLE / DT_TABLE) + 2;
    Ftab_.resize(static_cast<size_t>(npoints_) * nrows_);

    for (int i = 0; i < npoints_; ++i) {
        const double T = i * DT_TABLE;
        double* row = &Ftab_[static_cast<size_t>(i) * nrows_];

        // highest order from the series, lower orders by downward recurrence
        row[nrows_ - 1] = boys_series(nrows_ - 1, T);
        const double expmT = std::exp(-T);
        for (int nu = nrows_ - 1; nu > 0; --nu)
            row[nu - 1] = (2.0 * T * row[nu] + expmT) / (2.0 * nu - 1.0);
    }
}

// ==================== Taylor expansion ====================

/**
 * @brief Evaluate F_ν(T) by a Taylor expansion around the nearest
 *        tabulated point T_i, using dF_ν/dT = -F_ν+1:
 *
 *   F_ν(T) = Σ_k F_ν+k(T_i) (T_i - T)^k / k!
 *
 * With |T - T_i| ≤ DT_TABLE/2, the truncation error is below
 * (DT_TABLE/2)^NTAYLOR / NTAYLOR! ≈ 1.5e-17 relative to F_ν.
 *
 * @param nu  Boys order
 * @param T   argument
 *
 * @return Boys function value
 */
double BoysFunction::taylor(int nu, double T) const noexcept {
    const int i = static_cast<int>(T / DT_TABLE + 0.5);
    const double d = i * DT_TABLE - T;
    const double* row = &Ftab_[static_cast<size_t>(i) * nrows_ + nu];

    // Horner scheme for Σ_k row[k] d^k / k!
    static constexpr double inv[NTAYLOR] = {0.0, 1.0, 1.0 / 2, 1.0 / 3, 1.0 / 4,
                                            1.0 / 5, 1.0 / 6, 1.0 / 7};
    static_assert(NTAYLOR == 8, "inv[] holds 1/k for k < NTAYLOR");
    double sum = row[NTAYLOR - 1];
    for (int k = NTAYLOR - 1; k > 0; --k)
        sum = row[k - 1] + sum * (d * inv[k]);

    return sum;
}

// ==================== Exact reference ====================

/**
 * @brief Compute the Boys function without the table: by the positive
 *        series for moderate T, by upward recurrence where that is stable.
 *
 * @param nu  Boys order
 * @param T   argument
 *
 * @return exact value
 */
double BoysFunction::Fgamma_exact(int nu, double T) {
    if (T < 2.0 * nu + 40.0)
        return boys_series(nu, T);

    std::vector<double> F(nu + 1);
    boys_upward(nu, T, F.data());
    return F[nu];
}
