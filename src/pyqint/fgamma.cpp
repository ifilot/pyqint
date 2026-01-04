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

// -------------------- Machine constants --------------------
inline constexpr double kEps = std::numeric_limits<double>::epsilon();

// -------------------- Small-T series --------------------
inline double boys_series(int n, double T) {
    double term = 1.0 / (2.0 * n + 1.0);
    double sum  = term;

    for (int k = 0; k < 2000; ++k) {
        const double num = (2.0 * n + 2.0 * k + 1.0);
        const double den = (2.0 * n + 2.0 * k + 3.0);
        term *= (-T) / double(k + 1) * (num / den);
        const double newsum = sum + term;

        if (std::abs(term) <= std::abs(newsum) * (50.0 * kEps))
            return newsum;

        sum = newsum;
    }
    return sum;
}

// -------------------- Accurate F0 --------------------
inline double boys_F0(double T) {
    if (T == 0.0) return 1.0;

    if (T < 1e-10) {
        const double T2 = T * T;
        return 1.0 - T / 3.0 + T2 / 10.0 - (T2 * T) / 42.0;
    }

    const double r = std::sqrt(T);
    return 0.5 * std::sqrt(M_PI / T) * std::erf(r);
}

// ==================== Construction ====================

BoysFunction::BoysFunction(int nu_max)
    : nu_tab_max_(-1),
      logT_min_(0.0),
      logT_max_(0.0),
      inv_dlogT_(0.0) {
    init_table(nu_max);
}

// ==================== Public API ====================

double BoysFunction::operator()(int nu, double T) const
{
    if (T < 0.0)
        return std::numeric_limits<double>::quiet_NaN();

    if (T < T_SMALL || nu > nu_tab_max_)
        return Fgamma_exact(nu, T);

    if (T <= T_LARGE)
        return interpolate(nu, T);

    return Fgamma_exact(nu, T);
}

void BoysFunction::compute_block(int nu_max, double T, double* F) const
{
    if (T < T_SMALL) {
        Fgamma_block_exact(nu_max, T, F);
        return;
    }

    const double expmT = std::exp(-T);

    if (T > 2.0 * nu_max) {
        // -------- Upward recurrence --------
        F[0] = (*this)(0, T);
        const double inv2T = 1.0 / (2.0 * T);

        for (int nu = 0; nu < nu_max; ++nu)
            F[nu + 1] =
                ((2.0 * nu + 1.0) * F[nu] - expmT) * inv2T;
    }
    else {
        // -------- Downward recurrence --------
        F[nu_max] = (*this)(nu_max, T);

        for (int nu = nu_max; nu > 0; --nu)
            F[nu - 1] =
                (2.0 * T * F[nu] + expmT) / (2.0 * nu - 1.0);
    }
}

int BoysFunction::max_order() const noexcept
{
    return nu_tab_max_;
}

// ==================== Table initialization ====================

void BoysFunction::init_table(int nu_max) {
    nu_tab_max_ = std::max(0, nu_max);
    Htab_.resize((nu_tab_max_ + 1) * NTABLE);

    logT_min_  = std::log(T_SMALL);
    logT_max_  = std::log(T_LARGE);
    inv_dlogT_ = (NTABLE - 1) / (logT_max_ - logT_min_);

    for (int i = 0; i < NTABLE; ++i) {
        const double logT =
            logT_min_ + (logT_max_ - logT_min_) * i / (NTABLE - 1);

        const double T     = std::exp(logT);
        const double sqrtT = std::sqrt(T);

        double Tpow = sqrtT;

        for (int nu = 0; nu <= nu_tab_max_; ++nu) {
            const double Fnu = Fgamma_exact(nu, T);
            Htab_[idx(nu, i)] = Tpow * Fnu;
            Tpow *= T;
        }
    }
}

// ==================== Interpolation ====================

double BoysFunction::interpolate(int nu, double T) const noexcept {
    const double x = (std::log(T) - logT_min_) * inv_dlogT_;

    int i = static_cast<int>(x);
    if (i < 0)           i = 0;
    if (i > NTABLE - 2)  i = NTABLE - 2;

    const double w  = x - i;
    const double h0 = Htab_[idx(nu, i)];
    const double h1 = Htab_[idx(nu, i + 1)];

    const double H = (1.0 - w) * h0 + w * h1;

    const double sqrtT = std::sqrt(T);
    double Tnu = 1.0;
    for (int k = 0; k < nu; ++k)
        Tnu *= T;

    return H / (sqrtT * Tnu);
}

int BoysFunction::idx(int nu, int i) const noexcept {
    return nu * NTABLE + i;
}

// ==================== Exact reference ====================

double BoysFunction::Fgamma_exact(int nu, double T) {
    if (T < 0.1)
        return boys_series(nu, T);

    // erf-based + upward recurrence
    const double expmT = std::exp(-T);
    double F = boys_F0(T);

    const double inv2T = 0.5 / T;
    for (int k = 0; k < nu; ++k)
        F = ((2.0 * k + 1.0) * F - expmT) * inv2T;

    return F;
}

void BoysFunction::Fgamma_block_exact(int nu_max, double T, double* F) {
    if (T < 0.1) {
        for (int nu = 0; nu <= nu_max; ++nu)
            F[nu] = boys_series(nu, T);
        return;
    }

    const double expmT = std::exp(-T);
    F[0] = boys_F0(T);

    const double inv2T = 0.5 / T;
    for (int nu = 0; nu < nu_max; ++nu)
        F[nu + 1] =
            ((2.0 * nu + 1.0) * F[nu] - expmT) * inv2T;
}